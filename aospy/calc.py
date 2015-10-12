"""calc.py: classes for performing specified calculations on aospy data"""
import os
import pickle
import shutil
import subprocess
import tarfile
import time

import numpy as np
import pandas as pd
import xray

from . import Constant, Var
from .io import (_data_in_label, _data_out_label, _ens_label,
                 _month_indices, _yr_label, dmget, nc_name_gfdl,
                 _get_time_xray)
from .utils import (get_parent_attr, level_thickness_xray,
                    pfull_from_sigma_xray, dp_from_sigma_xray, int_dp_g,
                    int_dp_g_xray)

ps = Var(
    name='ps',
    units='Pa',
    domain='atmos',
    description='Surface pressure.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)


class CalcInterface(object):
    """Interface to Calc class."""
    def _set_nc_attrs(self):
        for attr in ('nc_start_date',
                     'nc_end_date',
                     'nc_dur',
                     'direc_nc',
                     'default_time_range',
                     'ens_mem_prefix',
                     'ens_mem_ext',
                     'nc_files',
                     'nc_dir_struc',
                     'default_yr_range',
                     'idealized'):

            attr_val = tuple([get_parent_attr(rn, attr, strict=False) for rn
                              in self.run])
            setattr(self, attr, attr_val)

    def _get_yr_range(self):
        """Set the object's span of years."""
        if self.yr_range == 'default':
            start_date, end_date = get_parent_attr(self.run[0],
                                                   'default_time_range')
        elif self.yr_range == 'all':
            start_date = get_parent_attr(self.run[0], 'nc_start_date')
            end_date = get_parent_attr(self.run[0], 'nc_end_date')
        else:
            start_date, end_date = self.yr_range
        return start_date, end_date

    def _get_num_yr(self):
        """Compute effective number of years in the input data."""
        if self.dtype_in_time in ('ts', 'inst'):
            num_yr = (self.end_yr['files'].year -
                      self.start_yr['files'].year + 1)
        else:
            num_yr = 1
        return num_yr

    def _make_time_chunks(self):
        """Create tuple of (start, end) pairs based on given year chunks."""
        if self.yr_chunk_len:
            # For now we will assume there are no year chunks. It won't be too hard to
            # extend this, but we will need to make parallel ranges (one in the xray set
            # of years and one in the files set of years). Then we will return a list of
            # dictionary tuples.
            dur = self.yr_chunk_len - 1
            st_yrs = list(pd.date_range(self.start_analysis_year,
                                   self.end_analysis_year,
                                   freq='%dAS-JAN' % self.yr_chunk_len))
            offset = self.yr_chunk_len*pd.tseries.offsets.YearEnd() + pd.tseries.offsets.DateOffset(days=1)
            end_yrs = list(pd.date_range(self.start_analysis_year + offset,
                                         self.end_analysis_year + offset,
                                         freq='%dAS-JAN' % self.yr_chunk_len))
        else:
            st_yrs, end_yrs = [self.start_yr], [self.end_yr]
        return zip(st_yrs, end_yrs)

    def __init__(self, proj=None, model=None, run=None, ens_mem=None, var=None,
                 yr_range=None, region=None, intvl_in=None, intvl_out=None,
                 dtype_in_time=None, dtype_in_vert=None, dtype_out_time=None,
                 dtype_out_vert=None, level=None, skip_time_inds=False,
                 yr_chunk_len=False, verbose=True):
        """Create the CalcInterface object with the given parameters."""
        if run not in model.runs.values():
            raise AttributeError("Model '%s' has no run '%s'.  Calc object "
                                 "will not be generated." % (model, run))
        # Everything gets turned into tuples.
        proj = tuple([proj])
        model = tuple([model])
        if not isinstance(run, (list, tuple)):
            run = tuple([run])
        # Make tuples the same length.
        if len(proj) == 1 and (len(model) > 1 or len(run) > 1):
            proj = tuple(list(proj)*len(run))
        if len(model) == 1 and len(run) > 1:
            model = tuple(list(model)*len(run))

        self.proj = proj
        self.model = model
        self.run = run

        self._set_nc_attrs()

        self.proj_str = '_'.join(set([p.name for p in self.proj]))
        self.model_str = '_'.join(set([m.name for m in self.model]))
        run_names = [r.name for r in self.run]
        self.run_str = '_'.join(set(run_names))
        self.run_str_full = '_'.join(run_names)

        self.var = var
        self.name = self.var.name
        self.domain = self.var.domain
        self.def_time = self.var.def_time
        self.def_vert = self.var.def_vert
        self.verbose = verbose

        try:
            self.function = self.var.func
        except AttributeError:
            self.function = lambda x: x
        if getattr(self.var, 'variables', False):
            self.variables = self.var.variables
        else:
            self.variables = (self.var,)

        self.ens_mem = ens_mem
        self.level = level
        self.intvl_in = intvl_in
        self.intvl_out = intvl_out
        self.dtype_in_time = dtype_in_time
        self.dtype_in_vert = dtype_in_vert
        self.ps = ps
        print(dtype_out_time)
        if isinstance(dtype_out_time, (list, tuple)):
            self.dtype_out_time = tuple(dtype_out_time)
        else:
            self.dtype_out_time = tuple([dtype_out_time])
        self.dtype_out_vert = dtype_out_vert
        self.region = region

        self.yr_range = yr_range
        self.start_yr, self.end_yr = self._get_yr_range()

        # If we use the xray mode then all dates will be a dictionary.
        # One will map to the date in the filename, the other will map
        # to the date used in analysis.
        if self.start_yr.year < 1678:
            offset = 1900
        else:
            offset = 0
        self.start_yr = {
            'files' : self.start_yr,
            'xray' : pd.to_datetime(np.datetime64(
                    '%04d-%02d-%02d' % (self.start_yr.year + offset,
                                        self.start_yr.month,
                                        self.start_yr.day)))
            }
        self.end_yr = {
            'files' : self.end_yr,
            'xray' :  pd.to_datetime(np.datetime64(
                    '%04d-%02d-%02d' % (self.end_yr.year + offset,
                                        self.end_yr.month,
                                        self.end_yr.day)))
            }
        self.num_yr = self._get_num_yr()
        self.months = _month_indices(intvl_out, iterable=True)
        self.skip_time_inds = skip_time_inds
        self.yr_chunk_len = yr_chunk_len
        self.chunk_ranges = self._make_time_chunks()


class Calc(object):
    """Class for executing, saving, and loading a single computation."""
    def __init__(self, calc_interface):
        self.__dict__ = vars(calc_interface)
        self._print_verbose("\nInitializing Calc instance: %s", self.__str__())

        [mod.set_grid_data() for mod in self.model]

        if isinstance(calc_interface.ens_mem, int):
            self.direc_nc = self.direc_nc[calc_interface.ens_mem]

        self.dt_set = False

        self.dir_scratch = self._dir_scratch()
        self.dir_archive = self._dir_archive()
        self.file_name = {d: self._file_name(d) for d in self.dtype_out_time}
        self.path_scratch = {
            d: '/'.join([self.dir_scratch,
                         self.file_name[d]]).replace('//', '/')
            for d in self.dtype_out_time
        }
        self.path_archive = (
            '/'.join([self.dir_archive, 'data.tar']).replace('//', '/')
        )

        self.data_out = {}

    def __str__(self):
        """String representation of the object."""
        return "Calc object: " + ', '.join(
            (self.name, self.proj_str, self.model_str, self.run_str_full)
        )

    __repr__ = __str__

    def _print_verbose(self, *args):
        """Print diagnostic message."""
        if not self.verbose:
            pass
        else:
            try:
                print(args[0] % args[1], '(%s)' % time.ctime())
            except IndexError:
                print(args[0], '(%s)' % time.ctime())

    def _set_time_dt(self):
        """Get time and dt arrays at needed time indices."""
        # Use the first var in the list that is an aospy.Var object.
        nc_var = self.variables[0]
        for var in self.variables:
            if isinstance(var, Var) and not var.in_nc_grid:
                nc_var = var
                break
        with self._get_nc(nc_var, self.start_yr, self.end_yr) as nc:
            time_obj = nc['time']
            inds, time = _get_time_xray(
                time_obj, self.start_yr['xray'], self.end_yr['xray'],
                self.months, indices=True
                )
            self.time = time
            self.time_inds = inds
            self.dt = self._get_dt(nc, inds)

            for name in ('level', 'lev', 'plev'):
                try:
                    self.pressure = nc.variables[name][:]
                except KeyError:
                    pass
                else:
                    break
            else:
                self.pressure = False

            for name in ('lat', 'latitude', 'LATITUDE', 'y', 'yto'):
                try:
                    self.lat = nc.variables[name][:]
                except KeyError:
                    pass
                else:
                    break
            else:
                self.lat = False

            for name in ('lon', 'longitude', 'LONGITUDE', 'x', 'xto'):
                try:
                    self.lon = nc.variables[name][:]
                except KeyError:
                    pass
                else:
                    break
            else:
                self.lon = False

    def _get_dt(self, nc, indices):
        """Get durations of the desired timesteps."""
        if self.dtype_in_time == 'inst':
            return np.ones(np.shape(indices))
        for dt_name in ('average_DT',):
            try:
                dt = nc[dt_name]
            except:
                pass
            else:
                return dt.sel(time=indices)
        for name in ('time_bounds', 'time_bnds'):
            try:
                time_bounds = nc.variables[name]
            except KeyError:
                pass
            else:
                assert time_bounds.ndim == 2
                assert time_bounds.dimensions[1] == 'bnds'
                dt = time_bounds[:,1] - time_bounds[:,0]
                return dt[indices]
        raise ValueError("dt array could not be created.")

    def _reshape_time_indices(self, array, start_yr, end_yr):
        """Reshape time array to have year- and within-year axes.

        2015-04-23: This might not work for sub-monthly data spanning leap
        years or other calendar idiosyncracies.  For example, suppose using
        3 hourly data for DJF spanning a leap year and non leap-year.  The
        leap year February will have 4 more timesteps than the non-leap year.
        """

        reshaped = np.reshape(array, (end_yr - start_yr + 1, -1))
        return reshaped[:,:,np.newaxis,np.newaxis,np.newaxis]

    def _dir_scratch(self):
        """Create string of the data directory on the scratch filesystem."""
        ens_label = _ens_label(self.ens_mem)
        dir_scratch = '/'.join(
            ['/work', os.getenv('USER'), self.proj_str, self.model_str,
             self.run_str, ens_label, self.var.name]
        ).replace('//', '/')
        return dir_scratch

    def _dir_archive(self):
        """Create string of the data directory on the archive filesystem."""
        ens_label = _ens_label(self.ens_mem)
        dir_archive = '/'.join(
            ['/archive', os.getenv('USER'), self.proj_str, 'data',
             self.model_str, self.run_str, ens_label]
        ).replace('//', '/')
        return dir_archive

    def _file_name(self, dtype_out_time):
        """Create the name of the aospy file."""
        out_lbl = _data_out_label(self.intvl_out, dtype_out_time,
                                  dtype_vert=self.dtype_out_vert)
        in_lbl = _data_in_label(self.intvl_in, self.dtype_in_time,
                                self.dtype_in_vert)
        ens_lbl = _ens_label(self.ens_mem)
        yr_lbl = _yr_label((self.start_yr['files'].year, self.end_yr['files'].year))
        file_name = '.'.join(
            [self.name, out_lbl, in_lbl, self.model_str, self.run_str_full,
             ens_lbl, yr_lbl, 'p']
        ).replace('..', '.')
        return file_name

    def _get_nc_one_dir_tar(self, name, direc_nc, n=0):
        """Get the names of the tar files when all in the same directory."""

        # tar may hold absolute or relative paths
        paths = []
        for tar in nc_files:
            full = '/'.join([direc_nc, tar]).replace('//', '/')
            if os.path.isfile(tar):
                paths.append(tar)
            elif os.path.isfile(full):
                paths.append(full)
            else:
                print("Warning: specified netCDF file `%s` not found" % nc)
        # Remove duplicate entries.
        files = list(set(paths))
        files.sort()
        return files

    def _get_nc_one_dir(self, name, direc_nc, n=0):
        """Get the names of netCDF files when all in same directory."""
        if isinstance(self.nc_files[n][name], str):
            nc_files = [self.nc_files[n][name]]
        else:
            nc_files = self.nc_files[n][name]
        # nc_files may hold absolute or relative paths
        paths = []
        for nc in nc_files:
            full = '/'.join([direc_nc, nc]).replace('//', '/')
            if os.path.isfile(nc):
                paths.append(nc)
            elif os.path.isfile(full):
                paths.append(full)
            else:
                print("Warning: specified netCDF file `%s` not found" % nc)
        # Remove duplicate entries.
        files = list(set(paths))
        files.sort()
        return files

    def _get_nc_gfdl_repo(self, name, n=0):
        """Get the names of netCDF files from a GFDL repo on /archive."""
        return self.model[n].find_nc_direc_repo(
            run_name=self.run[n].name, var_name=name
        )

    def _get_nc_gfdl_dir_struct(self, name, direc_nc, start_yr, end_yr, n=0):
        """
        Get the names of netCDF files stored in GFDL standard directory
        structure and names.
        """
        domain = self.domain
        dtype_lbl = self.dtype_in_time
        if self.dtype_in_vert == 'sigma' and name != 'ps':
            domain += '_level'
        if self.dtype_in_time == 'inst':
            domain += '_inst'
            dtype_lbl = 'ts'
        # 2015-05-12 using 'ts' or 'inst' data only for now
        # separator = '_' if self.dtype_in_time == 'av' else '/'
        separator = '/'
        if 'av_from_' in self.dtype_in_time:
            dtype = self.dtype_in_time.replace('av_from_', '')
            dtype_lbl = dtype
        else:
            dtype = self.dtype_in_time
        direc = (direc_nc + '/' + domain + '/' + dtype_lbl + '/' +
                 self.intvl_in + separator + str(self.nc_dur[n]) + 'yr/')
        files = [direc + nc_name_gfdl(name, domain, dtype,
                                      self.intvl_in, yr, self.intvl_out,
                                      self.nc_start_yr[n], self.nc_dur[n])
                 for yr in range(start_yr, end_yr)]
        # Remove duplicate entries.
        files = list(set(files))
        files.sort()
        return files

    def _get_nc(self, var, start_yr=False, end_yr=False, n=0):
        """
        Create an MFDataset for the variable spanning all desired timesteps.

        Files chosen depend on the specified variables and time intvl and the
        attributes of the netCDF files.
        """
        if isinstance(self.direc_nc, str):
            direc_nc = self.direc_nc
        elif isinstance(self.direc_nc, (list, tuple)):
            direc_nc = self.direc_nc[n]
        else:
            raise IOError("direc_nc must be string, list, or tuple: %s"
                          % self.direc_nc)
        # Cycle through possible names until the data is found.
        names = [var.name]
        if hasattr(var, 'alt_names'):
            names += list(var.alt_names)
        for name in names:
            if self.nc_dir_struc[n] == 'one_dir':
                try:
                    files = self._get_nc_one_dir(name, direc_nc, n=n)
                except KeyError:
                    pass
                else:
                    pass # Need to figure out what exception is being thrown.
            elif self.nc_dir_struc[0].lower() == 'gfdl':
                files = self._get_nc_gfdl_dir_struct(
                    name, direc_nc, start_yr['files'].year,
                    end_yr['files'].year, n=n
                )
            elif self.nc_dir_struc[0].lower() == 'gfdl_repo':
                try:
                    files = self._get_nc_gfdl_repo(name, n=n)
                except IOError:
                    pass
                else:
                    break
            else:
                raise ValueError("Specified directory type not supported: %s"
                                 % self.nc_dir_struc[n])

            try:
                dmget(files)
                ds = []
                for file in files:
                    test = xray.open_dataset(file,
                                             decode_cf=False,
                                             drop_variables=['time_bounds','nv'])
                    if start_yr['files'].year < 1678:
                        for v in ['time', 'average_T1', 'average_T2']:
                            test[v].attrs['units'] = 'days since 1900-01-01 00:00:00'
                        test['time'].attrs['calendar'] = 'noleap'
                    test = xray.decode_cf(test)
                    ds.append(test)
                return xray.concat(ds, dim='time')
            except (RuntimeError, KeyError):
                pass
            else:
                raise IOError("Could not find files for variable '%s'." % var)
                try:
                    files = self._get_nc_gfdl_dir_struct(name, direc_nc,
                                                         start_yr, end_yr, n=n)
                except:
                    raise
                else:
                    break
        else:
            raise IOError("netCDF files for variable `%s`, year range %s-%s, "
                          "in directory %s, not found" % (var, start_yr,
                                                          end_yr, direc_nc))

    def _get_pressure_vals(self, var, start_yr, end_yr, n=0):
        """Get pressure array, whether sigma or standard levels."""
        self._print_verbose("Getting pressure data: %s", var)
        if self.dtype_in_vert == 'pressure':
            if np.any(self.pressure):
                pressure = self.pressure
            else:
                pressure = self.model[n].level
            if var == 'p':
                data = pressure
            elif var == 'dp':
                data = level_thickness_xray(pressure)

        if self.dtype_in_vert == 'sigma':
            bk = self.model[n].bk
            pk = self.model[n].pk
            ps = self._get_var_data(self.ps, start_yr, end_yr)
            pfull_coord = self.model[n].pfull
            if var == 'p':
                data = pfull_from_sigma_xray(bk, pk, ps, pfull_coord)
            elif var == 'dp':
                data = dp_from_sigma_xray(bk, pk, ps, pfull_coord)
        return data

    def _get_pressure_vals_xray(self, var, start_yr, end_yr, n=0):
        return

    def _get_data_subset(self, data, region=False, time=False,
                         vert=False, lat=False, lon=False, n=0):
        """Subset the data array to the specified time/level/lat/lon, etc."""
        if region:
            # if type(region) is str:
                # data = data[region]
            # elif type(region) is Region:
            data = data[region.name]
        if np.any(time):
            data = data[time]
            if 'av_from_' in self.dtype_in_time:
                data = np.mean(data, axis=0)[np.newaxis,:]
        if np.any(vert):
            if self.dtype_in_vert != 'sigma':
                if np.max(self.model[n].level) > 1e4:
                    # Convert from Pa to hPa.
                    lev_hpa = self.model[n].level*1e-2
                else:
                    lev_hpa = self.model[n].level
                level_index = np.where(lev_hpa == self.level)
                if 'ts' in self.dtype_out_time:
                    data = np.squeeze(data[:,level_index])
                else:
                    data = np.squeeze(data[level_index])
        if np.any(lat):
            pass
        if np.any(lon):
            pass
        return data

    def _get_var_data(self, var, start_yr, end_yr, n=0, region=False,
                      vert=False, lat=False, lon=False):
        """Get the needed data from one aospy.Var."""
        self._print_verbose("\tGetting data from netCDF files: %s", var)
        with self._get_nc(var, start_yr, end_yr, n=n) as nc:
            # Variable names can differ.
            try:
                data = nc[var.name]
            except KeyError:
                for alt_name in var.alt_names:
                    try:
                        data = nc[alt_name]
                    except KeyError:
                        pass
                    else:
                        break
            try:
                time = nc['time']
            except:
                pass
            t_inds = _get_time_xray(time, start_yr['xray'], end_yr['xray'],
                                    self.months, indices='only')
            data = self._get_data_subset(data, region=region, time=t_inds,
                                         vert=self.level, lat=lat, lon=lon)
        # Interpolate data at sigma half levels to full levels.
        if self.dtype_in_vert == 'sigma' and var.def_vert == 'phalf':
            data = False
            # We'll need to work on this. We just need to make sure the
            # coordinates in pressure align to be able to add the arrays.
            # It will just take some care.
        return data

    def _get_all_vars_data(self, start_yr, end_yr):
        """Get the needed data from all of the vars in the calculation."""
        all_vals = []
        for n, var in enumerate(self.variables):
            # If only 1 run, use it to load all data.
            # Otherwise assume that # runs == # vars to load.
            if len(self.run) == 1:
                n = 0
            # Pressure handled specially due to complications from sigma vs. p.
            if var in ('p', 'dp'):
                data = self._get_pressure_vals(var, start_yr, end_yr)
            # Pass numerical constants as is.
            elif isinstance(var, (float, int, Constant)):
                data = var
            # aospy.Var objects remain.
            # lat and lon maybe are Calc attributes; if not get from model
            elif var.name == 'lat':
                if np.any(self.lat):
                    data = self.lat
                else:
                    data = getattr(self.model[0], var.name)
            elif var.name == 'lon':
                if np.any(self.lon):
                    data = self.lon
                else:
                    data = getattr(self.model[0], var.name)
            # Grab the time array as is.
            elif var.name == 'time':
                data = self.time
            # Get other grid, time, etc. arrays directly from model object
            elif var.name in ('level', 'pk', 'bk', 'sfc_area'):
                data = getattr(self.model[0], var.name)
            else:
                data = self._get_var_data(var, start_yr, end_yr, n=n)
            all_vals.append(data)
        return all_vals

    def _local_ts(self, dp, dt, num_yr, *data_in):
        """Compute the function on the given gridded input data."""
        result = self.function(*data_in)
        print('result', result.shape)
        # Apply spatial reductions methods.
        if self.def_vert and self.dtype_out_vert == 'vert_int':
            result = int_dp_g(result, dp)
        # If already averaged, pass data on.  Otherwise do time averaging.
        if 'av' in self.dtype_in_time or not self.def_time:
            return result
        try:
            result = result.reshape((num_yr, -1, result.shape[-3],
                                     result.shape[-2], result.shape[-1]))
        except ValueError:
            result = result.reshape((num_yr, -1, result.shape[-2],
                                     result.shape[-1]))
        # Average within each year, yielding a yearly time-series.
        try:
            return np.ma.multiply(result, dt).sum(axis=1) / dt.sum(axis=1)
        except ValueError:
            dt = dt[:,:,:,:,0]
            return np.ma.multiply(result, dt).sum(axis=1) / dt.sum(axis=1)

    def _local_ts_xray(self, dp, dt, *data_in):
        result = self.function(*data_in)
        # Apply spatial reduction methods.
        if self.def_vert and self.dtype_out_vert == 'vert_int':
            result = int_dp_g_xray(result, dp)
        # If already averaged, pass data on. Otherwise do time averaging.
        if 'av' in self.dtype_in_time or not self.def_time:
            return result

        if self.idealized[0]:
            return result

        # Otherwise do time averaging over the years.
        result *= dt
        # Group by year.
        result = result.groupby('time.year')
        dt = dt.groupby('time.year')
        return result.sum('time') / dt.sum('time')

    def _time_reduce(self, loc_ts):
        """Compute all desired calculations on a local time-series."""
        files = {}
        if self.idealized[0]:
            tdim = 'time'
        else:
            tdim = 'year'
        if 'ts' in self.dtype_out_time:
            files.update({'ts': loc_ts})
        if 'None' in self.dtype_out_time:
            # Some calcs (e.g. correlations) already do time reduction.
            files.update({'av': loc_ts})
        if 'av' in self.dtype_out_time:
            files.update({'av': loc_ts.mean(tdim)})
        if 'eddy.av' in self.dtype_out_time:
            files.update({'eddy.av': loc_ts.mean(tdim)})
        if 'std' in self.dtype_out_time:
            files.update({'std': loc_ts.std(tdim)})
        if 'eddy.std' in self.dtype_out_time:
            files.update({'eddy.std': loc_ts.std(tdim)})
        # Zonal asymmetry.
        if any('zasym' in out_type for out_type in self.dtype_out_time):
            # '.T'=transpose; makes numpy broadcasting syntax work.
            znl_ts = loc_ts.mean('lon')
            zasym_ts = (loc_ts - znl_ts)
            if 'zasym.ts' in self.dtype_out_time:
                files.update({'zasym.ts': zasym_ts})
            if 'zasym.av' in self.dtype_out_time:
                files.update({'zasym.av': zasym_ts.mean(tdim)})
            if 'zasym.std' in self.dtype_out_time:
                files.update({'zasym.std': zasym_ts.std(tdim)})
        return files

    def region_calcs(self, loc_ts, n=0):
        """Region-averaged computations.  Execute and save to external file."""
        calcs_reg = ('ts_xray', 'av_xray', 'std_xray')
        # Perform each calculation for each region.
        for calc in calcs_reg:
            calc_name = ('reg.' + calc)
            if calc_name in self.dtype_out_time:
                reg_dat = {}
                for reg in self.region.values():
                    # Just pass along the data if averaged already.
                    if 'av' in self.dtype_in_time:
                        data_out = reg.ts_xray(loc_ts, self.model[n])
                    # Otherwise perform the calculation.
                    else:
                        method = getattr(reg, calc)
                        data_out = method(loc_ts, self.model[n])
                    reg_dat.update({reg.name: data_out})
                self.save(reg_dat, calc_name)

    def _compute_chunk(self, start_yr, end_yr):
        """Perform the calculation on the given chunk of times."""
        self._print_verbose("\nComputing desired timeseries from netCDF data "
                            "for years %d-%d.",
                            (start_yr['files'].year, end_yr['files'].year))
        if not self.dt_set:
            self._set_time_dt()
            self.dt_set = True

        inds = _get_time_xray(
            self.time, start_yr['xray'], end_yr['xray'], self.months, indices='only'
            )
        dt = self.dt.sel(time=inds).astype(float)
        # Reshape time-indices basically makes it easier to group by year.
        # We should be able to do that in xray style a bit more transparently.
        #print(dt.astype(float)*1e-9/86400.0)

        #dt = (self.dt.astype(float)*1e-9/86400.0).groupby('time.year')
        data_in = self._get_all_vars_data(start_yr, end_yr)
        if self.dtype_out_vert == 'vert_int':
            dp = self._get_pressure_vals('dp', start_yr['files'], end_yr['files'])
        else:
            dp = False
        return self._local_ts_xray(dp, dt, *data_in)

    def compute(self):
        """Perform all desired calculations on the data and save externally."""
        # Compute the local time series for each chunk and then combine chunks.
        full_ts = [self._compute_chunk(start_yr, end_yr)
                       for start_yr, end_yr in self.chunk_ranges]
        if len(full_ts) == 1:
            full_ts = full_ts[0]
        else:
            self._print_verbose("Combining output data from all chunks.")
            full_ts = np.ma.concatenate(full_ts, axis=0)
        self.full_ts = full_ts
        # Compute the time series using monthly means if needed.
        # if any(['eddy' in do for do in self.dtype_out_time] +
        #        ['av_from' in do for do in self.dtype_out_time]):
        #     monthly_ts = [self._compute_chunk(start_yr, end_yr, monthly=True)
        #                   for start_yr, end_yr in self.chunk_ranges]
        #     if len(monthly_ts) == 1:
        #         monthly_ts = monthly_ts[0]
        #     else:
        #         self._print_verbose("Combining output data from all chunks.")
        #         monthly_ts = np.ma.concatenate(monthly_ts, axis=0)
        #     self.monthly_ts = monthly_ts
        # Apply time reduction methods on gridded data and save.
        self._print_verbose("Applying desired time-reduction methods.")
        if self.def_time:
            reduced = self._time_reduce(full_ts)
        else:
            reduced = {'': full_ts}
        self._print_verbose("Writing desired gridded outputs to disk.")
        for dtype_out_time, data in reduced.items():
            self.save(np.squeeze(data), dtype_out_time,
                      dtype_out_vert=self.dtype_out_vert)
        # Apply time reduction methods to regional averages and save.
        if any(['reg' in do for do in self.dtype_out_time]) and self.region:
            self._print_verbose("Computing and saving regional outputs.")
            self.region_calcs(full_ts)

    def _save_to_scratch(self, data, dtype_out_time, dtype_out_vert=False):
        """Save the data to the scratch filesystem."""
        print(dtype_out_time)
        path = self.path_scratch[dtype_out_time]
        if not os.path.isdir(self.dir_scratch):
            os.makedirs(self.dir_scratch)
        if 'reg' in dtype_out_time:
            try:
                reg_data = xray.open_dataset(path[:-2] + '.nc')
            except (EOFError, RuntimeError):
                reg_data = xray.Dataset() # Empty Dataset
            # Add the new data to the dictionary or Dataset.
            # Same method works for both.
            reg_data.update(data)
            data_out = reg_data
        else:
            data_out = data
        if isinstance(data_out, xray.DataArray):
            data_out = xray.Dataset({self.name : data_out})
        data_out.to_netcdf(path[:-2] + '.nc', engine='scipy')

    def _save_to_archive(self, dtype_out_time, dtype_out_vert=False):
        """Add the data to the tar file in /archive."""
        if not os.path.isdir(self.dir_archive):
            os.makedirs(self.dir_archive)
        # tarfile 'append' mode won't overwrite the old file, which we want.
        # So open in 'read' mode, extract the file, and then delete it.
        # But 'read' mode throws OSError if file doesn't exist: make it first.
        dmget([self.path_archive])
        # SAH 201-08-27: Within the last few weeks, both of the 'with' blocks
        # below are causing frequent hangups when trying to save to archive
        # for the first time during a shell session.  But then if I Ctrl-C to
        # kill the job and then restart, the hangup goes away.  Don't know
        # why it started happening all of a sudden.
        with tarfile.open(self.path_archive, 'a') as tar:
            pass
        with tarfile.open(self.path_archive, 'r') as tar:
            try:
                old_data_path = '/'.join(
                    [self.dir_archive, self.file_name[dtype_out_time]]
                ).replace('//', '/')

                tar.extract(self.file_name[dtype_out_time],
                            path=old_data_path)
            except KeyError:
                pass
            else:
                # The os module treats files on archive as non-empty
                # directories, so can't use os.remove or os.rmdir.
                shutil.rmtree(old_data_path)
                subprocess.call([
                    "tar", "--delete", "--file=%s" % self.path_archive,
                    self.file_name[dtype_out_time]
                ])
        with tarfile.open(self.path_archive, 'a') as tar:
            tar.add(self.path_scratch[dtype_out_time],
                    arcname=self.file_name[dtype_out_time])

    def _update_data_out(self, data, dtype):
        """Append the data of the given dtype_out to the data_out attr."""
        try:
            self.data_out.update({dtype: data})
        except AttributeError:
            self.data_out = {dtype: data}

    def save(self, data, dtype_out_time, dtype_out_vert=False,
             scratch=True, archive=False):
        """Save aospy data to data_out attr and to an external file."""
        self._update_data_out(data, dtype_out_time)
        #print(self._to_DataArray(data))
        if scratch:
            self._save_to_scratch(data, dtype_out_time,
                                  dtype_out_vert=dtype_out_vert)
        if archive:
            self._save_to_archive(dtype_out_time,
                                  dtype_out_vert=dtype_out_vert)
        print('\t%s' % self.path_scratch[dtype_out_time][:-2] + '.nc')

    def _load_from_scratch(self, dtype_out_time, dtype_out_vert=False):
        """Load aospy data saved on scratch file system."""
        with open(self.path_scratch[dtype_out_time], 'r') as data:
            data_vals = pickle.load(data)
        return data_vals

    def _load_from_archive(self, dtype_out_time, dtype_out_vert=False):
        """Load data save in tarball on archive file system."""
        path = '/'.join([self.dir_archive, 'data.tar']).replace('//', '/')
        dmget([path])
        with tarfile.open(path, 'r') as data_tar:
            data_vals = pickle.load(
                data_tar.extractfile(self.file_name[dtype_out_time])
            )
        return data_vals

    def load(self, dtype_out_time, dtype_out_vert=False, region=False,
             time=False, vert=False, lat=False, lon=False, plot_units=False,
             mask_unphysical=False):
        """Load the data from the object if possible or from disk."""
        # Grab from the object if its there.
        try:
            data = self.data_out[dtype_out_time]
        except (AttributeError, KeyError):
            # Otherwise get from disk.  Try scratch first, then archive.
            try:
                data = self._load_from_scratch(dtype_out_time, dtype_out_vert)
            except OSError:
                data = self._load_from_archive(dtype_out_time, dtype_out_vert)
        # Copy the array to self.data_out for ease of future access.
        self._update_data_out(data, dtype_out_time)
        # Subset the array and convert units as desired.
        if any((region, time, vert, lat, lon)):
            data = self._get_data_subset(data, region=region, time=time,
                                         vert=vert, lat=lat, lon=lon)
        # Apply desired plotting/cleanup methods.
        if mask_unphysical:
            data = self.var.mask_unphysical(data)
        if plot_units:
            data = self.var.to_plot_units(data, vert_int=dtype_out_vert)
        return data

    def _to_DataArray(self, data):
       """
       Converts a Calc instance to an xray DataArray.  Pulls grid information
       from instance.
       """
       if not self.pressure:
           print(data.shape)
           print(self.lat)
           print(self.lon)
           return xray.DataArray(data, coords=[self.lat, self.lon],
                                 dims=['lat','lon'],
                                 encoding={'lat' : 'f8', 'lon' : 'f8'})
       else:
           return None
