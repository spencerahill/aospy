"""calc.py: classes for performing specified calculations on aospy data"""
from __future__ import print_function
import os
import shutil
import subprocess
import tarfile
import time

import numpy as np
import xray

from . import Constant, Var
from .io import (_data_in_label, _data_out_label, _ens_label, _yr_label, dmget,
                 data_in_name_gfdl)
from .timedate import TimeManager, _get_time
from .utils import (get_parent_attr, level_thickness,
                    pfull_from_sigma, dp_from_sigma, int_dp_g)

TIME_STR = 'time'

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
    def _set_data_in_attrs(self):
        for attr in ('data_in_start_date',
                     'data_in_end_date',
                     'default_date_range',
                     'data_in_dur',
                     'data_in_direc',
                     'data_in_files',
                     'data_in_dir_struc',
                     'ens_mem_prefix',
                     'ens_mem_ext',
                     'idealized'):
            attr_val = tuple([get_parent_attr(rn, attr, strict=False)
                              for rn in self.run])
            setattr(self, attr, attr_val)

    def __init__(self, proj=None, model=None, run=None, ens_mem=None, var=None,
                 date_range=None, region=None, intvl_in=None, intvl_out=None,
                 dtype_in_time=None, dtype_in_vert=None, dtype_out_time=None,
                 dtype_out_vert=None, level=None, chunk_len=False,
                 verbose=True):
        """Create the CalcInterface object with the given parameters."""
        if run not in model.runs.values():
            raise AttributeError("Model '{}' has no run '{}'.  Calc object "
                                 "will not be generated.".format(model, run))
        # 2015-10-13 S. Hill: This tuple-izing is for support of calculations
        # where variables come from different runs.  However, this is a very
        # fragile way of implementing that functionality.  Eventually it will
        # be replaced with something better.
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

        self._set_data_in_attrs()

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
        if isinstance(dtype_out_time, (list, tuple)):
            self.dtype_out_time = tuple(dtype_out_time)
        else:
            self.dtype_out_time = tuple([dtype_out_time])
        self.dtype_out_vert = dtype_out_vert
        self.region = region

        self.months = TimeManager.month_indices(intvl_out)
        self.start_date = TimeManager.str_to_datetime(date_range[0])
        self.end_date = TimeManager.str_to_datetime(date_range[-1])
        tm = TimeManager(self.start_date, self.end_date, intvl_out)
        self.date_range = tm.create_time_array()

        self.start_date_xray = tm.apply_year_offset(self.start_date)
        self.end_date_xray = tm.apply_year_offset(self.end_date)


class Calc(object):
    """Class for executing, saving, and loading a single computation."""

    ARR_XRAY_NAME = 'aospy_result'

    def __str__(self):
        """String representation of the object."""
        return "Calc object: " + ', '.join(
            (self.name, self.proj_str, self.model_str, self.run_str_full)
        )

    __repr__ = __str__

    def _dir_scratch(self):
        """Create string of the data directory on the scratch filesystem."""
        ens_label = _ens_label(self.ens_mem)
        return os.path.join('/work', os.getenv('USER'), self.proj_str,
                            self.model_str, self.run_str, ens_label,
                            self.name)

    def _dir_archive(self):
        """Create string of the data directory on the archive filesystem."""
        ens_label = _ens_label(self.ens_mem)
        return os.path.join('/archive', os.getenv('USER'),
                            self.proj_str, 'data', self.model_str,
                            self.run_str, ens_label)

    def _file_name(self, dtype_out_time, extension='nc'):
        """Create the name of the aospy file."""
        out_lbl = _data_out_label(self.intvl_out, dtype_out_time,
                                  dtype_vert=self.dtype_out_vert)
        in_lbl = _data_in_label(self.intvl_in, self.dtype_in_time,
                                self.dtype_in_vert)
        ens_lbl = _ens_label(self.ens_mem)
        yr_lbl = _yr_label((self.start_date.year,
                            self.end_date.year))
        return '.'.join(
            [self.name, out_lbl, in_lbl, self.model_str, self.run_str_full,
             ens_lbl, yr_lbl, extension]
        ).replace('..', '.')

    def _path_scratch(self, dtype_out_time):
        return os.path.join(self.dir_scratch, self.file_name[dtype_out_time])

    def _path_archive(self):
        return os.path.join(self.dir_archive, 'data.tar')

    def _print_verbose(self, *args):
        """Print diagnostic message."""
        if not self.verbose:
            pass
        else:
            try:
                print('{} {}'.format(args[0], args[1]),
                      '({})'.format(time.ctime()))
            except IndexError:
                print('{}'.format(args[0]), '({})'.format(time.ctime()))

    def __init__(self, calc_interface):
        self.__dict__ = vars(calc_interface)
        print()
        self._print_verbose('Initializing Calc instance:', self.__str__())

        [mod.set_grid_data() for mod in self.model]

        if isinstance(calc_interface.ens_mem, int):
            self.data_in_direc = self.data_in_direc[calc_interface.ens_mem]

        self.dt_set = False

        self.dir_scratch = self._dir_scratch()
        self.dir_archive = self._dir_archive()
        self.file_name = {d: self._file_name(d) for d in self.dtype_out_time}
        self.path_scratch = {d: self._path_scratch(d)
                             for d in self.dtype_out_time}
        self.path_archive = self._path_archive()

        self.data_out = {}

    def _get_input_data_paths_one_dir(self, name, data_in_direc, n=0):
        """Get the names of netCDF files when all in same directory."""
        if isinstance(self.data_in_files[n][name], str):
            data_in_files = [self.data_in_files[n][name]]
        else:
            data_in_files = self.data_in_files[n][name]
        # data_in_files may hold absolute or relative paths
        paths = []
        for nc in data_in_files:
            full = '/'.join([data_in_direc, nc]).replace('//', '/')
            if os.path.isfile(nc):
                paths.append(nc)
            elif os.path.isfile(full):
                paths.append(full)
            else:
                print("Warning: specified netCDF file `{}` "
                      "not found".format(nc))
        # Remove duplicate entries.
        files = list(set(paths))
        files.sort()
        return files

    def _get_input_data_paths_gfdl_repo(self, name, n=0):
        """Get the names of netCDF files from a GFDL repo on /archive."""
        return self.model[n].find_data_in_direc_repo(
            run_name=self.run[n].name, var_name=name
        )

    def _get_input_data_paths_gfdl_dir_struct(self, name, data_in_direc,
                                              start_year, end_year, n=0):
        """Get paths to netCDF files save in GFDL standard output format."""
        domain = self.domain
        dtype_lbl = self.dtype_in_time
        if self.dtype_in_vert == 'sigma' and name != 'ps':
            domain += '_level'
        if self.dtype_in_time == 'inst':
            domain += '_inst'
            dtype_lbl = 'ts'
        if 'av_from_' in self.dtype_in_time:
            dtype = self.dtype_in_time.replace('av_from_', '')
            dtype_lbl = dtype
        else:
            dtype = self.dtype_in_time
        direc = os.path.join(data_in_direc, domain, dtype_lbl, self.intvl_in,
                             str(self.data_in_dur[n]) + 'yr')
        files = [os.path.join(direc, data_in_name_gfdl(name, domain, dtype,
                                                  self.intvl_in, year,
                                                  self.intvl_out,
                                                  self.data_in_start_date[n].year,
                                                  self.data_in_dur[n]))
                 for year in range(start_year, end_year + 1)]
        # Remove duplicate entries.
        files = list(set(files))
        files.sort()
        return files

    def _get_data_in_direc(self, n):
        if isinstance(self.data_in_direc, str):
            return self.data_in_direc
        if isinstance(self.data_in_direc, (list, tuple)):
            return self.data_in_direc[n]
        raise IOError("data_in_direc must be string, list, or tuple: "
                      "{}".format(self.data_in_direc))

    def _get_input_data_paths(self, var, start_date=False,
                              end_date=False, n=0):
        """Create xray.DataArray of the variable from its netCDF files on disk.

        Files chosen depend on the specified variables and time interval and
        the attributes of the netCDF files.
        """
        data_in_direc = self._get_data_in_direc(n)
        # Cycle through possible names until the data is found.
        for name in var.names:
            if self.data_in_dir_struc[n] == 'one_dir':
                try:
                    files = self._get_input_data_paths_one_dir(
                        name, data_in_direc, n=n
                    )
                except KeyError:
                    pass
                else:
                    break
            elif self.data_in_dir_struc[n].lower() == 'gfdl':
                try:
                    files = self._get_input_data_paths_gfdl_dir_struct(
                        name, data_in_direc, start_date.year,
                        end_date.year, n=n
                    )
                except:
                    raise
                else:
                    break
            elif self.data_in_dir_struc[n].lower() == 'gfdl_repo':
                try:
                    files = self._get_input_data_paths_gfdl_repo(name, n=n)
                except IOError:
                    pass
                else:
                    break
            else:
                raise ValueError("Specified directory type not supported"
                                 ": {}".format(self.data_in_dir_struc[n]))
        else:
            raise IOError("netCDF files for variable `{}`, year range {}-{}, "
                          "in directory {}, not found".format(var, start_date,
                                                              end_date,
                                                              data_in_direc))
        paths = list(set(files))
        paths.sort()
        return paths

    def _create_input_data_obj(self, var, start_date=False,
                               end_date=False, n=0, set_dt=False):
        """Create xray.DataArray for the Var from files on disk."""
        paths = self._get_input_data_paths(var, start_date, end_date, n)
        # 2015-10-15 S. Hill: This `dmget` call, which is unique to the
        # filesystem at the NOAA GFDL computing cluster, should be factored out
        # of this function.  A config setting or some other user input should
        # specify what method to call to access the files on the filesystem.
        dmget(paths)
        ds_chunks = []
        # 2015-10-16 S. Hill: Can we use the xray.open_mfdataset function here
        # instead of this logic of making individual datasets and then
        # calling xray.concat?  Or does the year<1678 logic make this not
        # possible?

        # 2015-10-16 19:06:00 S. Clark: The year<1678 logic is independent of using
        # xray.open_mfdataset. The main reason I held off on using it here
        # was that it opens a can of worms with regard to performance; we'd
        # need to add some logic to make sure the data were chunked in a
        # reasonable way (and logic to change the chunking if need be).
        for file_ in paths:
            test = xray.open_dataset(file_, decode_cf=False,
                                     drop_variables=['time_bounds', 'nv',
                                                     'average_T1',
                                                     'average_T2'])
            if start_date.year < 1678:
                for v in ['time']:
                    test[v].attrs['units'] = ('days since 1900-01-01 '
                                              '00:00:00')
                test['time'].attrs['calendar'] = 'noleap'
            test = xray.decode_cf(test)
            ds_chunks.append(test)
        ds = xray.concat(ds_chunks, dim='time')
        # 2015-10-16 S. Hill: Passing in each variable as a Dataset is causing
        # lots of problems in my functions, even ones as simple as just adding
        # the two variables together.  I think it makes most sense to just grab
        # the DataArray of the desired data from the Dataset.
        for name in var.names:
            try:
                arr = ds[name]
            except KeyError:
                pass
            else:
                break
        else:
            raise KeyError('Variable not found: {}'.format(var))
        # At least one variable has to get us the dt array also.
        if set_dt:
            for name in ['average_DT']:
                try:
                    dt = ds[name]
                except KeyError:
                    pass
                else:
                    self.dt = self._to_desired_dates(dt)
                    break
        return arr

    def _phalf_to_pfull(self, arr):
        """Interpolate data at sigma half levels to full levels."""
        raise NotImplementedError

    def _get_pressure_vals(self, var, start_date, end_date, n=0):
        """Get pressure array, whether sigma or standard levels."""
        self._print_verbose("Getting pressure data:", var)
        if self.dtype_in_vert == 'pressure':
            if np.any(self.pressure):
                pressure = self.pressure
            else:
                pressure = self.model[n].level
            if var == 'p':
                data = pressure
            elif var == 'dp':
                data = level_thickness(pressure)

        if self.dtype_in_vert == 'sigma':
            bk = self.model[n].bk
            pk = self.model[n].pk
            ps = self._create_input_data_obj(self.ps, start_date, end_date)
            pfull_coord = self.model[n].pfull
            if var == 'p':
                data = pfull_from_sigma(bk, pk, ps, pfull_coord)
            elif var == 'dp':
                data = dp_from_sigma(bk, pk, ps, pfull_coord)
        return data

    def _to_desired_dates(self, arr):
        """Restrict the xray DataArray or Dataset to the desired months."""
        times = _get_time(arr['time'], self.start_date_xray,
                          self.end_date_xray, self.months, indices=False)
        return arr.sel(time=times)

    def _get_input_data(self, var, start_date, end_date, n):
        # If only 1 run, use it to load all data.
        # Otherwise assume that # runs == # vars to load.
        if len(self.run) == 1:
            n = 0
        # Pressure handled specially due to complications from sigma vs. p.
        if var in ('p', 'dp'):
            data = self._get_pressure_vals(var, start_date, end_date)
        # Pass numerical constants as is.
        elif isinstance(var, (float, int)):
            data = var
        elif isinstance(var, Constant):
            data = var.value
        # aospy.Var objects remain.
        # Get grid, time, etc. arrays directly from model object
        elif var.name in ('lat', 'lon', 'time', 'level',
                          'pk', 'bk', 'sfc_area'):
            data = getattr(self.model[n], var.name)
        else:
            set_dt = True if not hasattr(self, 'dt') else False
            data = self._create_input_data_obj(var, start_date, end_date, n=n,
                                               set_dt=set_dt)
        # Force all data to be at full pressure levels, not half levels.
        phalf_bool = (self.dtype_in_vert == 'sigma' and not
                      isinstance(var, (Constant, str)) and
                      var.def_vert == 'phalf')
        if phalf_bool:
                data = self._phalf_to_pfull(data)
        # Restrict to the desired dates within each year.
        if isinstance(var, Constant):
            return data
        return self._to_desired_dates(data)

    def _get_all_data(self, start_date, end_date):
        """Get the needed data from all of the vars in the calculation."""
        return [self._get_input_data(v, start_date, end_date, n)
                for n, v in enumerate(self.variables)]

    def _prep_data(self, data, func_input_dtype):
        """Convert data to type needed by the given function.

        :param data: List of xray.DataArray objects.
        :param func_input_dtype: One of (None, 'DataArray', 'Dataset',
                                 'numpy'). Specifies which datatype to convert
                                 to.
        """
        if func_input_dtype is None:
            return data
        if func_input_dtype == 'DataArray':
            return data
        if func_input_dtype == 'Dataset':
            # S. Hill 2015-10-19: This should be filled in with logic that
            # creates a single Dataset comprising all of the DataArray objects
            # in `data`.
            return NotImplementedError
        if func_input_dtype == 'numpy':
            self.coords = data[0].coords
            return [d.values for d in data]
        raise ValueError("Unknown func_input_dtype "
                         "'{}'.".format(func_input_dtype))

    def _local_ts(self, *data_in):
        """Perform the computation at each gridpoint and time index."""
        arr = self.function(*data_in)
        if self.var.func_input_dtype == 'numpy':
            arr = xray.DataArray(arr, coords=self.coords)
        arr.name = self.name
        return arr

    def _compute(self, start_date, end_date):
        """Perform the calculation."""
        self._print_verbose('\n', 'Computing desired timeseries for years '
                            '{}-{}.'.format(start_date.year, end_date.year))
        data_in = self._prep_data(self._get_all_data(start_date, end_date),
                                  self.var.func_input_dtype)
        # 2015-10-16 S. Hill: At this step, maybe we want to combine the
        # Dataset objects, where there's one per variable, into a single one,
        # retaining the shared coordinate etc. arrays.  Then we can support
        # functions that accept a single dataset with all the data.  We could
        # also separate it out into individual DataArrays if the function takes
        # separate DataArrays, or into numpy arrays using .values if the
        # function requires numpy arrays.  However, this isn't implemented yet.
        # Right now, the Datasets for each Var are just combined into a single
        # list.
        local_ts = self._local_ts(*data_in)
        if self.dtype_in_time == 'inst':
            dt = xray.DataArray(np.ones(np.shape(local_ts[TIME_STR])),
                                dims=[TIME_STR], coords=[local_ts[TIME_STR]])
        else:
            dt = self.dt
        return local_ts, dt

    def _to_yearly_ts(self, arr, dt):
        """Average a sub-yearly time-series over each year."""
        # Convert from ns to days (prevent overflow)
        dt.values = dt.values.astype('timedelta64[D]').astype('float')
        return ((arr*dt).groupby('time.year').sum('time') /
                dt.groupby('time.year').sum('time'))

    def _vert_int(self, arr, dp):
        """Vertical integral"""
        return int_dp_g(arr, dp)

    def _time_reduce(self, loc_ts):
        """Compute all desired calculations on a local time-series."""
        tdim = 'time' if self.idealized[0] else 'year'
        files = {}
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
        calcs_reg = ('ts', 'av', 'std')
        # Perform each calculation for each region.
        for calc in calcs_reg:
            calc_name = ('reg.' + calc)
            if calc_name in self.dtype_out_time:
                reg_dat = {}
                for reg in self.region.values():
                    # Just pass along the data if averaged already.
                    if 'av' in self.dtype_in_time:
                        data_out = reg.ts(loc_ts, self.model[n])
                    # Otherwise perform the calculation.
                    else:
                        method = getattr(reg, calc)
                        data_out = method(loc_ts, self.model[n])
                    reg_dat.update({reg.name: data_out})
                self.save(reg_dat, calc_name)

    def compute(self):
        """Perform all desired calculations on the data and save externally."""
        # Get results at each desired timestep and spatial point.
        # Here we need to provide file read-in dates (NOT xray dates)
        full_ts, dt = self._compute(self.start_date, self.end_date)
        # Average within each year if time-defined.
        # (If we don't want to group by year not_def_time = True)
        not_def_time = ('av' in self.dtype_in_time or not self.def_time or
                        self.idealized[0])
        if not not_def_time:
            full_ts = self._to_yearly_ts(full_ts, dt)
        # Vertically integrate if vertically defined and specified.
        if self.dtype_out_vert == 'vert_int' and self.var.def_vert:
            # Here we need file read-in dates (NOT xray dates)
            dp = self._get_pressure_vals('dp', self.start_date,
                                         self.end_date)
            full_ts = self._vert_int(full_ts, dp)
        # Apply time reduction methods and save.
        if self.def_time:
            self._print_verbose("Applying desired time-reduction methods.")
            reduced = self._time_reduce(full_ts)
        else:
            reduced = {'': full_ts}
        self._print_verbose("Writing desired gridded outputs to disk.")
        for dtype_out_time, data in reduced.items():
            self.save(data, dtype_out_time, dtype_out_vert=self.dtype_out_vert)
        # Apply time reduction methods to regional averages and save.
        if any(['reg' in do for do in self.dtype_out_time]) and self.region:
            self._print_verbose("Computing and saving regional outputs.")
            self.region_calcs(full_ts)

    def _save_to_scratch(self, data, dtype_out_time):
        """Save the data to the scratch filesystem."""
        path = self.path_scratch[dtype_out_time]
        if not os.path.isdir(self.dir_scratch):
            os.makedirs(self.dir_scratch)
        if 'reg' in dtype_out_time:
            try:
                reg_data = xray.open_dataset(path)
            except (EOFError, RuntimeError):
                reg_data = xray.Dataset()
            # Add the new data to the dictionary or Dataset.
            # Same method works for both.
            reg_data.update(data)
            data_out = reg_data
        else:
            data_out = data
        if isinstance(data_out, xray.DataArray):
            data_out = xray.Dataset({self.name: data_out})
        data_out.to_netcdf(path, engine='scipy')

    def _save_to_archive(self, dtype_out_time):
        """Add the data to the tar file in /archive."""
        if not os.path.isdir(self.dir_archive):
            os.makedirs(self.dir_archive)
        # tarfile 'append' mode won't overwrite the old file, which we want.
        # So open in 'read' mode, extract the file, and then delete it.
        # But 'read' mode throws OSError if file doesn't exist: make it first.
        dmget([self.path_archive])
        with tarfile.open(self.path_archive, 'a') as tar:
            pass
        with tarfile.open(self.path_archive, 'r') as tar:
            old_data_path = os.path.join(self.dir_archive,
                                         self.file_name[dtype_out_time])
            try:
                tar.extract(self.file_name[dtype_out_time],
                            path=old_data_path)
            except KeyError:
                pass
            else:
                # The os module treats files on archive as non-empty
                # directories, so can't use os.remove or os.rmdir.
                shutil.rmtree(old_data_path)
                subprocess.call([
                    "tar", "--delete", "--file={}".format(self.path_archive),
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
        if scratch:
            self._save_to_scratch(data, dtype_out_time)
        if archive:
            self._save_to_archive(dtype_out_time)
        print('\t', '{}'.format(self.path_scratch[dtype_out_time]))

    def _load_from_scratch(self, dtype_out_time, dtype_out_vert=False):
        """Load aospy data saved on scratch file system."""
        ds = xray.open_dataset(self.path_scratch[dtype_out_time],
                               engine='scipy')
        return ds[self.name]

    def _load_from_archive(self, dtype_out_time, dtype_out_vert=False):
        """Load data save in tarball on archive file system."""
        path = os.path.join(self.dir_archive, 'data.tar')
        dmget([path])
        with tarfile.open(path, 'r') as data_tar:
            ds = xray.open_dataset(
                data_tar.extractfile(self.file_name[dtype_out_time]),
                engine='scipy'
            )
            return ds[self.name]

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
                data = np.mean(data, axis=0)[np.newaxis, :]
        if np.any(vert):
            if self.dtype_in_vert != 'sigma':
                if np.max(self.model[n].level) > 1e4:
                    # Convert from Pa to hPa.
                    lev_hpa = self.model[n].level*1e-2
                else:
                    lev_hpa = self.model[n].level
                level_index = np.where(lev_hpa == self.level)
                if 'ts' in self.dtype_out_time:
                    data = np.squeeze(data[:, level_index])
                else:
                    data = np.squeeze(data[level_index])
        if np.any(lat):
            raise NotImplementedError
        if np.any(lon):
            raise NotImplementedError
        return data

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
            except IOError:
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
