"""calc.py: classes for performing specified calculations on aospy data"""
from collections import OrderedDict
import logging
import os
import shutil
import subprocess
import tarfile
from time import ctime

import numpy as np
import pandas as pd
import xarray as xr

from .__config__ import (LAT_STR, LON_STR, LAT_BOUNDS_STR, LON_BOUNDS_STR,
                         PHALF_STR, PFULL_STR, PLEVEL_STR, TIME_STR, YEAR_STR,
                         ETA_STR, BOUNDS_STR)
from .constants import Constant, grav
from . import utils
from .var import Var


logging.basicConfig(level=logging.INFO)

dp = Var(
    name='dp',
    units='Pa',
    domain='atmos',
    description='Pressure thickness of model levels.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False,
)
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
    def _set_data_attrs(self):
        for attr in ('data_start_date',
                     'data_end_date',
                     'default_start_date',
                     'default_end_date',
                     'data_dur',
                     'data_direc',
                     'data_files',
                     'data_dir_struc',
                     'ens_mem_prefix',
                     'ens_mem_ext',
                     'idealized'):
            attr_val = tuple([
                utils.io.get_parent_attr(rn, attr, strict=False)
                for rn in self.run
            ])
            setattr(self, attr, attr_val)

    def __init__(self, proj=None, model=None, run=None, ens_mem=None, var=None,
                 date_range=None, region=None, intvl_in=None, intvl_out=None,
                 dtype_in_time=None, dtype_in_vert=None, dtype_out_time=None,
                 dtype_out_vert=None, level=None, chunk_len=False,
                 verbose=True):
        """Create the CalcInterface object with the given parameters."""
        # 2015-10-13 S. Hill: This tuple-izing is for support of calculations
        # where variables come from different runs.  However, this is a very
        # fragile way of implementing that functionality.  Eventually it will
        # be replaced with something better.
        if not isinstance(run, (list, tuple)):
            run = tuple([run])
        for r in run:
            msg = ("Model '{0}' has no run '{1}'.  Calc object "
                   "will not be generated.".format(model, run))
            if r not in model.runs.values():
                raise AttributeError(msg)
        proj = tuple([proj])
        model = tuple([model])

        # Make tuples the same length.
        if len(proj) == 1 and (len(model) > 1 or len(run) > 1):
            proj = tuple(list(proj)*len(run))
        if len(model) == 1 and len(run) > 1:
            model = tuple(list(model)*len(run))

        self.proj = proj
        self.model = model
        self.run = run

        self._set_data_attrs()

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

        self.months = utils.times.month_indices(intvl_out)
        if date_range == 'default':
            self.start_date = utils.times.ensure_datetime(
                self.run[0].default_start_date)
            self.end_date = utils.times.ensure_datetime(
                self.run[0].default_end_date)
        else:
            self.start_date = utils.times.ensure_datetime(date_range[0])
            self.end_date = utils.times.ensure_datetime(date_range[-1])
        self.date_range = utils.times.create_monthly_time_array(
            self.start_date, self.end_date, self.intvl_out
        )
        # Workaround for limited date range support due to nanosecond
        # precison.  See https://github.com/spencerahill/aospy/issues/98
        self.start_date_xarray = utils.times.numpy_datetime_range_workaround(
            self.start_date)
        self.end_date_xarray = (self.start_date_xarray +
                                (self.end_date - self.start_date))


class Calc(object):
    """Class for executing, saving, and loading a single computation."""

    ARR_XARRAY_NAME = 'aospy_result'

    _grid_attrs = OrderedDict([
        (LAT_STR,        ('lat', 'latitude', 'LATITUDE', 'y', 'yto')),
        (LAT_BOUNDS_STR, ('latb', 'lat_bnds', 'lat_bounds')),
        (LON_STR,        ('lon', 'longitude', 'LONGITUDE', 'x', 'xto')),
        (LON_BOUNDS_STR, ('lonb', 'lon_bnds', 'lon_bounds')),
        ('zsurf',        ('zsurf',)),
        ('sfc_area',     ('area', 'sfc_area')),
        ('land_mask',    ('land_mask',)),
        ('pk',           ('pk',)),
        ('bk',           ('bk',)),
        (PHALF_STR,      ('phalf',)),
        (PFULL_STR,      ('pfull',)),
        (PLEVEL_STR,     ('level', 'lev', 'plev')),
    ])

    def __str__(self):
        """String representation of the object."""
        return "Calc object: " + ', '.join(
            (self.name, self.proj_str, self.model_str, self.run_str_full)
        )

    __repr__ = __str__

    def _dir_out(self):
        """Create string of the data directory to save individual .nc files."""
        ens_label = utils.io.ens_label(self.ens_mem)
        return os.path.join(self.proj[0].direc_out, self.proj_str,
                            self.model_str, self.run_str, ens_label,
                            self.name)

    def _dir_tar_out(self):
        """Create string of the data directory to store a tar file."""
        ens_label = utils.io.ens_label(self.ens_mem)
        return os.path.join(self.proj[0].tar_direc_out,
                            self.proj_str, 'data', self.model_str,
                            self.run_str, ens_label)

    def _file_name(self, dtype_out_time, extension='nc'):
        """Create the name of the aospy file."""
        out_lbl = utils.io.data_out_label(self.intvl_out, dtype_out_time,
                                          dtype_vert=self.dtype_out_vert)
        in_lbl = utils.io.data_in_label(self.intvl_in, self.dtype_in_time,
                                        self.dtype_in_vert)
        ens_lbl = utils.io.ens_label(self.ens_mem)
        yr_lbl = utils.io.yr_label((self.start_date.year, self.end_date.year))
        return '.'.join(
            [self.name, out_lbl, in_lbl, self.model_str, self.run_str,
             ens_lbl, yr_lbl, extension]
        ).replace('..', '.')

    def _path_out(self, dtype_out_time):
        return os.path.join(self.dir_out, self.file_name[dtype_out_time])

    def _path_tar_out(self):
        return os.path.join(self.dir_tar_out, 'data.tar')

    @staticmethod
    def _print_verbose(*args):
        """Print diagnostic message."""
        try:
            return '{0} {1} ({2})'.format(args[0], args[1], ctime())
        except IndexError:
            return '{0} ({1})'.format(args[0], ctime())

    def __init__(self, calc_interface):
        self.__dict__ = vars(calc_interface)
        logging.info(self._print_verbose(
            'Initializing Calc instance:', self.__str__()
        ))

        [mod.set_grid_data() for mod in self.model]

        if isinstance(calc_interface.ens_mem, int):
            self.data_direc = self.data_direc[calc_interface.ens_mem]

        self.dt_set = False

        self.dir_out = self._dir_out()
        self.dir_tar_out = self._dir_tar_out()
        self.file_name = {d: self._file_name(d) for d in self.dtype_out_time}
        self.path_out = {d: self._path_out(d)
                         for d in self.dtype_out_time}
        self.path_tar_out = self._path_tar_out()

        self.data_out = {}

    def _data_files_one_dir(self, name, n=0):
        """Get the file names of the files in a single directory"""
        if self.intvl_in in self.data_files[n]:
            if isinstance(self.data_files[n][self.intvl_in][name], str):
                data_files = [self.data_files[n][self.intvl_in][name]]
            else:
                data_files = self.data_files[n][self.intvl_in][name]
        else:
            if isinstance(self.data_files[n][name], str):
                data_files = [self.data_files[n][name]]
            else:
                data_files = self.data_files[n][name]
        return data_files

    def _get_input_data_paths_one_dir(self, name, data_direc, n=0):
        """Get the names of netCDF files when all in same directory."""
        data_files = self._data_files_one_dir(name, n)
        # data_files may hold absolute or relative paths
        paths = []
        for nc in data_files:
            full = os.path.join(data_direc, nc)
            if os.path.isfile(nc):
                paths.append(nc)
            elif os.path.isfile(full):
                paths.append(full)
            else:
                logging.info("Specified netCDF file `{}` not found".format(nc))
        # Remove duplicate entries.
        files = list(set(paths))
        files.sort()
        return files

    def _get_input_data_paths_gfdl_repo(self, name, n=0):
        """Get the names of netCDF files from a GFDL repo on /archive."""
        return self.model[n].find_data_direc_repo(
            run_name=self.run[n].name, var_name=name
        )

    def _get_input_data_paths_gfdl_dir_struct(self, name, data_direc,
                                              start_year, end_year, n=0):
        """Get paths to netCDF files save in GFDL standard output format."""
        domain = self.domain
        dtype_lbl = self.dtype_in_time
        if self.intvl_in == 'daily':
            domain += '_daily'
        if self.dtype_in_vert == ETA_STR and name != 'ps':
            domain += '_level'
        if self.dtype_in_time == 'inst':
            domain += '_inst'
            dtype_lbl = 'ts'
        if 'monthly_from_' in self.dtype_in_time:
            dtype = self.dtype_in_time.replace('monthly_from_', '')
            dtype_lbl = dtype
        else:
            dtype = self.dtype_in_time
        dur_str = str(self.data_dur[n]) + 'yr'
        if self.dtype_in_time == 'av':
            subdir = self.intvl_in + '_' + dur_str
        else:
            subdir = os.path.join(self.intvl_in, dur_str)
        direc = os.path.join(data_direc, domain, dtype_lbl, subdir)
        files = [os.path.join(direc, utils.io.data_name_gfdl(
                 name, domain, dtype, self.intvl_in, year, self.intvl_out,
                 self.data_start_date[n].year, self.data_dur[n]
                 )) for year in range(start_year, end_year + 1)]
        # Remove duplicate entries.
        files = list(set(files))
        files.sort()
        return files

    def _get_data_direc(self, n):
        if isinstance(self.data_direc, str):
            return self.data_direc
        if isinstance(self.data_direc, (list, tuple)):
            return self.data_direc[n]
        raise IOError("data_direc must be string, list, or tuple: "
                      "{}".format(self.data_direc))

    def _get_input_data_paths(self, var, start_date=False,
                              end_date=False, n=0):
        """Create xarray.DataArray of the variable from its netCDF files.

        Files chosen depend on the specified variables and time interval and
        the attributes of the netCDF files.
        """
        data_direc = self._get_data_direc(n)
        # Cycle through possible names until the data is found.
        for name in var.names:
            if self.data_dir_struc[n] == 'one_dir':
                try:
                    files = self._get_input_data_paths_one_dir(
                        name, data_direc, n=n
                    )
                except KeyError as e:
                    logging.debug(str(repr(e)))
                else:
                    break
            elif self.data_dir_struc[n].lower() == 'gfdl':
                try:
                    files = self._get_input_data_paths_gfdl_dir_struct(
                        name, data_direc, start_date.year,
                        end_date.year, n=n
                    )
                except:
                    raise
                else:
                    break
            elif self.data_dir_struc[n].lower() == 'gfdl_repo':
                try:
                    files = self._get_input_data_paths_gfdl_repo(name, n=n)
                except IOError as e:
                    logging.debug(str(repr(e)))
                else:
                    break
            else:
                raise ValueError("Specified directory type not supported"
                                 ": {}".format(self.data_dir_struc[n]))
        else:
            msg = ("netCDF files for calc object `{0}`, variable `{1}`, year "
                   "range {2}-{3}, in directory {4}, not found")
            raise IOError(msg.format(self, var, start_date, end_date,
                                     data_direc))
        paths = list(set(files))
        paths.sort()
        return paths

    def _to_desired_dates(self, arr):
        """Restrict the xarray DataArray or Dataset to the desired months."""
        times = utils.times.extract_date_range_and_months(
            arr[TIME_STR], self.start_date_xarray, self.end_date_xarray,
            self.months
        )
        return arr.sel(time=times)

    def _add_grid_attributes(self, ds, n=0):
        """Add model grid attributes to a dataset"""
        for name_int, names_ext in self._grid_attrs.items():
            ds_coord_name = set(names_ext).intersection(set(ds.coords) |
                                                        set(ds.data_vars))
            model_attr = getattr(self.model[n], name_int, None)
            if ds_coord_name:
                # Force coords to have desired name.
                ds = ds.rename({list(ds_coord_name)[0]: name_int})
                ds = ds.set_coords(name_int)
                if not np.array_equal(ds[name_int], model_attr):
                    if np.allclose(ds[name_int], model_attr):
                        msg = ("Values for '{0}' are nearly (but not exactly) "
                               "the same in the Run {1} and the Model {2}.  "
                               "Therefore replacing Run's values with the "
                               "model's.".format(name_int, self.run[n],
                                                 self.model[n]))
                        logging.info(msg)
                        ds[name_int].values = model_attr.values
                    else:
                        msg = ("Model coordinates for '{0}' do not match those"
                               " in Run: {1} vs. {2}"
                               "".format(name_int, ds[name_int], model_attr))
                        logging.info(msg)

            else:
                # Bring in coord from model object if it exists.
                ds = ds.load()
                if model_attr is not None:
                    ds[name_int] = model_attr
                    ds = ds.set_coords(name_int)
            if self.dtype_in_vert == 'pressure' and PLEVEL_STR in ds.coords:
                self.pressure = ds.level
        return ds

    @staticmethod
    def dt_from_time_bnds(ds):
        """Compute the timestep durations from the time bounds array."""
        for name in ['time_bounds', 'time_bnds']:
            try:
                bounds = ds[name]
            except KeyError:
                pass
            else:
                dt = bounds.diff(BOUNDS_STR).squeeze().drop(BOUNDS_STR)
                # Convert from float # of days to np.timedelta64 in seconds.
                # TODO: Explicitly check that units are days.
                dt.values = np.array([np.timedelta64(int(d), 'D')
                                      for d in dt.values])
                return dt / np.timedelta64(1, 's')
        raise ValueError("Time bound data cannot be found in the dataset.\n"
                         "{0}".format(ds))

    def _get_dt(self, ds):
        """Find or create the array of timestep durations."""
        for name in ['average_DT']:
            try:
                dt = ds[name]
            except KeyError:
                logging.debug("dt array not found for nonexistent key name "
                              "`{0}`".format(name))
            else:
                # Convert to seconds
                return self._to_desired_dates(dt) / np.timedelta64(1, 's')
        return self._to_desired_dates(self.dt_from_time_bnds(ds))

    def _create_input_data_obj(self, var, start_date=False,
                               end_date=False, n=0, set_dt=False,
                               set_pfull=False):
        """Create xarray.DataArray for the Var from files on disk.

        """
        paths = self._get_input_data_paths(var, start_date, end_date, n)
        # TODO: refactor `dmget` to more general pre-processing step that
        #       user can specify in main or via a config.
        utils.io.dmget(paths)
        ds = xr.open_mfdataset(paths, decode_cf=False)
        # Workaround for years < 1678 causing overflows.
        if start_date < pd.Timestamp.min:
            ds = utils.times.numpy_datetime_workaround_encode_cf(ds)
        ds = xr.decode_cf(ds, decode_times=True)
        ds = self._add_grid_attributes(ds, n)
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
            try:
                self.dt = self._get_dt(ds)
            except ValueError:
                pass
        # At least one variable has to get us the pfull array, if it's needed.
        if set_pfull:
            try:
                self.pfull_coord = ds[PFULL_STR]
            except KeyError:
                pass
        return arr.load()

    def _get_pressure_from_p_coords(self, ps, name='p', n=0):
        """Get pressure or pressure thickness array for data on p-coords."""
        if np.any(self.pressure):
            pressure = self.pressure
        else:
            pressure = self.model[n].level
        if name == 'p':
            return pressure
        if name == 'dp':
            return utils.vertcoord.dp_from_p(pressure, ps)
        raise ValueError("name must be 'p' or 'dp':"
                         "'{}'".format(name))

    def _get_pressure_from_eta_coords(self, ps, name='p', n=0):
        """Get pressure (p) or p thickness array for data on model coords."""
        bk = self.model[n].bk
        pk = self.model[n].pk
        pfull_coord = self.model[n].pfull
        if name == 'p':
            return utils.vertcoord.pfull_from_ps(bk, pk, ps, pfull_coord)
        if name == 'dp':
            return utils.vertcoord.dp_from_ps(bk, pk, ps, pfull_coord)
        raise ValueError("name must be 'p' or 'dp':"
                         "'{}'".format(name))

    def _get_pressure_vals(self, var, start_date, end_date, n=0):
        """Get pressure array, whether sigma or standard levels."""
        try:
            ps = self._ps_data
        except AttributeError:
            self._ps_data = self._create_input_data_obj(self.ps, start_date,
                                                        end_date)
            ps = self._ps_data
        if self.dtype_in_vert == 'pressure':
            return self._get_pressure_from_p_coords(ps, name=var.name, n=n)
        if self.dtype_in_vert == ETA_STR:
            return self._get_pressure_from_eta_coords(ps, name=var.name, n=n)
        raise ValueError("`dtype_in_vert` must be either 'pressure' or "
                         "'sigma' for pressure data")

    def _correct_gfdl_inst_time(self, arr):
        """Correct off-by-one error in GFDL instantaneous model data."""
        if self.intvl_in.endswith('hr'):
            offset = -1*int(self.intvl_in[0])
        else:
            offset = 0
        time = utils.times.apply_time_offset(arr[TIME_STR], hours=offset)
        arr[TIME_STR] = time
        return arr

    def _get_input_data(self, var, start_date, end_date, n):
        """Get the data for a single variable over the desired date range."""
        logging.info(self._print_verbose("Getting input data:", var))
        # If only 1 run, use it to load all data.  Otherwise assume that num
        # runs equals num vars to load.
        if len(self.run) == 1:
            n = 0
        # Pass numerical constants as is.
        if isinstance(var, (float, int)):
            return var
        elif isinstance(var, Constant):
            return var.value
        # aospy.Var objects remain.
        # Pressure handled specially due to complications from sigma vs. p.
        elif var.name in ('p', 'dp'):
            data = self._get_pressure_vals(var, start_date, end_date)
            if self.dtype_in_vert == ETA_STR:
                if self.dtype_in_time == 'inst':
                    data = self._correct_gfdl_inst_time(data)
                return self._to_desired_dates(data)
            return data
        # Get grid, time, etc. arrays directly from model object
        elif var.name in (LAT_STR, LON_STR, TIME_STR, PLEVEL_STR,
                          'pk', 'bk', 'sfc_area'):
            data = getattr(self.model[n], var.name)
        else:
            set_dt = True if not hasattr(self, 'dt') else False
            cond_pfull = ((not hasattr(self, 'pfull')) and var.def_vert and
                          self.dtype_in_vert == ETA_STR)
            data = self._create_input_data_obj(var, start_date, end_date, n=n,
                                               set_dt=set_dt,
                                               set_pfull=cond_pfull)
            # Force all data to be at full pressure levels, not half levels.
            if self.dtype_in_vert == ETA_STR and var.def_vert == 'phalf':
                data = utils.vertcoord.to_pfull_from_phalf(data,
                                                           self.pfull_coord)
        # Correct GFDL instantaneous data time indexing problem.
        if var.def_time:
            if self.dtype_in_time == 'inst':
                data = self._correct_gfdl_inst_time(data)
            # Restrict to the desired dates within each year.
            if self.dtype_in_time != 'av':
                return self._to_desired_dates(data)
        else:
            return data

    def _prep_data(self, data, func_input_dtype):
        """Convert data to type needed by the given function.

        :param data: List of xarray.DataArray objects.
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

    def _get_all_data(self, start_date, end_date):
        """Get the needed data from all of the vars in the calculation."""
        return [self._prep_data(self._get_input_data(var, start_date,
                                                     end_date, n),
                                self.var.func_input_dtype)
                for n, var in enumerate(self.variables)]

    def _local_ts(self, *data):
        """Perform the computation at each gridpoint and time index."""
        arr = self.function(*data)
        if self.var.func_input_dtype == 'numpy':
            arr = xr.DataArray(arr, coords=self.coords)
        arr.name = self.name
        return arr

    def _compute(self, data, monthly_mean=False):
        """Perform the calculation."""
        if monthly_mean:
            data_monthly = []
            for d in data:
                try:
                    data_monthly.append(utils.times.monthly_mean_ts(d))
                except KeyError:
                    data_monthly.append(d)
            data = data_monthly
        local_ts = self._local_ts(*data)
        if self.dtype_in_time == 'inst':
            dt = xr.DataArray(np.ones_like(local_ts[TIME_STR]),
                              dims=[TIME_STR], coords=[local_ts[TIME_STR]])
            if not hasattr(self, 'dt'):
                self.dt = dt
        else:
            if hasattr(self, 'dt'):
                dt = self.dt
            else:
                logging.warning("dt array not found.  Assuming equally spaced "
                                "values in time, even though this may not be "
                                "the case")
                dt = xr.DataArray(np.ones(np.shape(local_ts[TIME_STR])),
                                  dims=[TIME_STR], coords=[local_ts[TIME_STR]])
                self.dt = dt
        if monthly_mean:
            dt = utils.times.monthly_mean_ts(dt)
        return local_ts, dt

    # TODO: Move to utils.vertcoord
    def _vert_int(self, arr, dp):
        """Vertical integral"""
        return utils.vertcoord.int_dp_g(arr, dp)

    def _compute_full_ts(self, data, monthly_mean=False, zonal_asym=False):
        """Perform calculation and create yearly timeseries at each point."""
        # Get results at each desired timestep and spatial point.
        # Here we need to provide file read-in dates (NOT xarray dates)
        full_ts, dt = self._compute(data, monthly_mean=monthly_mean)
        if zonal_asym:
            full_ts = full_ts - full_ts.mean(LON_STR)
        # Vertically integrate.
        vert_types = ('vert_int', 'vert_av')
        if self.dtype_out_vert in vert_types and self.var.def_vert:
            # Here we need file read-in dates (NOT xarray dates)
            full_ts = self._vert_int(full_ts, self._get_pressure_vals(
                dp, self.start_date, self.end_date
            ))
            if self.dtype_out_vert == 'vert_av':
                full_ts *= (grav.value / self._to_desired_dates(self._ps_data))
        return full_ts, dt

    def _avg_by_year(self, arr, dt):
        """Average a sub-yearly time-series over each year."""
        yr_str = TIME_STR + '.year'
        return ((arr*dt).groupby(yr_str).sum(TIME_STR) /
                dt.groupby(yr_str).sum(TIME_STR))

    def _full_to_yearly_ts(self, arr, dt):
        """Average the full timeseries within each year."""
        time_defined = self.def_time and not ('av' in self.dtype_in_time or
                                              self.idealized[0])
        if time_defined:
            arr = self._avg_by_year(arr, dt)
        return arr

    def _time_reduce(self, arr, reduction):
        """Perform the specified time reduction on a local time-series."""
        if self.dtype_in_time == 'av':
            return arr
        reductions = {
            'None': lambda xarr: xarr,
            'ts': lambda xarr: xarr,
            'av': lambda xarr: xarr.mean(YEAR_STR),
            'std': lambda xarr: xarr.std(YEAR_STR),
            }
        try:
            return reductions[reduction](arr)
        except KeyError:
            raise ValueError("Specified time-reduction method '{}' is not "
                             "supported".format(reduction))

    def region_calcs(self, arr, func, n=0):
        """Perform a calculation for all regions."""
        # Get pressure values for data output on hybrid vertical coordinates.
        bool_pfull = (self.def_vert and self.dtype_in_vert == ETA_STR and
                      self.dtype_out_vert is False)
        if bool_pfull:
            pfull = self._full_to_yearly_ts(self._prep_data(
                self._get_input_data(Var('p'), self.start_date, self.end_date,
                                     0), self.var.func_input_dtype
            ), self.dt).rename('pressure')
        # Loop over the regions, performing the calculation.
        reg_dat = {}
        for reg in self.region.values():
            # Just pass along the data if averaged already.
            if 'av' in self.dtype_in_time:
                data_out = reg.ts(arr)
            # Otherwise perform the calculation.
            else:
                method = getattr(reg, func)
                data_out = method(arr)
                if bool_pfull:
                    # Don't apply e.g. standard deviation to coordinates.
                    if func not in ['av', 'ts']:
                        method = reg.ts
                    # Convert Pa to hPa
                    coord = method(pfull) * 1e-2
                    data_out = data_out.assign_coords(
                        **{reg.name + '_pressure': coord}
                    )
            reg_dat.update(**{reg.name: data_out})
        return OrderedDict(sorted(reg_dat.items(), key=lambda t: t[0]))

    def _apply_all_time_reductions(self, full_ts, monthly_ts, eddy_ts):
        """Apply all requested time reductions to the data."""
        logging.info(self._print_verbose("Applying desired time-"
                                         "reduction methods."))
        # Determine which are regional, eddy, time-mean.
        reduc_specs = [r.split('.') for r in self.dtype_out_time]
        reduced = {}
        for reduc, specs in zip(self.dtype_out_time, reduc_specs):
            func = specs[-1]
            if 'eddy' in specs:
                data = eddy_ts
            elif 'time-mean' in specs:
                data = monthly_ts
            else:
                data = full_ts
            if 'reg' in specs:
                reduced.update({reduc: self.region_calcs(data, func)})
            else:
                reduced.update({reduc: self._time_reduce(data, func)})
        return OrderedDict(sorted(reduced.items(), key=lambda t: t[0]))

    def _make_full_mean_eddy_ts(self, data):
        """Create full, monthly-mean, and eddy timeseries of data."""
        bool_monthly = (['monthly_from' in self.dtype_in_time] +
                        ['time-mean' in dout for dout in self.dtype_out_time])
        bool_eddy = ['eddy' in dout for dout in self.dtype_out_time]
        if not all(bool_monthly):
            full, full_dt = self._compute_full_ts(data,
                                                  monthly_mean=False)
        else:
            full = False
        if any(bool_eddy) or any(bool_monthly):
            monthly, monthly_dt = self._compute_full_ts(data,
                                                        monthly_mean=True)
        else:
            monthly = False
        if any(bool_eddy):
            eddy = full - utils.times.monthly_mean_at_each_ind(monthly, full)
        else:
            eddy = False

        # Average within each year.
        if not all(bool_monthly):
            full = self._full_to_yearly_ts(full, full_dt)
        if any(bool_monthly):
            monthly = self._full_to_yearly_ts(monthly, monthly_dt)
        if any(bool_eddy):
            eddy = self._full_to_yearly_ts(eddy, full_dt)
        return full, monthly, eddy

    def compute(self, save_files=True, save_tar_files=True):
        """Perform all desired calculations on the data and save externally."""
        data = self._prep_data(self._get_all_data(self.start_date,
                                                  self.end_date),
                               self.var.func_input_dtype)
        logging.info('Computing timeseries for {0} -- '
                     '{1}.'.format(self.start_date, self.end_date))
        full, monthly, eddy = self._make_full_mean_eddy_ts(data)
        reduced = self._apply_all_time_reductions(full, monthly, eddy)
        logging.info("Writing desired gridded outputs to disk.")
        for dtype_time, data in reduced.items():
            self.save(data, dtype_time, dtype_out_vert=self.dtype_out_vert,
                      save_files=save_files, save_tar_files=save_tar_files)

    def _save_files(self, data, dtype_out_time):
        """Save the data to netcdf files in direc_out."""
        path = self.path_out[dtype_out_time]
        if not os.path.isdir(self.dir_out):
            os.makedirs(self.dir_out)
        if 'reg' in dtype_out_time:
            try:
                reg_data = xr.open_dataset(path)
            except (EOFError, RuntimeError, IOError):
                reg_data = xr.Dataset()
            # Add the new data to the dictionary or Dataset.
            # Same method works for both.
            reg_data.update(data)
            data_out = reg_data
        else:
            data_out = data
        if isinstance(data_out, xr.DataArray):
            data_out = xr.Dataset({self.name: data_out})
        data_out.to_netcdf(path, engine='scipy')

    def _save_tar_files(self, dtype_out_time):
        """Add the data to the tar file in tar_out_direc."""
        if not os.path.isdir(self.dir_tar_out):
            os.makedirs(self.dir_tar_out)
        # tarfile 'append' mode won't overwrite the old file, which we want.
        # So open in 'read' mode, extract the file, and then delete it.
        # But 'read' mode throws OSError if file doesn't exist: make it first.
        utils.io.dmget([self.path_tar_out])
        with tarfile.open(self.path_tar_out, 'a') as tar:
            pass
        with tarfile.open(self.path_tar_out, 'r') as tar:
            old_data_path = os.path.join(self.dir_tar_out,
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
                    "tar", "--delete", "--file={}".format(self.path_tar_out),
                    self.file_name[dtype_out_time]
                ])
        with tarfile.open(self.path_tar_out, 'a') as tar:
            tar.add(self.path_out[dtype_out_time],
                    arcname=self.file_name[dtype_out_time])

    def _update_data_out(self, data, dtype):
        """Append the data of the given dtype_out to the data_out attr."""
        try:
            self.data_out.update({dtype: data})
        except AttributeError:
            self.data_out = {dtype: data}

    def save(self, data, dtype_out_time, dtype_out_vert=False,
             save_files=True, save_tar_files=False):
        """Save aospy data to data_out attr and to an external file."""
        self._update_data_out(data, dtype_out_time)
        if save_files:
            self._save_files(data, dtype_out_time)
        if save_tar_files and self.proj[0].tar_direc_out:
            self._save_tar_files(dtype_out_time)
        logging.info('\t{}'.format(self.path_out[dtype_out_time]))

    def _load_from_disk(self, dtype_out_time, dtype_out_vert=False,
                        region=False):
        """Load aospy data saved as netcdf files on the file system."""
        ds = xr.open_dataset(self.path_out[dtype_out_time],
                             engine='scipy')
        if region:
            arr = ds[region.name]
            # Use region-specific pressure values if available.
            if self.dtype_in_vert == ETA_STR and not dtype_out_vert:
                reg_pfull_str = region.name + '_pressure'
                arr = arr.drop([r for r in arr.coords.iterkeys()
                                if r not in (PFULL_STR, reg_pfull_str)])
                # Rename pfull to pfull_ref always.
                arr = arr.rename({PFULL_STR: PFULL_STR + '_ref'})
                # Rename region_pfull to pfull if its there.
                if hasattr(arr, reg_pfull_str):
                    return arr.rename({reg_pfull_str: PFULL_STR})
                return arr
            return arr
        return ds[self.name]

    def _load_from_tar(self, dtype_out_time, dtype_out_vert=False):
        """Load data save in tarball form on the file system."""
        path = os.path.join(self.dir_tar_out, 'data.tar')
        utils.io.dmget([path])
        with tarfile.open(path, 'r') as data_tar:
            ds = xr.open_dataset(
                data_tar.extractfile(self.file_name[dtype_out_time]),
                engine='scipy'
            )
            return ds[self.name]

    def _get_data_subset(self, data, region=False, time=False,
                         vert=False, lat=False, lon=False, n=0):
        """Subset the data array to the specified time/level/lat/lon, etc."""
        if region:
            raise NotImplementedError
        if np.any(time):
            data = data[time]
            if 'monthly_from_' in self.dtype_in_time:
                data = np.mean(data, axis=0)[np.newaxis, :]
        if np.any(vert):
            if self.dtype_in_vert == ETA_STR:
                data = data[{PFULL_STR: vert}]
            else:
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
        msg = ("Loading data from disk for object={0}, dtype_out_time={1}, "
               "dtype_out_vert={2}, and region="
               "{3}".format(self, dtype_out_time, dtype_out_vert, region))
        logging.info(msg + ' ({})'.format(ctime()))
        # Grab from the object if its there.
        try:
            data = self.data_out[dtype_out_time]
        except (AttributeError, KeyError):
            # Otherwise get from disk.  Try scratch first, then archive.
            try:
                data = self._load_from_disk(dtype_out_time, dtype_out_vert,
                                            region=region)
            except IOError:
                data = self._load_from_tar(dtype_out_time, dtype_out_vert)
        # Copy the array to self.data_out for ease of future access.
        self._update_data_out(data, dtype_out_time)
        # Subset the array and convert units as desired.
        if any((time, vert, lat, lon)):
            data = self._get_data_subset(data, region=False, time=time,
                                         vert=vert, lat=lat, lon=lon)
        # Apply desired plotting/cleanup methods.
        if mask_unphysical:
            data = self.var.mask_unphysical(data)
        if plot_units:
            data = self.var.to_plot_units(data, dtype_vert=dtype_out_vert)
        return data
