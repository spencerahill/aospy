"""Functionality for performing user-specified calculations on aospy data."""
from collections import OrderedDict
import logging
import os
import shutil
import subprocess
import tarfile
from time import ctime

import numpy as np
import xarray as xr

from ._constants import GRAV_EARTH
from .var import Var
from . import internal_names
from . import utils

logging.basicConfig(level=logging.INFO)


_P_VARS = {internal_names.ETA_STR: utils.vertcoord.p_eta,
           'pressure': utils.vertcoord.p_level}
_DP_VARS = {internal_names.ETA_STR: utils.vertcoord.dp_eta,
            'pressure': utils.vertcoord.dp_level}
_TIME_DEFINED_REDUCTIONS = ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts']


def _replace_pressure(arguments, dtype_in_vert):
    """Replace p and dp Vars with appropriate Var objects specific to
    the dtype_in_vert."""
    arguments_out = []
    for arg in arguments:
        if isinstance(arg, Var):
            if arg.name == 'p':
                arguments_out.append(_P_VARS[dtype_in_vert])
            elif arg.name == 'dp':
                arguments_out.append(_DP_VARS[dtype_in_vert])
            else:
                arguments_out.append(arg)
        else:
            arguments_out.append(arg)
    return arguments_out


class Calc(object):
    """Class for executing, saving, and loading a single computation."""

    ARR_XARRAY_NAME = 'aospy_result'

    _grid_coords = [internal_names.LAT_STR, internal_names.LAT_BOUNDS_STR,
                    internal_names.LON_STR, internal_names.LON_BOUNDS_STR,
                    internal_names.ZSURF_STR, internal_names.SFC_AREA_STR,
                    internal_names.LAND_MASK_STR, internal_names.PK_STR,
                    internal_names.BK_STR, internal_names.PHALF_STR,
                    internal_names.PFULL_STR, internal_names.PLEVEL_STR]
    _grid_attrs = OrderedDict([(key, internal_names.GRID_ATTRS[key])
                               for key in _grid_coords])

    def __str__(self):
        """String representation of the object."""
        return "<aospy.Calc instance: " + ', '.join(
            (self.name, self.proj.name, self.model.name, self.run.name)
        ) + ">"

    __repr__ = __str__

    def _dir_out(self):
        """Create string of the data directory to save individual .nc files."""
        return os.path.join(self.proj.direc_out, self.proj.name,
                            self.model.name, self.run.name, self.name)

    def _dir_tar_out(self):
        """Create string of the data directory to store a tar file."""
        return os.path.join(self.proj.tar_direc_out, self.proj.name,
                            self.model.name, self.run.name)

    def _file_name(self, dtype_out_time, extension='nc'):
        """Create the name of the aospy file."""
        if dtype_out_time is None:
            dtype_out_time = ''
        out_lbl = utils.io.data_out_label(self.intvl_out, dtype_out_time,
                                          dtype_vert=self.dtype_out_vert)
        in_lbl = utils.io.data_in_label(self.intvl_in, self.dtype_in_time,
                                        self.dtype_in_vert)
        start_year = utils.times.infer_year(self.start_date)
        end_year = utils.times.infer_year(self.end_date)
        yr_lbl = utils.io.yr_label((start_year, end_year))
        return '.'.join(
            [self.name, out_lbl, in_lbl, self.model.name,
             self.run.name, yr_lbl, extension]
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

    def __init__(self, proj=None, model=None, run=None, var=None,
                 date_range=None, region=None, intvl_in=None, intvl_out=None,
                 dtype_in_time=None, dtype_in_vert=None, dtype_out_time=None,
                 dtype_out_vert=None, level=None, time_offset=None):
        """Instantiate a Calc object.

        Parameters
        ----------
        proj : aospy.Proj object
            The project for this calculation.
        model : aospy.Model object
            The model for this calculation.
        run : aospy.Run object
            The run for this calculation.
        var : aospy.Var object
            The variable for this calculation.
        region : sequence of aospy.Region objects
            The region(s) over which any regional reductions will be performed.
        date_range : tuple of datetime.datetime objects
            The range of dates over which to perform calculations.
        intvl_in : {None, 'annual', 'monthly', 'daily', '6hr', '3hr'}, optional
            The time resolution of the input data.
        dtype_in_time : {None, 'inst', 'ts', 'av', 'av_ts'}, optional
            What the time axis of the input data represents:

            - 'inst' : Timeseries of instantaneous values
            - 'ts' : Timeseries of averages over the period of each time-index
            - 'av' : A single value averaged over a date range

        dtype_in_vert : {None, 'pressure', 'sigma'}, optional
            The vertical coordinate system used by the input data:

            - None : not defined vertically
            - 'pressure' : pressure coordinates
            - 'sigma' : hybrid sigma-pressure coordinates

        intvl_out : {'ann', season-string, month-integer}
            The sub-annual time interval over which to compute:

            - 'ann' : Annual mean
            - season-string : E.g. 'JJA' for June-July-August
            - month-integer : 1 for January, 2 for February, etc.

        dtype_out_time : tuple with elements being one or more of:

            - Gridpoint-by-gridpoint output:

              - 'av' : Gridpoint-by-gridpoint time-average
              - 'std' : Gridpoint-by-gridpoint temporal standard deviation
              - 'ts' : Gridpoint-by-gridpoint time-series

            - Averages over each region specified via `region`:

              - 'reg.av', 'reg.std', 'reg.ts' : analogous to 'av', 'std', 'ts'

        dtype_out_vert : {None, 'vert_av', 'vert_int'}, optional
            How to reduce the data vertically:

            - None : no vertical reduction (i.e. output is defined vertically)
            - 'vert_av' : mass-weighted vertical average
            - 'vert_int' : mass-weighted vertical integral

        time_offset : {None, dict}, optional
            How to offset input data in time to correct for metadata errors

            - None : no time offset applied
            - dict : e.g. ``{'hours': -3}`` to offset times by -3 hours
              See :py:meth:`aospy.utils.times.apply_time_offset`.

        """
        if run not in model.runs:
            raise AttributeError("Model '{0}' has no run '{1}'.  Calc object "
                                 "will not be generated.".format(model, run))
        self.proj = proj
        self.model = model

        self.run = run
        self.default_start_date = self.run.default_start_date
        self.default_end_date = self.run.default_end_date
        self.data_loader = self.run.data_loader

        self.var = var
        self.name = self.var.name
        self.domain = self.var.domain
        self.def_time = self.var.def_time
        self.def_vert = self.var.def_vert

        logging.debug(self._print_verbose('Initializing Calc '
                                          'instance:', self.__str__()))

        try:
            self.function = self.var.func
        except AttributeError:
            self.function = lambda x: x
        if getattr(self.var, 'variables', False):
            self.variables = self.var.variables
        else:
            self.variables = (self.var,)

        self.intvl_in = intvl_in
        self.intvl_out = intvl_out
        self.dtype_in_time = dtype_in_time
        self.dtype_in_vert = dtype_in_vert

        if isinstance(dtype_out_time, (list, tuple)):
            self.dtype_out_time = tuple(dtype_out_time)
        else:
            self.dtype_out_time = tuple([dtype_out_time])
        if not self.def_time:
            for reduction in self.dtype_out_time:
                if reduction in _TIME_DEFINED_REDUCTIONS:
                    msg = ("Var {0} has no time dimension "
                           "for the given time reduction "
                           "{1}".format(self.name, reduction))
                    raise ValueError(msg)

        self.level = level
        self.dtype_out_vert = dtype_out_vert
        self.region = region

        self.months = utils.times.month_indices(intvl_out)
        if date_range == 'default':
            self.start_date = utils.times.ensure_datetime(
                self.run.default_start_date)
            self.end_date = utils.times.ensure_datetime(
                self.run.default_end_date)
        else:
            self.start_date = utils.times.ensure_datetime(date_range[0])
            self.end_date = utils.times.ensure_datetime(date_range[-1])

        self.time_offset = time_offset
        self.data_loader_attrs = dict(
            domain=self.domain, intvl_in=self.intvl_in,
            dtype_in_vert=self.dtype_in_vert,
            dtype_in_time=self.dtype_in_time, intvl_out=self.intvl_out)

        self.model.set_grid_data()

        self.dir_out = self._dir_out()
        self.dir_tar_out = self._dir_tar_out()
        self.file_name = {d: self._file_name(d) for d in self.dtype_out_time}
        self.path_out = {d: self._path_out(d)
                         for d in self.dtype_out_time}
        self.path_tar_out = self._path_tar_out()

        self.data_out = {}

    def _to_desired_dates(self, arr):
        """Restrict the xarray DataArray or Dataset to the desired months."""
        times = utils.times.extract_months(
            arr[internal_names.TIME_STR], self.months
        )
        return arr.sel(time=times)

    def _add_grid_attributes(self, ds):
        """Add model grid attributes to a dataset"""
        for name_int, names_ext in self._grid_attrs.items():
            ds_coord_name = set(names_ext).intersection(set(ds.coords) |
                                                        set(ds.data_vars))
            model_attr = getattr(self.model, name_int, None)
            if ds_coord_name and (model_attr is not None):
                # Force coords to have desired name.
                ds = ds.rename({list(ds_coord_name)[0]: name_int})
                ds = ds.set_coords(name_int)
                if not np.array_equal(ds[name_int], model_attr):
                    if np.allclose(ds[name_int], model_attr):
                        msg = ("Values for '{0}' are nearly (but not exactly) "
                               "the same in the Run {1} and the Model {2}.  "
                               "Therefore replacing Run's values with the "
                               "model's.".format(name_int, self.run,
                                                 self.model))
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
            if (self.dtype_in_vert == 'pressure' and
                    internal_names.PLEVEL_STR in ds.coords):
                self.pressure = ds.level
        return ds

    def _get_input_data(self, var, start_date, end_date):
        """Get the data for a single variable over the desired date range."""
        logging.info(self._print_verbose("Getting input data:", var))

        if isinstance(var, (float, int)):
            return var
        else:
            cond_pfull = ((not hasattr(self, internal_names.PFULL_STR))
                          and var.def_vert and
                          self.dtype_in_vert == internal_names.ETA_STR)
            data = self.data_loader.recursively_compute_variable(
                var, start_date, end_date, self.time_offset, self.model,
                **self.data_loader_attrs)
            name = data.name
            data = self._add_grid_attributes(data.to_dataset(name=data.name))
            data = data[name]
            if cond_pfull:
                try:
                    self.pfull_coord = data[internal_names.PFULL_STR]
                except KeyError:
                    pass
            # Force all data to be at full pressure levels, not half levels.
            bool_to_pfull = (self.dtype_in_vert == internal_names.ETA_STR and
                             var.def_vert == internal_names.PHALF_STR)
            if bool_to_pfull:
                data = utils.vertcoord.to_pfull_from_phalf(data,
                                                           self.pfull_coord)
        if var.def_time:
            # Restrict to the desired dates within each year.
            if self.dtype_in_time != 'av':
                return self._to_desired_dates(data)
        else:
            return data

    def _get_all_data(self, start_date, end_date):
        """Get the needed data from all of the vars in the calculation."""
        return [self._get_input_data(var, start_date, end_date)
                for var in _replace_pressure(self.variables,
                                             self.dtype_in_vert)]

    def _local_ts(self, *data):
        """Perform the computation at each gridpoint and time index."""
        return self.function(*data).rename(self.name)

    def _compute(self, data):
        """Perform the calculation."""
        local_ts = self._local_ts(*data)
        dt = local_ts[internal_names.TIME_WEIGHTS_STR]
        # Convert dt to units of days to prevent overflow
        dt = dt / np.timedelta64(1, 'D')
        return local_ts, dt

    def _compute_full_ts(self, data):
        """Perform calculation and create yearly timeseries at each point."""
        # Get results at each desired timestep and spatial point.
        full_ts, dt = self._compute(data)
        # Vertically integrate.
        vert_types = ('vert_int', 'vert_av')
        if self.dtype_out_vert in vert_types and self.var.def_vert:
            dp = self._get_input_data(_DP_VARS[self.dtype_in_vert],
                                      self.start_date, self.end_date)
            full_ts = utils.vertcoord.int_dp_g(full_ts, dp)
            if self.dtype_out_vert == 'vert_av':
                ps = self._get_input_data(utils.vertcoord.ps,
                                          self.start_date, self.end_date)
                full_ts *= (GRAV_EARTH / ps)
        return full_ts, dt

    def _full_to_yearly_ts(self, arr, dt):
        """Average the full timeseries within each year."""
        time_defined = self.def_time and not ('av' in self.dtype_in_time)
        if time_defined:
            arr = utils.times.yearly_average(arr, dt)
        return arr

    def _time_reduce(self, arr, reduction):
        """Perform the specified time reduction on a local time-series."""
        if self.dtype_in_time == 'av' or not self.def_time:
            return arr
        reductions = {
            'ts': lambda xarr: xarr,
            'av': lambda xarr: xarr.mean(internal_names.YEAR_STR),
            'std': lambda xarr: xarr.std(internal_names.YEAR_STR),
            }
        try:
            return reductions[reduction](arr)
        except KeyError:
            raise ValueError("Specified time-reduction method '{}' is not "
                             "supported".format(reduction))

    def region_calcs(self, arr, func):
        """Perform a calculation for all regions."""
        # Get pressure values for data output on hybrid vertical coordinates.
        bool_pfull = (self.def_vert and self.dtype_in_vert ==
                      internal_names.ETA_STR and self.dtype_out_vert is False)
        if bool_pfull:
            pfull_data = self._get_input_data(_P_VARS[self.dtype_in_vert],
                                              self.start_date,
                                              self.end_date)
            pfull = self._full_to_yearly_ts(
                pfull_data, arr[internal_names.TIME_WEIGHTS_STR]
            ).rename('pressure')
        # Loop over the regions, performing the calculation.
        reg_dat = {}
        for reg in self.region:
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
        return xr.Dataset(reg_dat)

    def _apply_all_time_reductions(self, data):
        """Apply all requested time reductions to the data."""
        logging.info(self._print_verbose("Applying desired time-"
                                         "reduction methods."))
        reduc_specs = [r.split('.') for r in self.dtype_out_time]
        reduced = {}
        for reduc, specs in zip(self.dtype_out_time, reduc_specs):
            func = specs[-1]
            if 'reg' in specs:
                reduced.update({reduc: self.region_calcs(data, func)})
            else:
                reduced.update({reduc: self._time_reduce(data, func)})
        return OrderedDict(sorted(reduced.items(), key=lambda t: t[0]))

    def compute(self, write_to_tar=True):
        """Perform all desired calculations on the data and save externally."""
        data = self._get_all_data(self.start_date, self.end_date)
        logging.info('Computing timeseries for {0} -- '
                     '{1}.'.format(self.start_date, self.end_date))
        full, full_dt = self._compute_full_ts(data)
        full_out = self._full_to_yearly_ts(full, full_dt)
        reduced = self._apply_all_time_reductions(full_out)
        logging.info("Writing desired gridded outputs to disk.")
        for dtype_time, data in reduced.items():
            data = _add_metadata_as_attrs(data, self.var.units,
                                          self.var.description,
                                          self.dtype_out_vert)
            self.save(data, dtype_time, dtype_out_vert=self.dtype_out_vert,
                      save_files=True, write_to_tar=write_to_tar)
        return self

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
            reg_data.update(data)
            data_out = reg_data
        else:
            data_out = data
        if isinstance(data_out, xr.DataArray):
            data_out = xr.Dataset({self.name: data_out})
        data_out.to_netcdf(path, engine='netcdf4', format='NETCDF3_64BIT')

    def _write_to_tar(self, dtype_out_time):
        """Add the data to the tar file in tar_out_direc."""
        # When submitted in parallel and the directory does not exist yet
        # multiple processes may try to create a new directory; this leads
        # to an OSError for all processes that tried to make the
        # directory, but were later than the first.
        try:
            os.makedirs(self.dir_tar_out)
        except OSError:
            pass
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
                retcode = subprocess.call([
                    "tar", "--delete", "--file={}".format(self.path_tar_out),
                    self.file_name[dtype_out_time]
                ])
                if retcode:
                    msg = ("The 'tar' command to save your aospy output "
                           "exited with an error.  Most likely, this is due "
                           "to using an old version of 'tar' (especially if "
                           "you are on a Mac).  Consider installing a newer "
                           "version of 'tar' or disabling tar output by "
                           "setting `write_to_tar=False` in the "
                           "`calc_exec_options` argument of "
                           "`submit_mult_calcs`.")
                    logging.warn(msg)
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
             save_files=True, write_to_tar=False):
        """Save aospy data to data_out attr and to an external file."""
        self._update_data_out(data, dtype_out_time)
        if save_files:
            self._save_files(data, dtype_out_time)
        if write_to_tar and self.proj.tar_direc_out:
            self._write_to_tar(dtype_out_time)
        logging.info('\t{}'.format(self.path_out[dtype_out_time]))

    def _load_from_disk(self, dtype_out_time, dtype_out_vert=False,
                        region=False):
        """Load aospy data saved as netcdf files on the file system."""
        ds = xr.open_dataset(self.path_out[dtype_out_time])
        if region:
            arr = ds[region.name]
            # Use region-specific pressure values if available.
            if (self.dtype_in_vert == internal_names.ETA_STR
                    and not dtype_out_vert):
                reg_pfull_str = region.name + '_pressure'
                arr = arr.drop([r for r in arr.coords.iterkeys()
                                if r not in (internal_names.PFULL_STR,
                                             reg_pfull_str)])
                # Rename pfull to pfull_ref always.
                arr = arr.rename({internal_names.PFULL_STR:
                                  internal_names.PFULL_STR + '_ref'})
                # Rename region_pfull to pfull if its there.
                if hasattr(arr, reg_pfull_str):
                    return arr.rename({reg_pfull_str:
                                       internal_names.PFULL_STR})
                return arr
            return arr
        return ds[self.name]

    def _load_from_tar(self, dtype_out_time, dtype_out_vert=False):
        """Load data save in tarball form on the file system."""
        path = os.path.join(self.dir_tar_out, 'data.tar')
        utils.io.dmget([path])
        with tarfile.open(path, 'r') as data_tar:
            ds = xr.open_dataset(
                data_tar.extractfile(self.file_name[dtype_out_time])
            )
            return ds[self.name]

    def load(self, dtype_out_time, dtype_out_vert=False, region=False,
             plot_units=False, mask_unphysical=False):
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
        # Apply desired plotting/cleanup methods.
        if mask_unphysical:
            data = self.var.mask_unphysical(data)
        if plot_units:
            data = self.var.to_plot_units(data, dtype_vert=dtype_out_vert)
        return data


def _add_metadata_as_attrs(data, units, description, dtype_out_vert):
    """Add metadata attributes to Dataset or DataArray"""
    if isinstance(data, xr.DataArray):
        return _add_metadata_as_attrs_da(data, units, description,
                                         dtype_out_vert)
    else:
        for name, arr in data.data_vars.items():
            _add_metadata_as_attrs_da(arr, units, description,
                                      dtype_out_vert)
        return data


def _add_metadata_as_attrs_da(data, units, description, dtype_out_vert):
    """Add metadata attributes to DataArray"""
    if dtype_out_vert == 'vert_int':
        if units != '':
            units = '(vertical integral of {0}): {0} kg m^-2)'.format(units)
        else:
            units = '(vertical integral of quantity with unspecified units)'
    data.attrs['units'] = units
    data.attrs['description'] = description
    return data
