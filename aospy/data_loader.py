"""aospy DataLoader objects"""
import logging
import os

import numpy as np
import xarray as xr

from .internal_names import (ETA_STR, GRID_ATTRS, TIME_STR, TIME_BOUNDS_STR,
                             TIME_VAR_STRS)
from .utils import times, io


def _preprocess_and_rename_grid_attrs(func, **kwargs):
    """Call a custom preprocessing method first then rename grid attrs

    This wrapper is needed to generate a single function to pass to the
    ``preprocesss`` of xr.open_mfdataset.  It makes sure that the
    user-specified preprocess function is called on the loaded Dataset before
    aospy's is applied.  An example for why this might be needed is output from
    the WRF model; one needs to add a CF-compliant units attribute to the time
    coordinate of all input files, because it is not present by default.

    Parameters
    ----------
    func : function
       An arbitrary function to call before calling
       ``grid_attrs_to_aospy_names`` in ``_load_data_from_disk``.  Must take
       an xr.Dataset as an argument as well as ``**kwargs``.

    Returns
    -------
    function
        A function that calls the provided function ``func`` on the Dataset
        before calling ``grid_attrs_to_aospy_names``; this is meant to be
        passed as a ``preprocess`` argument to ``xr.open_mfdataset``.
    """

    def func_wrapper(ds):
        return grid_attrs_to_aospy_names(func(ds, **kwargs))
    return func_wrapper


def grid_attrs_to_aospy_names(data):
    """Rename grid attributes to be consistent with aospy conventions.

    Search all of the dataset's coords and dims looking for matches to known
    grid attribute names; any that are found subsequently get renamed to the
    aospy name as specified in ``aospy.internal_names.GRID_ATTRS``.

    Also forces any renamed grid attribute that is saved as a dim without a
    coord to have a coord, which facilitates subsequent slicing/subsetting.

    This function does not compare to Model coordinates or add missing
    coordinates from Model objects.

    Parameters
    ----------
    data : xr.Dataset

    Returns
    -------
    xr.Dataset
        Data returned with coordinates consistent with aospy
        conventions
    """
    for name_int, names_ext in GRID_ATTRS.items():
        dims_and_vars = set(data.variables).union(set(data.dims))
        data_coord_name = set(names_ext).intersection(dims_and_vars)
        if data_coord_name:
            data = data.rename({data_coord_name.pop(): name_int})
    return set_grid_attrs_as_coords(data, set_time_vars=False)


def set_grid_attrs_as_coords(ds, set_time_vars=True):
    """Set available grid attributes as coordinates in a given Dataset.

    Grid attributes are assumed to have their internal aospy names. Grid
    attributes are set as coordinates, such that they are carried by all
    selected DataArrays with overlapping index dimensions.

    Parameters
    ----------
    ds : Dataset
        Input data

    Returns
    -------
    Dataset
        Dataset with grid attributes set as coordinates
    """
    if set_time_vars:
        int_names = GRID_ATTRS.keys()
    else:
        int_names = (set(GRID_ATTRS.keys()) -
                     set(TIME_VAR_STRS))
    grid_attrs_in_ds = set(int_names).intersection(set(ds.coords) |
                                                   set(ds.data_vars))
    ds.set_coords(grid_attrs_in_ds, inplace=True)
    return ds


def _maybe_cast_to_float64(da):
    """Cast DataArrays to np.float64 if they are of type np.float32

    Parameters
    ----------
    da : xr.DataArray
        Input DataArray

    Returns
    -------
    DataArray
    """
    if da.dtype == np.float32:
        logging.warning('Datapoints were stored using the np.float32 datatype.'
                        'For accurate reduction operations using bottleneck, '
                        'datapoints are being cast to the np.float64 datatype.'
                        ' For more information see: https://github.com/pydata/'
                        'xarray/issues/1346')
        return da.astype(np.float64)
    else:
        return da


def _sel_var(ds, var, upcast_float32=True):
    """Select the specified variable by trying all possible alternative names.

    Parameters
    ----------
    ds : Dataset
        Dataset possibly containing var
    var : aospy.Var
        Variable to find data for
    upcast_float32 : bool (default True)
        Whether to cast a float32 DataArray up to float64

    Returns
    -------
    DataArray

    Raises
    ------
    KeyError
        If the variable is not in the Dataset
    """
    for name in var.names:
        try:
            da = ds[name].rename(var.name)
            if upcast_float32:
                return _maybe_cast_to_float64(da)
            else:
                return da
        except KeyError:
            pass
    msg = '{0} not found among names: {1} in\n{2}'.format(var, var.names, ds)
    raise LookupError(msg)


def _prep_time_data(ds):
    """Prepare time coord. information in Dataset for use in aospy.

    1. Edit units attribute of time variable if it contains a Timestamp invalid
       date
    2. If the Dataset contains a time bounds coordinate, add attributes
       representing the true beginning and end dates of the time interval used
       to construct the Dataset
    3. If the Dataset contains a time bounds coordinate, overwrite the time
       coordinate values with the averages of the time bounds at each timestep
    4. Decode the times into np.datetime64 objects for time indexing

    Parameters
    ----------
    ds : Dataset
        Pre-processed Dataset with time coordinate renamed to
        internal_names.TIME_STR

    Returns
    -------
    Dataset, int, int
        The processed Dataset and minimum and maximum years in the loaded data
    """
    ds = times.ensure_time_as_dim(ds)
    ds, min_year, max_year = times.numpy_datetime_workaround_encode_cf(ds)
    if TIME_BOUNDS_STR in ds:
        ds = times.ensure_time_avg_has_cf_metadata(ds)
        ds[TIME_STR] = times.average_time_bounds(ds)
    else:
        logging.warning("dt array not found.  Assuming equally spaced "
                        "values in time, even though this may not be "
                        "the case")
        ds = times.add_uniform_time_weights(ds)
    ds = xr.decode_cf(ds, decode_times=True, decode_coords=False,
                      mask_and_scale=True)
    return ds, min_year, max_year


def _load_data_from_disk(file_set, preprocess_func=lambda ds: ds, **kwargs):
    """Load a Dataset from a list or glob-string of files.

    Datasets from files are concatenated along time,
    and all grid attributes are renamed to their aospy internal names.

    Parameters
    ----------
    file_set : list or str
        List of paths to files or glob-string
    preprocess_func : function (optional)
        Custom function to call before applying any aospy logic
        to the loaded dataset

    Returns
    -------
    Dataset
    """
    apply_preload_user_commands(file_set)
    func = _preprocess_and_rename_grid_attrs(preprocess_func, **kwargs)
    return xr.open_mfdataset(file_set, preprocess=func, concat_dim=TIME_STR,
                             decode_times=False, decode_coords=False,
                             mask_and_scale=True)


def apply_preload_user_commands(file_set, cmd=io.dmget):
    """Call desired functions on file list before loading.

    For example, on the NOAA Geophysical Fluid Dynamics Laboratory
    computational cluster, data that is saved on their tape archive
    must be accessed via a `dmget` (or `hsmget`) command before being used.
    """
    if cmd is not None:
        cmd(file_set)


class DataLoader(object):
    """A fundamental DataLoader object"""
    def load_variable(self, var=None, start_date=None, end_date=None,
                      time_offset=None, **DataAttrs):
        """Load a DataArray for requested variable and time range.

        Automatically renames all grid attributes to match aospy conventions.

        Parameters
        ----------
        var : Var
            aospy Var object
        start_date : datetime.datetime
            start date for interval
        end_date : datetime.datetime
            end date for interval
        time_offset : dict
            Option to add a time offset to the time coordinate to correct for
            incorrect metadata.
        **DataAttrs
            Attributes needed to identify a unique set of files to load from

        Returns
        -------
        da : DataArray
             DataArray for the specified variable, date range, and interval in
        """
        file_set = self._generate_file_set(var=var, start_date=start_date,
                                           end_date=end_date, **DataAttrs)
        ds = _load_data_from_disk(file_set, self.preprocess_func,
                                  start_date=start_date, end_date=end_date,
                                  time_offset=time_offset, **DataAttrs)
        ds, min_year, max_year = _prep_time_data(ds)
        ds = set_grid_attrs_as_coords(ds)
        da = _sel_var(ds, var, self.upcast_float32)
        da = self._maybe_apply_time_shift(da, time_offset, **DataAttrs)

        start_date_xarray = times.numpy_datetime_range_workaround(
            start_date, min_year, max_year)
        end_date_xarray = start_date_xarray + (end_date - start_date)
        return times.sel_time(da, np.datetime64(start_date_xarray),
                              np.datetime64(end_date_xarray)).load()

    @staticmethod
    def _maybe_apply_time_shift(da, time_offset=None, **DataAttrs):
        """Apply specified time shift to DataArray"""
        if time_offset is not None:
            time = times.apply_time_offset(da[TIME_STR], **time_offset)
            da[TIME_STR] = time
        return da

    def _generate_file_set(self, var=None, start_date=None, end_date=None,
                           domain=None, intvl_in=None, dtype_in_vert=None,
                           dtype_in_time=None, intvl_out=None):
        raise NotImplementedError(
            'All DataLoaders require a _generate_file_set method')


class DictDataLoader(DataLoader):
    """A DataLoader that uses a dict mapping lists of files to string tags

    This is the simplest DataLoader; it is useful for instance if one is
    dealing with raw model history files, which tend to group all variables
    of a single output interval into single filesets. The
    intvl_in parameter is a string description of the time frequency of the
    data one is referencing (e.g. 'monthly', 'daily', '3-hourly').  In
    principle, one can give it any string value.

    Parameters
    ----------
    file_map : dict
        A dict mapping an input interval to a list of files
    upcast_float32 : bool (default True)
        Whether to cast loaded DataArrays with the float32 datatype to float64
        before doing calculations
    preprocess_func : function (optional)
        A function to apply to every Dataset before processing in aospy.  Must
        take a Dataset and ``**kwargs`` as its two arguments.

    Examples
    --------
    Case of two sets of files, one with monthly average output, and one with
    3-hourly output.

    >>> file_map = {'monthly': '000[4-6]0101.atmos_month.nc',
    ...             '3hr': '000[4-6]0101.atmos_8xday.nc'}
    >>> data_loader = DictDataLoader(file_map)

    If one wanted to correct a CF-incompliant units attribute on each Dataset
    read in, which depended on the ``intvl_in`` of the fileset one could
    define a ``preprocess_func`` which took into account the ``intvl_in``
    keyword argument.

    >>> def preprocess(ds, **kwargs):
    ...     if kwargs['intvl_in'] == 'monthly':
    ...         ds['time'].attrs['units'] = 'days since 0001-01-0000'
    ...     if kwargs['intvl_in'] == '3hr':
    ...         ds['time'].attrs['units'] = 'hours since 0001-01-0000'
    ...     return ds
    >>> data_loader = DictDataLoader(file_map, preprocess)
    """
    def __init__(self, file_map=None, upcast_float32=True,
                 preprocess_func=lambda ds, **kwargs: ds):
        """Create a new DictDataLoader"""
        self.file_map = file_map
        self.upcast_float32 = upcast_float32
        self.preprocess_func = preprocess_func

    def _generate_file_set(self, var=None, start_date=None, end_date=None,
                           domain=None, intvl_in=None, dtype_in_vert=None,
                           dtype_in_time=None, intvl_out=None):
        """Returns the file_set for the given interval in."""
        try:
            return self.file_map[intvl_in]
        except KeyError:
            raise KeyError('File set does not exist for the specified'
                           ' intvl_in {0}'.format(intvl_in))


class NestedDictDataLoader(DataLoader):
    """DataLoader that uses a nested dictionary mapping to load files

    This is the most flexible existing type of DataLoader; it allows for the
    specification of different sets of files for different variables.  The
    intvl_in parameter is a string description of the time frequency of the
    data one is referencing (e.g. 'monthly', 'daily', '3-hourly').  In
    principle, one can give it any string value.  The variable name
    can be any variable name in your aospy object library (including
    alternative names).

    Parameters
    ----------
    file_map : dict
        A dict mapping intvl_in to dictionaries mapping Var
        objects to lists of files
    upcast_float32 : bool (default True)
        Whether to cast loaded DataArrays with the float32 datatype to float64
        before doing calculations
    preprocess_func : function (optional)
        A function to apply to every Dataset before processing in aospy.  Must
        take a Dataset and ``**kwargs`` as its two arguments.

    Examples
    --------
    Case of a set of monthly average files for large scale precipitation,
    and another monthly average set of files for convective precipitation.

    >>> file_map = {'monthly': {'precl': '000[4-6]0101.precl.nc',
    ...                         'precc': '000[4-6]0101.precc.nc'}}
    >>> data_loader = NestedDictDataLoader(file_map)

    See :py:class:`aospy.data_loader.DictDataLoader` for an example of a
    possible function to pass as a ``preprocess_func``.
    """
    def __init__(self, file_map=None, upcast_float32=True,
                 preprocess_func=lambda ds, **kwargs: ds):
        """Create a new NestedDictDataLoader"""
        self.file_map = file_map
        self.upcast_float32 = upcast_float32
        self.preprocess_func = preprocess_func

    def _generate_file_set(self, var=None, start_date=None, end_date=None,
                           domain=None, intvl_in=None, dtype_in_vert=None,
                           dtype_in_time=None, intvl_out=None):
        for name in var.names:
            try:
                return self.file_map[intvl_in][name]
            except KeyError:
                pass
        raise KeyError('Files for the var {0} cannot be found in for the '
                       'intvl_in {1} in this'
                       ' OneDirDataLoader'.format(var, intvl_in))


class GFDLDataLoader(DataLoader):
    """DataLoader for NOAA GFDL model output

    This is an example of a domain-specific custom DataLoader, designed
    specifically for finding files output by the Geophysical Fluid Dynamics
    Laboratory's model history file post-processing tools.

    Parameters
    ----------
    template : GFDLDataLoader
        Optional argument to specify a base GFDLDataLoader to inherit
        parameters from
    data_direc : str
        Root directory of data files
    data_dur : int
        Number of years included per post-processed file
    data_start_date : datetime.datetime
        Start date of data files
    data_end_date : datetime.datetime
        End date of data files
    upcast_float32 : bool (default True)
        Whether to cast loaded DataArrays with the float32 datatype to float64
        before doing calculations
    preprocess_func : function (optional)
        A function to apply to every Dataset before processing in aospy.  Must
        take a Dataset and ``**kwargs`` as its two arguments.

    Examples
    --------
    Case without a template to start from.

    >>> base = GFDLDataLoader(data_direc='/archive/control/pp', data_dur=5,
    ...                       data_start_date=datetime(2000, 1, 1),
    ...                       data_end_date=datetime(2010, 12, 31))

    Case with a starting template.

    >>> data_loader = GFDLDataLoader(base, data_direc='/archive/2xCO2/pp')

    See :py:class:`aospy.data_loader.DictDataLoader` for an example of a
    possible function to pass as a ``preprocess_func``.
    """
    def __init__(self, template=None, data_direc=None, data_dur=None,
                 data_start_date=None, data_end_date=None,
                 upcast_float32=None,
                 preprocess_func=lambda ds, **kwargs: ds):
        """Create a new GFDLDataLoader"""
        attrs = ['data_direc', 'data_dur', 'data_start_date', 'data_end_date',
                 'preprocess_func']
        if template:
            for attr in attrs:
                setattr(self, attr, getattr(template, attr))

            if data_direc is not None:
                self.data_direc = data_direc
            else:
                self.data_direc = template.data_direc
            if data_dur is not None:
                self.data_dur = data_dur
            else:
                self.data_dur = template.data_dur
            if data_start_date is not None:
                self.data_start_date = data_start_date
            else:
                self.data_start_date = template.data_start_date
            if data_end_date is not None:
                self.data_end_date = data_end_date
            else:
                self.data_end_date = template.data_end_date
            if preprocess_func is not None:
                self.preprocess_func = preprocess_func
            else:
                self.preprocess_func = template.preprocess_func
            if upcast_float32 is not None:
                self.upcast_float32 = upcast_float32
            else:
                self.upcast_float32 = template.upcast_float32
        else:
            self.data_direc = data_direc
            self.data_dur = data_dur
            self.data_start_date = data_start_date
            self.data_end_date = data_end_date
            self.preprocess_func = preprocess_func
            if upcast_float32 is None:
                self.upcast_float32 = True
            else:
                self.upcast_float32 = upcast_float32

    @staticmethod
    def _maybe_apply_time_shift(da, time_offset=None, **DataAttrs):
        """Correct off-by-one error in GFDL instantaneous model data.

        Instantaneous data that is outputted by GFDL models is generally off by
        one timestep.  For example, a netCDF file that is supposed to
        correspond to 6 hourly data for the month of January, will have its
        last time value be in February.
        """
        if time_offset is not None:
            time = times.apply_time_offset(da[TIME_STR], **time_offset)
            da[TIME_STR] = time
        else:
            if DataAttrs['dtype_in_time'] == 'inst':
                if DataAttrs['intvl_in'].endswith('hr'):
                    offset = -1 * int(DataAttrs['intvl_in'][0])
                else:
                    offset = 0
                time = times.apply_time_offset(da[TIME_STR], hours=offset)
                da[TIME_STR] = time
        return da

    def _generate_file_set(self, var=None, start_date=None, end_date=None,
                           domain=None, intvl_in=None, dtype_in_vert=None,
                           dtype_in_time=None, intvl_out=None):
        for name in var.names:
            file_set = self._input_data_paths_gfdl(
                name, start_date, end_date, domain, intvl_in, dtype_in_vert,
                dtype_in_time, intvl_out)
            if all([os.path.isfile(filename) for filename in file_set]):
                return file_set
        raise IOError('Files for the var {0} cannot be located'
                      'using GFDL post-processing conventions'.format(var))

    def _input_data_paths_gfdl(self, name, start_date, end_date, domain,
                               intvl_in, dtype_in_vert, dtype_in_time,
                               intvl_out):
        dtype_lbl = dtype_in_time
        if intvl_in == 'daily':
            domain += '_daily'
        if dtype_in_vert == ETA_STR and name != 'ps':
            domain += '_level'
        if dtype_in_time == 'inst':
            domain += '_inst'
            dtype_lbl = 'ts'
        if 'monthly_from_' in dtype_in_time:
            dtype = dtype_in_time.replace('monthly_from_', '')
            dtype_lbl = dtype
        else:
            dtype = dtype_in_time
        dur_str = str(self.data_dur) + 'yr'
        if dtype_in_time == 'av':
            subdir = intvl_in + '_' + dur_str
        else:
            subdir = os.path.join(intvl_in, dur_str)
        direc = os.path.join(self.data_direc, domain, dtype_lbl, subdir)
        files = [os.path.join(direc, io.data_name_gfdl(
                    name, domain, dtype, intvl_in, year, intvl_out,
                    self.data_start_date.year, self.data_dur))
                 for year in range(start_date.year, end_date.year + 1)]
        files = list(set(files))
        files.sort()
        return files
