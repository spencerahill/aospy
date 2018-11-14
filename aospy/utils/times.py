"""Utility functions for handling times, dates, etc."""
import datetime
import logging
import re
import warnings

import cftime
import numpy as np
import pandas as pd
import xarray as xr

from ..internal_names import (
    BOUNDS_STR, RAW_END_DATE_STR, RAW_START_DATE_STR,
    SUBSET_END_DATE_STR, SUBSET_START_DATE_STR, TIME_BOUNDS_STR, TIME_STR,
    TIME_WEIGHTS_STR
)


def apply_time_offset(time, years=0, months=0, days=0, hours=0):
    """Apply a specified offset to the given time array.

    This is useful for GFDL model output of instantaneous values.  For example,
    3 hourly data postprocessed to netCDF files spanning 1 year each will
    actually have time values that are offset by 3 hours, such that the first
    value is for 1 Jan 03:00 and the last value is 1 Jan 00:00 of the
    subsequent year.  This causes problems in xarray, e.g. when trying to group
    by month.  It is resolved by manually subtracting off those three hours,
    such that the dates span from 1 Jan 00:00 to 31 Dec 21:00 as desired.

    Parameters
    ----------
    time : xarray.DataArray representing a timeseries
    years, months, days, hours : int, optional
        The number of years, months, days, and hours, respectively, to offset
        the time array by.  Positive values move the times later.

    Returns
    -------
    pandas.DatetimeIndex

    Examples
    --------
    Case of a length-1 input time array:

    >>> times = xr.DataArray(datetime.datetime(1899, 12, 31, 21))
    >>> apply_time_offset(times)
    Timestamp('1900-01-01 00:00:00')

    Case of input time array with length greater than one:

    >>> times = xr.DataArray([datetime.datetime(1899, 12, 31, 21),
    ...                       datetime.datetime(1899, 1, 31, 21)])
    >>> apply_time_offset(times) # doctest: +NORMALIZE_WHITESPACE
    DatetimeIndex(['1900-01-01', '1899-02-01'], dtype='datetime64[ns]',
                  freq=None)
    """
    return (pd.to_datetime(time.values) +
            pd.tseries.offsets.DateOffset(years=years, months=months,
                                          days=days, hours=hours))


def average_time_bounds(ds):
    """Return the average of each set of time bounds in the Dataset.

    Useful for creating a new time array to replace the Dataset's native time
    array, in the case that the latter matches either the start or end bounds.
    This can cause errors in grouping (akin to an off-by-one error) if the
    timesteps span e.g. one full month each.  Note that the Dataset's times
    must not have already undergone "CF decoding", wherein they are converted
    from floats using the 'units' attribute into datetime objects.

    Parameters
    ----------
    ds : xarray.Dataset
        A Dataset containing a time bounds array with name matching
        internal_names.TIME_BOUNDS_STR.  This time bounds array must have two
        dimensions, one of which's coordinates is the Dataset's time array, and
        the other is length-2.

    Returns
    -------
    xarray.DataArray
        The mean of the start and end times of each timestep in the original
        Dataset.

    Raises
    ------
    ValueError
        If the time bounds array doesn't match the shape specified above.

    """
    bounds = ds[TIME_BOUNDS_STR]
    new_times = bounds.mean(dim=BOUNDS_STR, keep_attrs=True)
    new_times = new_times.drop(TIME_STR).rename(TIME_STR)
    new_times[TIME_STR] = new_times
    return new_times


def monthly_mean_ts(arr):
    """Convert a sub-monthly time-series into one of monthly means.

    Also drops any months with no data in the original DataArray.

    Parameters
    ----------
    arr : xarray.DataArray
        Timeseries of sub-monthly temporal resolution data

    Returns
    -------
    xarray.DataArray
        Array resampled to comprise monthly means

    See Also
    --------
    monthly_mean_at_each_ind : Copy monthly means to each submonthly time

    """
    return arr.resample(**{TIME_STR: '1M'}).mean(TIME_STR).dropna(TIME_STR)


def monthly_mean_at_each_ind(monthly_means, sub_monthly_timeseries):
    """Copy monthly mean over each time index in that month.

    Parameters
    ----------
    monthly_means : xarray.DataArray
        array of monthly means
    sub_monthly_timeseries : xarray.DataArray
        array of a timeseries at sub-monthly time resolution

    Returns
    -------
    xarray.DataArray with eath monthly mean value from `monthly_means` repeated
    at each time within that month from `sub_monthly_timeseries`

    See Also
    --------
    monthly_mean_ts : Create timeseries of monthly mean values
    """
    time = monthly_means[TIME_STR]
    start = time.indexes[TIME_STR][0].replace(day=1, hour=0)
    end = time.indexes[TIME_STR][-1]
    new_indices = pd.DatetimeIndex(start=start, end=end, freq='MS')
    arr_new = monthly_means.reindex(time=new_indices, method='backfill')
    return arr_new.reindex_like(sub_monthly_timeseries, method='pad')


def yearly_average(arr, dt):
    """Average a sub-yearly time-series over each year.

    Resulting timeseries comprises one value for each year in which the
    original array had valid data.  Accounts for (i.e. ignores) masked values
    in original data when computing the annual averages.

    Parameters
    ----------
    arr : xarray.DataArray
        The array to be averaged
    dt : xarray.DataArray
        Array of the duration of each timestep

    Returns
    -------
    xarray.DataArray
        Has the same shape and mask as the original ``arr``, except for the
        time dimension, which is truncated to one value for each year that
        ``arr`` spanned

    """
    assert_matching_time_coord(arr, dt)
    yr_str = TIME_STR + '.year'
    # Retain original data's mask.
    dt = dt.where(np.isfinite(arr))
    return ((arr*dt).groupby(yr_str).sum(TIME_STR) /
            dt.groupby(yr_str).sum(TIME_STR))


def ensure_datetime(obj):
    """Return the object if it is a datetime-like object

    Parameters
    ----------
    obj : Object to be tested.

    Returns
    -------
    The original object if it is a datetime-like object

    Raises
    ------
    TypeError if `obj` is not datetime-like
    """
    _VALID_TYPES = (str, datetime.datetime, cftime.datetime,
                    np.datetime64)
    if isinstance(obj, _VALID_TYPES):
        return obj
    raise TypeError("datetime-like object required.  "
                    "Type given: {}".format(type(obj)))


def datetime_or_default(date, default):
    """Return a datetime-like object or a default.

    Parameters
    ----------
    date : `None` or datetime-like object or str
    default : The value to return if `date` is `None`

    Returns
    -------
    `default` if `date` is `None`, otherwise returns the result of
    `utils.times.ensure_datetime(date)`

    """
    if date is None:
        return default
    else:
        return ensure_datetime(date)


def month_indices(months):
    """Convert string labels for months to integer indices.

    Parameters
    ----------
     months : str, int
         If int, number of the desired month, where January=1, February=2,
         etc.  If str, must match either 'ann' or some subset of
         'jfmamjjasond'.  If 'ann', use all months.  Otherwise, use the
         specified months.

    Returns
    -------
    np.ndarray of integers corresponding to desired month indices

    Raises
    ------
    TypeError : If `months` is not an int or str

    See also
    --------
    _month_conditional
    """
    if not isinstance(months, (int, str)):
        raise TypeError("`months` must be of type int or str: "
                        "type(months) == {}".format(type(months)))
    if isinstance(months, int):
        return [months]
    if months.lower() == 'ann':
        return np.arange(1, 13)
    first_letter = 'jfmamjjasond' * 2
    # Python indexing starts at 0; month indices start at 1 for January.
    count = first_letter.count(months)
    if (count == 0) or (count > 2):
        message = ("The user must provide a unique pattern of consecutive "
                   "first letters of months within '{}'. The provided "
                   "string '{}' does not comply."
                   "  For individual months use integers."
                   "".format(first_letter, months))
        raise ValueError(message)
    st_ind = first_letter.find(months.lower())
    return np.arange(st_ind, st_ind + len(months)) % 12 + 1


def _month_conditional(time, months):
    """Create a conditional statement for selecting data in a DataArray.

    Parameters
    ----------
    time : xarray.DataArray
         Array of times for which to subsample for specific months.
    months : int, str, or xarray.DataArray of times
        If int or str, passed to `month_indices`
    Returns
    -------
    Array of bools specifying which months to keep

    See Also
    --------
    month_indices
    """
    if isinstance(months, (int, str)):
        months_array = month_indices(months)
    else:
        months_array = months
    cond = False
    for month in months_array:
        cond |= (time['{}.month'.format(TIME_STR)] == month)
    return cond


def extract_months(time, months):
    """Extract times within specified months of the year.

    Parameters
    ----------
    time : xarray.DataArray
         Array of times that can be represented by numpy.datetime64 objects
         (i.e. the year is between 1678 and 2262).
    months : Desired months of the year to include

    Returns
    -------
    xarray.DataArray of the desired times
    """
    inds = _month_conditional(time, months)
    return time.sel(time=inds)


def ensure_time_avg_has_cf_metadata(ds):
    """Add time interval length and bounds coordinates for time avg data.

    If the Dataset or DataArray contains time average data, enforce
    that there are coordinates that track the lower and upper bounds of
    the time intervals, and that there is a coordinate that tracks the
    amount of time per time average interval.

    CF conventions require that a quantity stored as time averages
    over time intervals must have time and time_bounds coordinates [1]_.
    aospy further requires AVERAGE_DT for time average data, for accurate
    time-weighted averages, which can be inferred from the CF-required
    time_bounds coordinate if needed.  This step should be done
    prior to decoding CF metadata with xarray to ensure proper
    computed timedeltas for different calendar types.

    .. [1] http://cfconventions.org/cf-conventions/v1.6.0/cf-conventions.html#_data_representative_of_cells

    Parameters
    ----------
    ds : Dataset or DataArray
        Input data

    Returns
    -------
    Dataset or DataArray
        Time average metadata attributes added if needed.
    """  # noqa: E501
    if TIME_WEIGHTS_STR not in ds:
        time_weights = ds[TIME_BOUNDS_STR].diff(BOUNDS_STR)
        time_weights = time_weights.rename(TIME_WEIGHTS_STR).squeeze()
        if BOUNDS_STR in time_weights.coords:
            time_weights = time_weights.drop(BOUNDS_STR)
        ds[TIME_WEIGHTS_STR] = time_weights

    raw_start_date = ds[TIME_BOUNDS_STR].isel(**{TIME_STR: 0, BOUNDS_STR: 0})
    ds[RAW_START_DATE_STR] = raw_start_date.reset_coords(drop=True)
    raw_end_date = ds[TIME_BOUNDS_STR].isel(**{TIME_STR: -1, BOUNDS_STR: 1})
    ds[RAW_END_DATE_STR] = raw_end_date.reset_coords(drop=True)

    for coord in [TIME_BOUNDS_STR, RAW_START_DATE_STR, RAW_END_DATE_STR]:
        ds[coord].attrs['units'] = ds[TIME_STR].attrs['units']
        if 'calendar' in ds[TIME_STR].attrs:
            ds[coord].attrs['calendar'] = ds[TIME_STR].attrs['calendar']

    unit_interval = ds[TIME_STR].attrs['units'].split('since')[0].strip()
    ds[TIME_WEIGHTS_STR].attrs['units'] = unit_interval
    return ds


def add_uniform_time_weights(ds):
    """Append uniform time weights to a Dataset.

    All DataArrays with a time coordinate require a time weights coordinate.
    For Datasets read in without a time bounds coordinate or explicit
    time weights built in, aospy adds uniform time weights at each point
    in the time coordinate.

    Parameters
    ----------
    ds : Dataset
        Input data

    Returns
    -------
    Dataset
    """
    time = ds[TIME_STR]
    unit_interval = time.attrs['units'].split('since')[0].strip()
    time_weights = xr.ones_like(time)
    time_weights.attrs['units'] = unit_interval
    del time_weights.attrs['calendar']
    ds[TIME_WEIGHTS_STR] = time_weights
    return ds


def _assert_has_data_for_time(da, start_date, end_date):
    """Check to make sure data is in Dataset for the given time range.

    Parameters
    ----------
    da : DataArray
         DataArray with a time variable
    start_date : datetime-like object or str
         start date
    end_date : datetime-like object or str
         end date

    Raises
    ------
    AssertionError
         If the time range is not within the time range of the DataArray

    """
    if isinstance(start_date, str) and isinstance(end_date, str):
        logging.warning(
            'When using strings to specify start and end dates, the check '
            'to determine if data exists for the full extent of the desired '
            'interval is not implemented.  Therefore it is possible that '
            'you are doing a calculation for a lesser interval than you '
            'specified.  If you would like this check to occur, use explicit '
            'datetime-like objects for bounds instead.')
        return

    if RAW_START_DATE_STR in da.coords:
        with warnings.catch_warnings(record=True):
            da_start = da[RAW_START_DATE_STR].values
            da_end = da[RAW_END_DATE_STR].values
    else:
        times = da.time.isel(**{TIME_STR: [0, -1]})
        da_start, da_end = times.values

    message = ('Data does not exist for requested time range: {0} to {1};'
               ' found data from time range: {2} to {3}.')
    # Add tolerance of one second, due to precision of cftime.datetimes
    tol = datetime.timedelta(seconds=1)
    if isinstance(da_start, np.datetime64):
        tol = np.timedelta64(tol, 'ns')
    range_exists = ((da_start - tol) <= start_date and
                    (da_end + tol) >= end_date)
    assert (range_exists), message.format(start_date, end_date,
                                          da_start, da_end)


def sel_time(da, start_date, end_date):
    """Subset a DataArray or Dataset for a given date range.

    Ensures that data are present for full extent of requested range.
    Appends start and end date of the subset to the DataArray.

    Parameters
    ----------
    da : DataArray or Dataset
        data to subset
    start_date : np.datetime64
        start of date interval
    end_date : np.datetime64
        end of date interval

    Returns
    ----------
    da : DataArray or Dataset
        subsetted data

    Raises
    ------
    AssertionError
        if data for requested range do not exist for part or all of
        requested range
    """
    _assert_has_data_for_time(da, start_date, end_date)
    da[SUBSET_START_DATE_STR] = xr.DataArray(start_date)
    da[SUBSET_END_DATE_STR] = xr.DataArray(end_date)
    return da.sel(**{TIME_STR: slice(start_date, end_date)})


def assert_matching_time_coord(arr1, arr2):
    """Check to see if two DataArrays have the same time coordinate.

    Parameters
    ----------
    arr1 : DataArray or Dataset
        First DataArray or Dataset
    arr2 : DataArray or Dataset
        Second DataArray or Dataset

    Raises
    ------
    ValueError
        If the time coordinates are not identical between the two Datasets
    """
    message = ('Time weights not indexed by the same time coordinate as'
               ' computed data.  This will lead to an improperly computed'
               ' time weighted average.  Exiting.\n'
               'arr1: {}\narr2: {}')
    if not (arr1[TIME_STR].identical(arr2[TIME_STR])):
        raise ValueError(message.format(arr1[TIME_STR], arr2[TIME_STR]))


def ensure_time_as_index(ds):
    """Ensures that time is an indexed coordinate on relevant quantites.

    Sometimes when the data we load from disk has only one timestep, the
    indexing of time-defined quantities in the resulting xarray.Dataset gets
    messed up, in that the time bounds array and data variables don't get
    indexed by time, even though they should.  Therefore, we need this helper
    function to (possibly) correct this.

    Note that this must be applied before CF-conventions are decoded; otherwise
    it casts ``np.datetime64[ns]`` as ``int`` values.

    Parameters
    ----------
    ds : Dataset
        Dataset with a time coordinate

    Returns
    -------
    Dataset

    """
    time_indexed_coords = {TIME_WEIGHTS_STR, TIME_BOUNDS_STR}
    time_indexed_vars = set(ds.data_vars).union(time_indexed_coords)
    time_indexed_vars = time_indexed_vars.intersection(ds.variables)
    for name in time_indexed_vars:
        if TIME_STR not in ds[name].indexes:
            da = ds[name].expand_dims(TIME_STR)
            da[TIME_STR] = ds[TIME_STR]
            ds[name] = da
    return ds


def infer_year(date):
    """Given a datetime-like object or string infer the year.

    Parameters
    ----------
    date : datetime-like object or str
        Input date

    Returns
    -------
    int

    Examples
    --------
    >>> infer_year('2000')
    2000
    >>> infer_year('2000-01')
    2000
    >>> infer_year('2000-01-31')
    2000
    >>> infer_year(datetime.datetime(2000, 1, 1))
    2000
    >>> infer_year(np.datetime64('2000-01-01'))
    2000
    >>> infer_year(DatetimeNoLeap(2000, 1, 1))
    2000
    >>>
    """
    if isinstance(date, str):
        # Look for a string that begins with four numbers; the first four
        # numbers found are the year.
        pattern = r'(?P<year>\d{4})'
        result = re.match(pattern, date)
        if result:
            return int(result.groupdict()['year'])
        else:
            raise ValueError('Invalid date string provided: {}'.format(date))
    elif isinstance(date, np.datetime64):
        return date.item().year
    else:
        return date.year


def maybe_convert_to_index_date_type(index, date):
    """Convert a datetime-like object to the index's date type.

    Datetime indexing in xarray can be done using either a pandas
    DatetimeIndex or a CFTimeIndex.  Both support partial-datetime string
    indexing regardless of the calendar type of the underlying data;
    therefore if a string is passed as a date, we return it unchanged.  If a
    datetime-like object is provided, it will be converted to the underlying
    date type of the index.  For a DatetimeIndex that is np.datetime64; for a
    CFTimeIndex that is an object of type cftime.datetime specific to the
    calendar used.

    Parameters
    ----------
    index : pd.Index
        Input time index
    date : datetime-like object or str
        Input datetime

    Returns
    -------
    date of the type appropriate for the time index of the Dataset
    """
    if isinstance(date, str):
        return date

    if isinstance(index, pd.DatetimeIndex):
        if isinstance(date, np.datetime64):
            return date
        else:
            return np.datetime64(str(date))
    else:
        date_type = index.date_type
        if isinstance(date, date_type):
            return date
        else:
            if isinstance(date, np.datetime64):
                # Convert to datetime.date or datetime.datetime object
                date = date.item()

            if isinstance(date, datetime.date):
                # Convert to a datetime.datetime object
                date = datetime.datetime.combine(
                    date, datetime.datetime.min.time())

            return date_type(date.year, date.month, date.day, date.hour,
                             date.minute, date.second, date.microsecond)
