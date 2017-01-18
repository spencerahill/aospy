"""Utility functions for handling times, dates, etc."""
import datetime

import numpy as np
import pandas as pd
import xarray as xr

from .. import internal_names


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
    return arr.resample(
        '1M', internal_names.TIME_STR,
        how='mean').dropna(internal_names.TIME_STR)


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
    time = monthly_means[internal_names.TIME_STR]
    start = time.indexes[internal_names.TIME_STR][0].replace(day=1, hour=0)
    end = time.indexes[internal_names.TIME_STR][-1]
    new_indices = pd.DatetimeIndex(start=start, end=end, freq='MS')
    arr_new = monthly_means.reindex(time=new_indices, method='backfill')
    return arr_new.reindex_like(sub_monthly_timeseries, method='pad')


def ensure_datetime(obj):
    """Return the object if it is of type datetime.datetime; else raise.

    Parameters
    ----------
    obj : Object to be tested.

    Returns
    -------
    The original object if it is a datetime.datetime object.

    Raises
    ------
    TypeError if `obj` is not of type `datetime.datetime`.
    """
    if isinstance(obj, datetime.datetime):
        return obj
    raise TypeError("`datetime.datetime` object required.  "
                    "Type given: {}".format(type(obj)))


def datetime_or_default(date, default):
    """Return a datetime.datetime object or a default.

    Parameters
    ----------
    date : `None` or datetime-like object
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


def numpy_datetime_range_workaround(date):
    """Reset a date to earliest allowable year if outside of valid range.

    Hack to address np.datetime64, and therefore pandas and xarray, not
    supporting dates outside the range 1677-09-21 and 2262-04-11 due to
    nanosecond precision.  See e.g.
    https://github.com/spencerahill/aospy/issues/96.


    Parameters
    ----------
    date : datetime.datetime object

    Returns
    -------
    datetime.datetime object
        Original datetime.datetime object if the original date is within the
        permissible dates, otherwise a datetime.datetime object with the year
        offset to the earliest allowable year.

    """
    if date < pd.Timestamp.min:
        return datetime.datetime(pd.Timestamp.min.year + 1, date.month,
                                 date.day)
    return date


def numpy_datetime_workaround_encode_cf(ds):
    """Generate CF-compliant units for out-of-range dates.

    Hack to address np.datetime64, and therefore pandas and xarray, not
    supporting dates outside the range 1677-09-21 and 2262-04-11 due to
    nanosecond precision.  See e.g.
    https://github.com/spencerahill/aospy/issues/96.

    Specifically, we coerce the data such that, when decoded, the earliest
    value starts in 1678 but with its month, day, and shorter timescales
    (hours, minutes, seconds, etc.) intact and with the time-spacing between
    values intact.

    Parameters
    ----------
    ds : xarray.Dataset

    Returns
    -------
    xarray.Dataset

    """
    time = ds[internal_names.TIME_STR]
    units = time.attrs['units']
    units_yr = units.split(' since ')[1].split('-')[0]
    min_yr_decoded = xr.decode_cf(time.to_dataset('dummy'))
    min_yr = min_yr_decoded[internal_names.TIME_STR].values[0].year
    new_units_yr = pd.Timestamp.min.year + 2 - min_yr
    new_units = units.replace(units_yr, str(new_units_yr))

    for VAR_STR in internal_names.TIME_VAR_STRS:
        if VAR_STR in ds:
            var = ds[VAR_STR]
            var.attrs['units'] = new_units
    return ds


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
        cond |= (time['{}.month'.format(internal_names.TIME_STR)] == month)
    return cond


def create_monthly_time_array(start_date, end_date, months):
    """Create an array of months compliant with numpy datetime limited range.

    Parameters
    ----------
    start_date, end_date : datetime.datetime
    months : int, str, or xarray.DataArray of times
        Specifies which months in each year to include via `_month_conditonal`

    Returns
    -------
    xarray.DataArray
        Array spans the specified date range and includes only the specified
        months of the year.  Year is also reset to start at the beginning of
        the range of dates that np.datetime64 can support if originally out of
        range.

    See Also
    --------
    _month_conditional : Generate boolean array specifying desired months

    """
    start_compliant = numpy_datetime_range_workaround(start_date)
    end_compliant = start_compliant + (end_date - start_date)
    all_months = pd.date_range(start=start_compliant, end=end_compliant,
                               freq='M')
    time = xr.DataArray(all_months, dims=[internal_names.TIME_STR])
    return time[_month_conditional(time, months)]


def extract_date_range_and_months(time, start_date, end_date, months):
    """Extract times within a specified date range and months of the year.

    Parameters
    ----------
    time : xarray.DataArray
         Array of times that can be represented by numpy.datetime64 objects
         (i.e. the year is between 1678 and 2262).
    start_date, end_date : datetime.datetime
        Desired start and end date
    end_date : Desired end date
    months : Desired months of the year to include

    Returns
    -------
    xarray.DataArray of the desired times
    """
    inds = _month_conditional(time, months)
    inds &= (time[internal_names.TIME_STR] <= np.datetime64(end_date))
    inds &= (time[internal_names.TIME_STR] >= np.datetime64(start_date))
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
    """
    AVG_START_DATE_STR = internal_names.AVG_START_DATE_STR
    AVG_END_DATE_STR = internal_names.AVG_END_DATE_STR
    TIME_BOUNDS_STR = internal_names.TIME_BOUNDS_STR
    TIME_STR = internal_names.TIME_STR
    NV_STR = internal_names.NV_STR
    AVERAGE_DT_STR = internal_names.AVERAGE_DT_STR

    if AVERAGE_DT_STR not in ds:
        average_DT = ds[TIME_BOUNDS_STR].diff(NV_STR)
        average_DT = average_DT.rename(AVERAGE_DT_STR).squeeze()
        ds[AVERAGE_DT_STR] = average_DT.drop(NV_STR)

    avg_start_date = ds[TIME_BOUNDS_STR].isel(**{TIME_STR: 0, NV_STR: 0})
    ds[AVG_START_DATE_STR] = avg_start_date.drop([TIME_STR, NV_STR])
    avg_end_date = ds[TIME_BOUNDS_STR].isel(**{TIME_STR: -1, NV_STR: 1})
    ds[AVG_END_DATE_STR] = avg_end_date.drop([TIME_STR, NV_STR])

    for coord in [TIME_BOUNDS_STR, AVG_START_DATE_STR, AVG_END_DATE_STR]:
        ds[coord].attrs['units'] = ds[TIME_STR].attrs['units']
        if 'calendar' in ds[TIME_STR].attrs:
            ds[coord].attrs['calendar'] = ds[TIME_STR].attrs['calendar']
    return ds


def _assert_has_data_for_time(da, start_date, end_date):
    """Check to make sure data is in Dataset for the given time range.

    Parameters
    ----------
    da : DataArray
         DataArray with a time variable
    start_date : netCDF4.netcdftime or np.datetime64
         start date
    end_date : netCDF4.netcdftime or np.datetime64
         end date

    Raises
    ------
    AssertionError
         if the time range is not within the time range of the DataArray
    """
    if internal_names.AVG_START_DATE_STR in da:
        da_start = da[internal_names.AVG_START_DATE_STR].values
        da_end = da[internal_names.AVG_END_DATE_STR].values
    else:
        times = da.time.isel(**{internal_names.TIME_STR: [0, -1]})
        da_start, da_end = times.values
    message = 'Data does not exist for requested time range: {0} to {1}'
    range_exists = start_date >= da_start and end_date <= da_end
    assert (range_exists), message.format(start_date, end_date)


def sel_time(da, start_date, end_date):
    """Subset a DataArray or Dataset for a given date range.

    Ensures that data are present for full extent of requested range.

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
    return da.sel(**{internal_names.TIME_STR: slice(start_date, end_date)})
