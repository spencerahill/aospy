#!/usr/bin/env python
"""Test suite for aospy.timedate module."""
import datetime
import warnings

import cftime
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from itertools import product

from aospy.data_loader import set_grid_attrs_as_coords
from aospy.internal_names import (
    TIME_STR, TIME_BOUNDS_STR, BOUNDS_STR, TIME_WEIGHTS_STR,
    RAW_START_DATE_STR, RAW_END_DATE_STR, SUBSET_START_DATE_STR,
    SUBSET_END_DATE_STR
)
from aospy.automate import _merge_dicts
from aospy.utils.times import (
    apply_time_offset,
    average_time_bounds,
    monthly_mean_ts,
    monthly_mean_at_each_ind,
    ensure_datetime,
    datetime_or_default,
    month_indices,
    _month_conditional,
    extract_months,
    ensure_time_avg_has_cf_metadata,
    _assert_has_data_for_time,
    add_uniform_time_weights,
    assert_matching_time_coord,
    ensure_time_as_index,
    sel_time,
    yearly_average,
    infer_year,
    maybe_convert_to_index_date_type
)


_INVALID_DATE_OBJECTS = [1985, True, None]


def test_apply_time_offset():
    start = datetime.datetime(1900, 5, 10)
    years, months, days, hours = -2, 1, 7, 3
    # test lengths 0, 1, and >1 of input time array
    for periods in range(3):
        times = pd.date_range(start=start, freq='M', periods=periods)
        actual = apply_time_offset(xr.DataArray(times), years=years,
                                   months=months, days=days, hours=hours)
        desired = (times + pd.tseries.offsets.DateOffset(
            years=years, months=months, days=days, hours=hours
        ))
        assert actual.identical(desired)


def test_monthly_mean_ts_single_month():
    time = pd.date_range('2000-01-01', freq='6H', periods=4 * 31)
    arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                       coords={TIME_STR: time})
    desired = arr.mean(TIME_STR)
    actual = monthly_mean_ts(arr)
    np.testing.assert_allclose(actual, desired)


def test_monthly_mean_ts_submonthly():
    time = pd.date_range('2000-01-01', freq='1D', periods=365 * 3)
    arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                       coords={TIME_STR: time})
    desired = arr.resample(**{TIME_STR: '1M'}).mean(TIME_STR)
    actual = monthly_mean_ts(arr)
    assert desired.identical(actual)


def test_monthly_mean_ts_monthly():
    time = pd.date_range('2000-01-01', freq='1M', periods=120)
    arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                       coords={TIME_STR: time})
    actual = monthly_mean_ts(arr)
    assert arr.identical(actual)


@pytest.mark.filterwarnings('ignore:Mean of empty slice')
def test_monthly_mean_ts_na():
    time = pd.to_datetime(['2000-06-01', '2001-06-01'])
    arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                       coords={TIME_STR: time})
    arr = arr.resample(**{TIME_STR: '1M'}).mean(TIME_STR)
    actual = monthly_mean_ts(arr)
    desired = arr.dropna(TIME_STR)
    assert desired.identical(actual)


def test_monthly_mean_at_each_ind():
    times_submonthly = pd.to_datetime(['2000-06-01', '2000-06-15',
                                       '2000-07-04', '2000-07-19'])
    times_means = pd.to_datetime(['2000-06-01', '2000-07-01'])
    len_other_dim = 2
    arr_submonthly = xr.DataArray(
        np.random.random((len(times_submonthly), len_other_dim)),
        dims=[TIME_STR, 'dim0'], coords={TIME_STR: times_submonthly}
    )
    arr_means = xr.DataArray(
        np.random.random((len(times_means), len_other_dim)),
        dims=arr_submonthly.dims, coords={TIME_STR: times_means}
    )
    actual = monthly_mean_at_each_ind(arr_means, arr_submonthly)
    desired_values = np.stack([arr_means.values[0]] * len_other_dim +
                              [arr_means.values[1]] * len_other_dim,
                              axis=0)
    desired = xr.DataArray(desired_values, dims=arr_submonthly.dims,
                           coords=arr_submonthly.coords)
    assert actual.identical(desired)


@pytest.mark.parametrize('date', [np.datetime64('2000-01-01'),
                                  cftime.DatetimeNoLeap(1, 1, 1),
                                  datetime.datetime(1, 1, 1),
                                  '2000-01-01'])
def test_ensure_datetime_valid_input(date):
    assert ensure_datetime(date) == date


def test_ensure_datetime_invalid_input():
    with pytest.raises(TypeError):
        for obj in _INVALID_DATE_OBJECTS:
            ensure_datetime(obj)


def test_datetime_or_default():
    date = np.datetime64('2000-01-01')
    assert datetime_or_default(None, 'dummy') == 'dummy'
    assert datetime_or_default(date, 'dummy') == ensure_datetime(date)


def test_month_indices():
    np.testing.assert_array_equal(month_indices('ann'), range(1, 13))
    np.testing.assert_array_equal(month_indices('jja'),
                                  np.array([6, 7, 8]))
    with pytest.raises(ValueError):
        month_indices('dfm')
        month_indices('j')
        month_indices('q')
    assert month_indices('s') == [9]
    np.testing.assert_array_equal(month_indices('djf'),
                                  np.array([12, 1, 2]))
    assert month_indices(12) == [12]


def test_month_conditional():
    test = pd.date_range('2000-01-01', '2000-03-01', freq='M')
    test = xr.DataArray(test, dims=[TIME_STR], coords=[test])
    result_jan = _month_conditional(test, [1])
    np.testing.assert_array_equal(result_jan, np.array([True, False]))

    result_jan_feb = _month_conditional(test, [1, 2])
    np.testing.assert_array_equal(result_jan_feb, np.array([True, True]))

    result_march = _month_conditional(test, [3])
    np.testing.assert_array_equal(result_march, np.array([False, False]))

    test = pd.date_range('1999-12-31 18:00:00', '2000-01-01 00:00:00',
                         freq='6H')
    test = xr.DataArray(test, dims=[TIME_STR])
    result_jan = _month_conditional(test, [1])
    np.testing.assert_array_equal(result_jan,
                                  np.array([False, True]))

    result_jd = _month_conditional(test, [1, 12])
    np.testing.assert_array_equal(result_jd,
                                  np.array([True, True]))

    # Test month not in range
    result_march = _month_conditional(test, [3])
    np.testing.assert_array_equal(result_march,
                                  np.array([False, False]))


def test_extract_months():
    time = xr.DataArray(pd.date_range(start='2001-02-18', end='2002-07-12',
                                      freq='1D'), dims=[TIME_STR])
    months = 'mam'  # March-April-May
    desired = xr.concat([
        xr.DataArray(pd.date_range(start='2001-03-01', end='2001-05-31',
                                   freq='1D'), dims=[TIME_STR]),
        xr.DataArray(pd.date_range(start='2002-03-01', end='2002-05-31',
                                   freq='1D'), dims=[TIME_STR])
    ], dim=TIME_STR)
    actual = extract_months(time, months)
    xr.testing.assert_identical(actual, desired)


def test_extract_months_single_month():
    time = xr.DataArray(pd.date_range(start='1678-01-01', end='1678-01-31',
                                      freq='1M'), dims=[TIME_STR])
    months = 1
    desired = time
    actual = extract_months(time, months)
    xr.testing.assert_identical(actual, desired)


@pytest.fixture
def ds_time_encoded_cf():
    time_bounds = np.array([[0, 31], [31, 59], [59, 90]])
    nv = np.array([0, 1])
    time = np.array([15, 46, 74])
    data = np.zeros((3))
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name='a').to_dataset()
    ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                       coords=[time, nv],
                                       dims=[TIME_STR, BOUNDS_STR],
                                       name=TIME_BOUNDS_STR)
    units_str = 'days since 2000-01-01 00:00:00'
    cal_str = 'noleap'
    ds[TIME_STR].attrs['units'] = units_str
    ds[TIME_STR].attrs['calendar'] = cal_str
    return ds


def test_ensure_time_avg_has_cf_metadata(ds_time_encoded_cf):
    ds = ds_time_encoded_cf
    time = ds[TIME_STR].values
    time_bounds = ds[TIME_BOUNDS_STR].values
    units_str = ds[TIME_STR].attrs['units']
    cal_str = ds[TIME_STR].attrs['calendar']

    with pytest.raises(KeyError):
        ds[TIME_BOUNDS_STR].attrs['units']
    with pytest.raises(KeyError):
        ds[TIME_BOUNDS_STR].attrs['calendar']

    ds = ensure_time_avg_has_cf_metadata(ds)

    result = ds[TIME_BOUNDS_STR].attrs['units']
    assert result == units_str
    result = ds[TIME_BOUNDS_STR].attrs['calendar']
    assert result == cal_str

    avg_DT_data = np.diff(time_bounds, axis=1).squeeze()
    average_DT_expected = xr.DataArray(avg_DT_data,
                                       coords=[time],
                                       dims=[TIME_STR],
                                       name=TIME_WEIGHTS_STR)
    average_DT_expected[TIME_STR].attrs['units'] = units_str
    average_DT_expected.attrs['units'] = 'days'
    average_DT_expected[TIME_STR].attrs['calendar'] = cal_str
    assert ds[TIME_WEIGHTS_STR].identical(average_DT_expected)

    assert ds[RAW_START_DATE_STR].values == [0]
    assert ds[RAW_START_DATE_STR].attrs['units'] == units_str
    assert ds[RAW_START_DATE_STR].attrs['calendar'] == cal_str

    assert ds[RAW_END_DATE_STR].values == [90]
    assert ds[RAW_END_DATE_STR].attrs['units'] == units_str
    assert ds[RAW_END_DATE_STR].attrs['calendar'] == cal_str


def test_add_uniform_time_weights():
    time = np.array([15, 46, 74])
    data = np.zeros((3))
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name='a').to_dataset()
    units_str = 'days since 2000-01-01 00:00:00'
    cal_str = 'noleap'
    ds[TIME_STR].attrs['units'] = units_str
    ds[TIME_STR].attrs['calendar'] = cal_str

    with pytest.raises(KeyError):
        ds[TIME_WEIGHTS_STR]

    ds = add_uniform_time_weights(ds)
    time_weights_expected = xr.DataArray(
        [1, 1, 1], coords=ds[TIME_STR].coords, name=TIME_WEIGHTS_STR)
    time_weights_expected.attrs['units'] = 'days'
    assert ds[TIME_WEIGHTS_STR].identical(time_weights_expected)


def test_assert_has_data_for_time():
    time_bounds = np.array([[0, 31], [31, 59], [59, 90]])
    nv = np.array([0, 1])
    time = np.array([15, 46, 74])
    data = np.zeros((3))
    var_name = 'a'
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name=var_name).to_dataset()
    ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                       coords=[time, nv],
                                       dims=[TIME_STR, BOUNDS_STR],
                                       name=TIME_BOUNDS_STR)
    units_str = 'days since 2000-01-01 00:00:00'
    ds[TIME_STR].attrs['units'] = units_str
    ds = ensure_time_avg_has_cf_metadata(ds)
    ds = set_grid_attrs_as_coords(ds)
    ds = xr.decode_cf(ds)
    da = ds[var_name]

    start_date = np.datetime64('2000-01-01')
    end_date = np.datetime64('2000-03-31')
    _assert_has_data_for_time(da, start_date, end_date)

    start_date_bad = np.datetime64('1999-12-31')
    end_date_bad = np.datetime64('2000-04-01')

    with pytest.raises(AssertionError):
        _assert_has_data_for_time(da, start_date_bad, end_date)

    with pytest.raises(AssertionError):
        _assert_has_data_for_time(da, start_date, end_date_bad)

    with pytest.raises(AssertionError):
        _assert_has_data_for_time(da, start_date_bad, end_date_bad)


_CFTIME_DATE_TYPES = {
    'noleap': cftime.DatetimeNoLeap,
    '365_day': cftime.DatetimeNoLeap,
    '360_day': cftime.Datetime360Day,
    'julian': cftime.DatetimeJulian,
    'all_leap': cftime.DatetimeAllLeap,
    '366_day': cftime.DatetimeAllLeap,
    'gregorian': cftime.DatetimeGregorian,
    'proleptic_gregorian': cftime.DatetimeProlepticGregorian
}


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
@pytest.mark.filterwarnings('ignore:Unable to decode')
@pytest.mark.parametrize(['calendar', 'date_type'],
                         list(_CFTIME_DATE_TYPES.items()))
def test_assert_has_data_for_time_cftime_datetimes(calendar, date_type):
    time_bounds = np.array([[0, 2], [2, 4], [4, 6]])
    nv = np.array([0, 1])
    time = np.array([1, 3, 5])
    data = np.zeros((3))
    var_name = 'a'
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name=var_name).to_dataset()
    ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                       coords=[time, nv],
                                       dims=[TIME_STR, BOUNDS_STR],
                                       name=TIME_BOUNDS_STR)
    units_str = 'days since 0002-01-02 00:00:00'
    ds[TIME_STR].attrs['units'] = units_str
    ds[TIME_STR].attrs['calendar'] = calendar
    ds = ensure_time_avg_has_cf_metadata(ds)
    ds = set_grid_attrs_as_coords(ds)

    with warnings.catch_warnings(record=True):
        with xr.set_options(enable_cftimeindex=True):
            ds = xr.decode_cf(ds)
    da = ds[var_name]

    start_date = date_type(2, 1, 2)
    end_date = date_type(2, 1, 8)

    _assert_has_data_for_time(da, start_date, end_date)

    start_date_bad = date_type(2, 1, 1)
    end_date_bad = date_type(2, 1, 9)

    with pytest.raises(AssertionError):
        _assert_has_data_for_time(da, start_date_bad, end_date)

    with pytest.raises(AssertionError):
        _assert_has_data_for_time(da, start_date, end_date_bad)

    with pytest.raises(AssertionError):
        _assert_has_data_for_time(da, start_date_bad, end_date_bad)


def test_assert_has_data_for_time_str_input():
    time_bounds = np.array([[0, 31], [31, 59], [59, 90]])
    nv = np.array([0, 1])
    time = np.array([15, 46, 74])
    data = np.zeros((3))
    var_name = 'a'
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name=var_name).to_dataset()
    ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                       coords=[time, nv],
                                       dims=[TIME_STR, BOUNDS_STR],
                                       name=TIME_BOUNDS_STR)
    units_str = 'days since 2000-01-01 00:00:00'
    ds[TIME_STR].attrs['units'] = units_str
    ds = ensure_time_avg_has_cf_metadata(ds)
    ds = set_grid_attrs_as_coords(ds)
    ds = xr.decode_cf(ds)
    da = ds[var_name]

    start_date = '2000-01-01'
    end_date = '2000-03-31'
    _assert_has_data_for_time(da, start_date, end_date)

    start_date_bad = '1999-12-31'
    end_date_bad = '2000-04-01'

    # With strings these checks are disabled
    _assert_has_data_for_time(da, start_date_bad, end_date)
    _assert_has_data_for_time(da, start_date, end_date_bad)
    _assert_has_data_for_time(da, start_date_bad, end_date_bad)


def test_assert_matching_time_coord():
    rng = pd.date_range('2000-01-01', '2001-01-01', freq='M')
    arr1 = xr.DataArray(rng, coords=[rng], dims=[TIME_STR])
    arr2 = xr.DataArray(rng, coords=[rng], dims=[TIME_STR])
    assert_matching_time_coord(arr1, arr2)

    arr2 = arr2.sel(**{TIME_STR: slice('2000-03', '2000-05')})
    with pytest.raises(ValueError):
        assert_matching_time_coord(arr1, arr2)


def test_ensure_time_as_index_no_change():
    # Already properly indexed, so shouldn't be modified.
    arr = xr.DataArray([-23, 42.4], coords=[[1, 2]], dims=[TIME_STR])
    arr[TIME_STR].attrs['units'] = 'days since 2000-01-01 00:00:00'
    arr[TIME_STR].attrs['calendar'] = 'standard'
    ds = arr.to_dataset(name='a')
    ds.coords[TIME_WEIGHTS_STR] = xr.DataArray(
        [1, 1], dims=[TIME_STR], coords={TIME_STR: arr[TIME_STR]}
    )
    ds.coords[TIME_BOUNDS_STR] = xr.DataArray(
        [[0.5, 1.5], [1.5, 2.5]], dims=[TIME_STR, BOUNDS_STR],
        coords={TIME_STR: arr[TIME_STR]}
    )
    xr.testing.assert_identical(ds, ensure_time_as_index(ds))


def test_ensure_time_as_index_with_change():
    # Time bounds array doesn't index time initially, which gets fixed.
    arr = xr.DataArray([-93], dims=[TIME_STR], coords={TIME_STR: [3]})
    arr[TIME_STR].attrs['units'] = 'days since 2000-01-01 00:00:00'
    arr[TIME_STR].attrs['calendar'] = 'standard'
    ds = arr.to_dataset(name='a')
    ds.coords[TIME_WEIGHTS_STR] = xr.DataArray(
        [1], dims=[TIME_STR], coords={TIME_STR: arr[TIME_STR]}
    )
    ds.coords[TIME_BOUNDS_STR] = xr.DataArray(
        [[3.5, 4.5]], dims=[TIME_STR, BOUNDS_STR],
        coords={TIME_STR: arr[TIME_STR]}
    )
    ds = ds.isel(**{TIME_STR: 0}).expand_dims(TIME_STR)
    actual = ensure_time_as_index(ds)
    expected = arr.to_dataset(name='a')
    expected.coords[TIME_WEIGHTS_STR] = xr.DataArray(
        [1], dims=[TIME_STR], coords={TIME_STR: arr[TIME_STR]}
    )
    expected.coords[TIME_BOUNDS_STR] = xr.DataArray(
        [[3.5, 4.5]], dims=[TIME_STR, BOUNDS_STR],
        coords={TIME_STR: arr[TIME_STR]}
        )
    xr.testing.assert_identical(actual, expected)


def test_sel_time():
    time_bounds = np.array([[0, 31], [31, 59], [59, 90]])
    nv = np.array([0, 1])
    time = np.array([15, 46, 74])
    data = np.zeros((3))
    var_name = 'a'
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name=var_name).to_dataset()
    ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                       coords=[time, nv],
                                       dims=[TIME_STR, BOUNDS_STR],
                                       name=TIME_BOUNDS_STR)
    units_str = 'days since 2000-01-01 00:00:00'
    ds[TIME_STR].attrs['units'] = units_str
    ds = ensure_time_avg_has_cf_metadata(ds)
    ds = set_grid_attrs_as_coords(ds)
    ds = xr.decode_cf(ds)
    da = ds[var_name]

    start_date = np.datetime64('2000-02-01')
    end_date = np.datetime64('2000-03-31')
    result = sel_time(da, start_date, end_date)
    assert result[SUBSET_START_DATE_STR].values == start_date
    assert result[SUBSET_END_DATE_STR].values == end_date


def test_yearly_average_no_mask():
    times = pd.to_datetime(['2000-06-01', '2000-06-15',
                            '2001-07-04', '2001-10-01', '2001-12-31',
                            '2004-01-01'])
    arr = xr.DataArray(np.random.random((len(times),)),
                       dims=[TIME_STR], coords={TIME_STR: times})
    dt = arr.copy(deep=True)
    dt.values = np.random.random((len(times),))

    actual = yearly_average(arr, dt)

    yr2000 = (arr[0]*dt[0] + arr[1]*dt[1]) / (dt[0] + dt[1])
    yr2001 = ((arr[2]*dt[2] + arr[3]*dt[3] + arr[4]*dt[4]) /
              (dt[2] + dt[3] + dt[4]))
    yr2004 = arr[-1]
    yrs_coord = [2000, 2001, 2004]
    yr_avgs = np.array([yr2000, yr2001, yr2004])
    desired = xr.DataArray(yr_avgs, dims=['year'], coords={'year': yrs_coord})
    xr.testing.assert_allclose(actual, desired)


def test_yearly_average_masked_data():
    times = pd.to_datetime(['2000-06-01', '2000-06-15',
                            '2001-07-04', '2001-10-01', '2001-12-31',
                            '2004-01-01'])
    arr = xr.DataArray(np.random.random((len(times),)),
                       dims=[TIME_STR], coords={TIME_STR: times})
    arr[0] = -999
    arr = arr.where(arr != -999)
    dt = arr.copy(deep=True)
    dt.values = np.random.random((len(times),))

    actual = yearly_average(arr, dt)

    yr2000 = arr[1]
    yr2001 = ((arr[2]*dt[2] + arr[3]*dt[3] + arr[4]*dt[4]) /
              (dt[2] + dt[3] + dt[4]))
    yr2004 = arr[-1]
    yrs_coord = [2000, 2001, 2004]
    yr_avgs = np.array([yr2000, yr2001, yr2004])
    desired = xr.DataArray(yr_avgs, dims=['year'], coords={'year': yrs_coord})
    xr.testing.assert_allclose(actual, desired)


def test_average_time_bounds(ds_time_encoded_cf):
    ds = ds_time_encoded_cf
    actual = average_time_bounds(ds)[TIME_STR]

    desired_values = ds[TIME_BOUNDS_STR].mean(dim=BOUNDS_STR).values
    desired = xr.DataArray(desired_values, dims=[TIME_STR],
                           coords={TIME_STR: desired_values}, name=TIME_STR)

    xr.testing.assert_identical(actual, desired)


_INFER_YEAR_TESTS = [
    (np.datetime64('2000-01-01'), 2000),
    (datetime.datetime(2000, 1, 1), 2000),
    ('2000', 2000),
    ('2000-01', 2000),
    ('2000-01-01', 2000)
]
_INFER_YEAR_TESTS = _INFER_YEAR_TESTS + [
    (date_type(2000, 1, 1), 2000) for date_type in _CFTIME_DATE_TYPES.values()]


@pytest.mark.parametrize(
    ['date', 'expected'],
    _INFER_YEAR_TESTS)
def test_infer_year(date, expected):
    assert infer_year(date) == expected


@pytest.mark.parametrize('date', ['-0001', 'A001', '01'])
def test_infer_year_invalid(date):
    with pytest.raises(ValueError):
        infer_year(date)


_DATETIME_INDEX = pd.date_range('2000-01-01', freq='M', periods=1)
_DATETIME_CONVERT_TESTS = {}
for date_label, date_type in _CFTIME_DATE_TYPES.items():
    key = 'DatetimeIndex-{}'.format(date_label)
    _DATETIME_CONVERT_TESTS[key] = (_DATETIME_INDEX, date_type(2000, 1, 1),
                                    np.datetime64('2000-01'))
_NON_CFTIME_DATES = {
    'datetime.datetime': datetime.datetime(2000, 1, 1),
    'np.datetime64': np.datetime64('2000-01-01'),
    'str': '2000'
}
for date_label, date in _NON_CFTIME_DATES.items():
    key = 'DatetimeIndex-{}'.format(date_label)
    if isinstance(date, str):
        _DATETIME_CONVERT_TESTS[key] = (_DATETIME_INDEX, date, date)
    else:
        _DATETIME_CONVERT_TESTS[key] = (_DATETIME_INDEX, date,
                                        np.datetime64('2000-01'))

_CFTIME_INDEXES = {
    'CFTimeIndex[{}]'.format(key): xr.CFTimeIndex([value(1, 1, 1)]) for
    key, value in _CFTIME_DATE_TYPES.items()
}
_CFTIME_CONVERT_TESTS = {}
for ((index_label, index),
     (date_label, date_type)) in product(_CFTIME_INDEXES.items(),
                                         _CFTIME_DATE_TYPES.items()):
    key = '{}-{}'.format(index_label, date_label)
    _CFTIME_CONVERT_TESTS[key] = (index, date_type(1, 1, 1),
                                  index.date_type(1, 1, 1))
_NON_CFTIME_DATES_0001 = {
    'datetime.datetime': datetime.datetime(1, 1, 1),
    'np.datetime64': np.datetime64('0001-01-01'),
    'str': '0001'
}
for ((idx_label, index),
     (date_label, date)) in product(_CFTIME_INDEXES.items(),
                                    _NON_CFTIME_DATES_0001.items()):
    key = '{}-{}'.format(index_label, date_label)
    if isinstance(date, str):
        _CFTIME_CONVERT_TESTS[key] = (index, date, date)
    else:
        _CFTIME_CONVERT_TESTS[key] = (index, date, index.date_type(1, 1, 1))

_CONVERT_DATE_TYPE_TESTS = _merge_dicts(_DATETIME_CONVERT_TESTS,
                                        _CFTIME_CONVERT_TESTS)


@pytest.mark.parametrize(['index', 'date', 'expected'],
                         list(_CONVERT_DATE_TYPE_TESTS.values()),
                         ids=list(_CONVERT_DATE_TYPE_TESTS.keys()))
def test_maybe_convert_to_index_date_type(index, date, expected):
    result = maybe_convert_to_index_date_type(index, date)
    assert result == expected
