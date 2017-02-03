#!/usr/bin/env python
"""Test suite for aospy.timedate module."""
import datetime
import sys
import unittest

import numpy as np
import pandas as pd
import xarray as xr

from aospy.data_loader import set_grid_attrs_as_coords
from aospy.internal_names import (
    TIME_STR, TIME_BOUNDS_STR, BOUNDS_STR, TIME_WEIGHTS_STR,
    RAW_START_DATE_STR, RAW_END_DATE_STR, SUBSET_START_DATE_STR,
    SUBSET_END_DATE_STR
)
from aospy.utils.times import (
    apply_time_offset,
    monthly_mean_ts,
    monthly_mean_at_each_ind,
    ensure_datetime,
    datetime_or_default,
    numpy_datetime_range_workaround,
    numpy_datetime_workaround_encode_cf,
    month_indices,
    _month_conditional,
    extract_months,
    ensure_time_avg_has_cf_metadata,
    _assert_has_data_for_time,
    add_uniform_time_weights,
    assert_matching_time_coord,
    ensure_time_as_dim,
    convert_scalar_to_indexable_coord,
    sel_time
)


_INVALID_DATE_OBJECTS = [1985, True, None, '2016-04-07', np.datetime64(1, 'Y')]


class UtilsTimesTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestUtilsTimes(UtilsTimesTestCase):
    def test_apply_time_offset(self):
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

    def test_monthly_mean_ts_single_month(self):
        time = pd.date_range('2000-01-01', freq='6H', periods=4 * 31)
        arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                           coords={TIME_STR: time})
        desired = arr.mean(TIME_STR)
        actual = monthly_mean_ts(arr)
        np.testing.assert_allclose(actual, desired)

    def test_monthly_mean_ts_submonthly(self):
        time = pd.date_range('2000-01-01', freq='1D', periods=365 * 3)
        arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                           coords={TIME_STR: time})
        desired = arr.resample('1M', TIME_STR, how='mean')
        actual = monthly_mean_ts(arr)
        assert desired.identical(actual)

    def test_monthly_mean_ts_monthly(self):
        time = pd.date_range('2000-01-01', freq='1M', periods=120)
        arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                           coords={TIME_STR: time})
        actual = monthly_mean_ts(arr)
        assert arr.identical(actual)

    def test_monthly_mean_ts_na(self):
        time = pd.to_datetime(['2000-06-01', '2001-06-01'])
        arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                           coords={TIME_STR: time}).resample('1M', TIME_STR)
        actual = monthly_mean_ts(arr)
        desired = arr.dropna(TIME_STR)
        assert desired.identical(actual)

    def test_monthly_mean_at_each_ind(self):
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

    def test_ensure_datetime_valid_input(self):
        for date in [datetime.datetime(1981, 7, 15),
                     datetime.datetime(1, 1, 1)]:
            self.assertEqual(ensure_datetime(date), date)

    def test_ensure_datetime_invalid_input(self):
        with self.assertRaises(TypeError):
            for obj in _INVALID_DATE_OBJECTS:
                ensure_datetime(obj)

    def test_datetime_or_default(self):
        date = datetime.datetime(1, 2, 3)
        assert datetime_or_default(None, 'dummy') == 'dummy'
        assert datetime_or_default(date, 'dummy') == ensure_datetime(date)

    def test_numpy_datetime_range_workaround(self):
        self.assertEqual(numpy_datetime_range_workaround(
            datetime.datetime(pd.Timestamp.min.year + 1, 1, 1),
            pd.Timestamp.min.year + 1
        ), datetime.datetime(pd.Timestamp.min.year + 1, 1, 1))

        self.assertEqual(
            numpy_datetime_range_workaround(datetime.datetime(3, 1, 1), 1),
            datetime.datetime(pd.Timestamp.min.year + 3, 1, 1)
        )

        self.assertEqual(
            numpy_datetime_range_workaround(datetime.datetime(5, 1, 1), 4),
            datetime.datetime(pd.Timestamp.min.year + 2, 1, 1)
        )

    def test_numpy_datetime_workaround_encode_cf(self):
        # 255169 days from 0001-01-01 corresponds to date 700-02-04.
        days = 255169.
        time = xr.DataArray([days], dims=[TIME_STR])
        ds = xr.Dataset(coords={TIME_STR: time})
        ds[TIME_STR].attrs['units'] = 'days since 0001-01-01 00:00:00'
        ds[TIME_STR].attrs['calendar'] = 'noleap'
        actual, min_year = numpy_datetime_workaround_encode_cf(ds)

        time_desired = xr.DataArray([days], dims=[TIME_STR])
        desired = xr.Dataset(coords={TIME_STR: time_desired})
        desired[TIME_STR].attrs['units'] = (
            'days since {0}-01-01 00:00:00'.format(979)
        )
        desired[TIME_STR].attrs['calendar'] = 'noleap'

        assert actual.identical(desired)
        self.assertEqual(xr.decode_cf(actual).time.values[0],
                         np.datetime64('1678-02-04'))
        self.assertEqual(min_year, 700)

        # Test a case where times are in the Timestamp-valid range
        time = xr.DataArray([10], dims=[TIME_STR])
        ds = xr.Dataset(coords={TIME_STR: time})
        ds[TIME_STR].attrs['units'] = 'days since 2000-01-01 00:00:00'
        ds[TIME_STR].attrs['calendar'] = 'noleap'
        actual, min_year = numpy_datetime_workaround_encode_cf(ds)
        self.assertEqual(xr.decode_cf(actual).time.values[0],
                         np.datetime64('2000-01-11'))
        self.assertEqual(min_year, 2000)

    def test_month_indices(self):
        np.testing.assert_array_equal(month_indices('ann'), range(1, 13))
        np.testing.assert_array_equal(month_indices('jja'),
                                      np.array([6, 7, 8]))
        with self.assertRaises(ValueError):
            month_indices('dfm')
            month_indices('j')
            month_indices('q')
        self.assertEqual(month_indices('s'), [9])
        np.testing.assert_array_equal(month_indices('djf'),
                                      np.array([12, 1, 2]))
        self.assertEqual(month_indices(12), [12])

    def test_month_conditional(self):
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

    def test_extract_months(self):
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

    def test_extract_months_single_month(self):
        time = xr.DataArray(pd.date_range(start='1678-01-01', end='1678-01-31',
                                          freq='1M'), dims=[TIME_STR])
        months = 1
        desired = time
        actual = extract_months(time, months)
        xr.testing.assert_identical(actual, desired)

    def test_ensure_time_avg_has_cf_metadata(self):
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

        with self.assertRaises(KeyError):
            ds[TIME_BOUNDS_STR].attrs['units']
        with self.assertRaises(KeyError):
            ds[TIME_BOUNDS_STR].attrs['calendar']

        ds = ensure_time_avg_has_cf_metadata(ds)

        result = ds[TIME_BOUNDS_STR].attrs['units']
        self.assertEqual(result, units_str)
        result = ds[TIME_BOUNDS_STR].attrs['calendar']
        self.assertEqual(result, cal_str)

        avg_DT_data = np.diff(time_bounds, axis=1).squeeze()
        average_DT_expected = xr.DataArray(avg_DT_data,
                                           coords=[time],
                                           dims=[TIME_STR],
                                           name=TIME_WEIGHTS_STR)
        average_DT_expected[TIME_STR].attrs['units'] = units_str
        average_DT_expected.attrs['units'] = 'days'
        average_DT_expected[TIME_STR].attrs['calendar'] = cal_str
        assert ds[TIME_WEIGHTS_STR].identical(average_DT_expected)

        self.assertEqual(ds[RAW_START_DATE_STR].values, [0])
        self.assertEqual(ds[RAW_START_DATE_STR].attrs['units'], units_str)
        self.assertEqual(ds[RAW_START_DATE_STR].attrs['calendar'], cal_str)

        self.assertEqual(ds[RAW_END_DATE_STR].values, [90])
        self.assertEqual(ds[RAW_END_DATE_STR].attrs['units'], units_str)
        self.assertEqual(ds[RAW_END_DATE_STR].attrs['calendar'], cal_str)

    def test_add_uniform_time_weights(self):
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

        with self.assertRaises(KeyError):
            ds[TIME_WEIGHTS_STR]

        ds = add_uniform_time_weights(ds)
        time_weights_expected = xr.DataArray(
            [1, 1, 1], coords=ds[TIME_STR].coords, name=TIME_WEIGHTS_STR)
        time_weights_expected.attrs['units'] = 'days'
        assert ds[TIME_WEIGHTS_STR].identical(time_weights_expected)

    def test_assert_has_data_for_time(self):
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

        with self.assertRaises(AssertionError):
            _assert_has_data_for_time(da, start_date_bad, end_date)

        with self.assertRaises(AssertionError):
            _assert_has_data_for_time(da, start_date, end_date_bad)

        with self.assertRaises(AssertionError):
            _assert_has_data_for_time(da, start_date_bad, end_date_bad)

    def test_assert_matching_time_coord(self):
        rng = pd.date_range('2000-01-01', '2001-01-01', freq='M')
        arr1 = xr.DataArray(rng, coords=[rng], dims=[TIME_STR])
        arr2 = xr.DataArray(rng, coords=[rng], dims=[TIME_STR])
        assert_matching_time_coord(arr1, arr2)

        arr2 = arr2.sel(**{TIME_STR: slice('2000-03', '2000-05')})
        with self.assertRaises(ValueError):
            assert_matching_time_coord(arr1, arr2)

    def test_ensure_time_as_dim(self):
        arr = xr.DataArray([3, 4], coords=[[1, 2]], dims=[TIME_STR])
        arr[TIME_STR].attrs['units'] = 'days since 2000-01-01 00:00:00'
        arr[TIME_STR].attrs['calendar'] = 'standard'
        ds = arr.to_dataset(name='a')
        assert TIME_STR in ds.dims
        assert ds.identical(ensure_time_as_dim(ds))

        scalar_time_in_ds = ds.isel(**{TIME_STR: 0})
        assert TIME_STR not in scalar_time_in_ds.dims
        result = ensure_time_as_dim(scalar_time_in_ds)

        arr = xr.DataArray([3], coords=[[1]], dims=[TIME_STR])
        arr[TIME_STR].attrs['units'] = 'days since 2000-01-01 00:00:00'
        arr[TIME_STR].attrs['calendar'] = 'standard'
        expected = arr.to_dataset(name='a')
        xr.testing.assert_identical(result, expected)

    def test_convert_scalar_to_indexable_coord(self):
        da = xr.DataArray([3, 4], coords=[[1, 2]], dims=['a'], name='b')
        da['a'].attrs['test'] = 'c'
        scalar_coord = da.isel(a=0)['a']
        assert 'a' not in scalar_coord.dims

        indexable_coord = convert_scalar_to_indexable_coord(scalar_coord)
        assert 'a' in indexable_coord.dims

        expected = xr.DataArray([1], coords=[[1]], dims=['a'], name='a')
        expected.attrs['test'] = 'c'
        xr.testing.assert_identical(indexable_coord, expected)

    def test_sel_time(self):
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
        self.assertEqual(result[SUBSET_START_DATE_STR].values, start_date)
        self.assertEqual(result[SUBSET_END_DATE_STR].values, end_date)


if __name__ == '__main__':
    sys.exit(unittest.main())
