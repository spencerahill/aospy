#!/usr/bin/env python
"""Test suite for aospy.timedate module."""
import datetime
import sys
import unittest

import numpy as np
import xarray as xr
import pandas as pd

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
    create_monthly_time_array,
    extract_date_range_and_months,
)
from aospy import TIME_STR


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
        time = pd.date_range('2000-01-01', freq='6H', periods=4*31)
        arr = xr.DataArray(np.random.random(time.shape), dims=[TIME_STR],
                           coords={TIME_STR: time})
        desired = arr.mean(TIME_STR)
        actual = monthly_mean_ts(arr)
        np.testing.assert_allclose(actual, desired)

    def test_monthly_mean_ts_submonthly(self):
        time = pd.date_range('2000-01-01', freq='1D', periods=365*3)
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
        desired_values = np.stack([arr_means.values[0]]*len_other_dim +
                                  [arr_means.values[1]]*len_other_dim,
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
            datetime.datetime(pd.Timestamp.min.year + 1, 1, 1)
        ), datetime.datetime(pd.Timestamp.min.year + 1, 1, 1))
        self.assertEqual(
            numpy_datetime_range_workaround(datetime.datetime(1, 1, 1)),
            datetime.datetime(pd.Timestamp.min.year + 1, 1, 1)
        )

    def test_numpy_datetime_workaround_encode_cf(self):
        # 255169 days from 0001-01-01 corresponds to date 700-02-04.
        days = 255169.
        time = xr.DataArray([days], dims=[TIME_STR])
        ds = xr.Dataset(coords={TIME_STR: time})
        ds[TIME_STR].attrs['units'] = 'days since 0001-01-01 00:00:00'
        ds[TIME_STR].attrs['calendar'] = 'noleap'
        actual = numpy_datetime_workaround_encode_cf(ds)

        time_desired = xr.DataArray([days], dims=[TIME_STR])
        desired = xr.Dataset(coords={TIME_STR: time_desired})
        desired[TIME_STR].attrs['units'] = (
            'days since {0}-01-01 00:00:00'.format(979)
        )
        desired[TIME_STR].attrs['calendar'] = 'noleap'

        assert actual.identical(desired)
        self.assertEqual(xr.decode_cf(actual).time.values[0],
                         np.datetime64('1678-02-04'))

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

        # Test sub-monthly intervals
        test = pd.date_range('1999-12-31 18:00:00', '2000-01-01 00:00:00',
                             freq='6H')
        test = xr.DataArray(test, dims=[TIME_STR])
        result_jan = _month_conditional(test, [1])
        np.testing.assert_array_equal(result_jan,
                                      np.array([False, True]))

        # Test wrap-around year
        result_jd = _month_conditional(test, [1, 12])
        np.testing.assert_array_equal(result_jd,
                                      np.array([True, True]))

        # Test month not in range
        result_march = _month_conditional(test, [3])
        np.testing.assert_array_equal(result_march,
                                      np.array([False, False]))

    def test_create_monthly_time_array(self):
        start = datetime.datetime(1, 2, 2)
        end = datetime.datetime(2, 6, 5)
        # using feb and march
        months_bool = [True]*2 + [False]*10 + [True]*2 + [False]*3
        actual = create_monthly_time_array(start, end, 'fm')
        all_months = pd.date_range(
            start=datetime.datetime(pd.Timestamp.min.year + 1, start.month,
                                    start.day),
            end=datetime.datetime(pd.Timestamp.min.year + 2, end.month,
                                  end.day), freq='M'
        )
        desired = xr.DataArray(all_months, dims=[TIME_STR])[months_bool]
        assert actual.identical(desired)

    def test_extract_date_range_and_months(self):
        time = xr.DataArray(pd.date_range(start='2000-02-18', end='2002-07-12',
                                          freq='1D'), dims=[TIME_STR])
        start_date = datetime.datetime(2000, 8, 1)
        end_date = datetime.datetime(2002, 6, 30)
        months = 'mam'  # March-April-May
        desired = xr.concat([
            xr.DataArray(pd.date_range(start='2001-03-01', end='2001-05-31',
                                       freq='1D'), dims=[TIME_STR]),
            xr.DataArray(pd.date_range(start='2002-03-01', end='2002-05-31',
                                       freq='1D'), dims=[TIME_STR])
            ], dim=TIME_STR)
        actual = extract_date_range_and_months(time, start_date, end_date,
                                               months)
        assert actual.identical(desired)

    def test_extract_date_range_and_months_single_month(self):
        time = xr.DataArray(pd.date_range(start='1678-01-01', end='1678-01-31',
                                          freq='1M'), dims=[TIME_STR])
        start_date = datetime.datetime(1678, 1, 1)
        end_date = datetime.datetime(1678, 1, 31)
        months = 1
        desired = time
        actual = extract_date_range_and_months(time, start_date, end_date,
                                               months)
        assert actual.identical(desired)


if __name__ == '__main__':
    sys.exit(unittest.main())
