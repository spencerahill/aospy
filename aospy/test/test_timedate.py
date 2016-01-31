#!/usr/bin/env python
"""Test suite for aospy.timedate module."""
import unittest

import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime

from aospy.timedate import TimeManager
from aospy import TIME_STR


class AospyTimeManagerTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestTimeManager(AospyTimeManagerTestCase):
    def test_to_datetime_bool(self):
        self.assertEqual(TimeManager.to_datetime(True), True)
        self.assertEqual(TimeManager.to_datetime(False), False)

    def test_to_datetime_datetime(self):
        self.assertEqual(TimeManager.to_datetime(datetime(2000, 1, 1)),
                         datetime(2000, 1, 1))

    def test_to_datetime_year_only(self):
        self.assertEqual(TimeManager.to_datetime(2000), datetime(2000, 1, 1))

    def test_to_datetime_str(self):
        self.assertEqual(TimeManager.to_datetime('2000-01-01'),
                         datetime(2000, 1, 1))

    def test_str_to_datetime(self):
        self.assertEqual(TimeManager.str_to_datetime('2000-01-01'),
                         datetime(2000, 1, 1))
        self.assertEqual(TimeManager.str_to_datetime('2000-01-0199999'),
                         datetime(2000, 1, 1))

    def test_apply_year_offset(self):
        self.assertEqual(TimeManager.apply_year_offset(datetime(1678, 1, 1)),
                         datetime(1678, 1, 1))
        self.assertEqual(TimeManager.apply_year_offset(datetime(1, 1, 1)),
                         datetime(1900, 1, 1))

    def test_construct_month_conditional(self):
        test = pd.date_range('2000-01-01', '2000-03-01', freq='M')
        test = xr.DataArray(test, dims=[TIME_STR], coords=[test])
        result_jan = TimeManager._construct_month_conditional(test, [1])
        np.testing.assert_array_equal(result_jan, np.array([True, False]))

        result_jan_feb = TimeManager._construct_month_conditional(test, [1, 2])
        np.testing.assert_array_equal(result_jan_feb, np.array([True, True]))

        result_march = TimeManager._construct_month_conditional(test, [3])
        np.testing.assert_array_equal(result_march, np.array([False, False]))

        # Test sub-monthly intervals
        test = pd.date_range('1999-12-31 18:00:00', '2000-01-01 00:00:00',
                             freq='6H')
        test = xr.DataArray(test, dims=[TIME_STR])
        result_jan = TimeManager._construct_month_conditional(test, [1])
        np.testing.assert_array_equal(result_jan,
                                      np.array([False, True]))

        # Test wrap-around year
        result_jd = TimeManager._construct_month_conditional(test, [1, 12])
        np.testing.assert_array_equal(result_jd,
                                      np.array([True, True]))

        # Test month not in range
        result_march = TimeManager._construct_month_conditional(test, [3])
        np.testing.assert_array_equal(result_march,
                                      np.array([False, False]))

    def test_month_indices(self):
        np.testing.assert_array_equal(TimeManager.month_indices('ann'),
                                      range(1, 13))
        np.testing.assert_array_equal(TimeManager.month_indices('jja'),
                                      np.array([6, 7, 8]))
        with self.assertRaises(RuntimeError):
            TimeManager.month_indices('dfm')
        with self.assertRaises(RuntimeError):
            TimeManager.month_indices('j')
        with self.assertRaises(RuntimeError):
            TimeManager.month_indices('q')
        self.assertEqual(TimeManager.month_indices('s'), [9])
        np.testing.assert_array_equal(TimeManager.month_indices('djf'),
                                      np.array([12, 1, 2]))
        self.assertEqual(TimeManager.month_indices(12), [12])

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTimeManager)
    unittest.TextTestRunner(verbosity=2).run(suite)
