#!/usr/bin/env python
"""Test suite for aospy.timedate module."""
import sys
import unittest

import numpy as np
import xarray as xr
from datetime import datetime

from aospy.timedate import TimeManager


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

    def test_apply_year_offset(self):
        self.assertEqual(TimeManager.apply_year_offset(datetime(1678, 1, 1)),
                         datetime(1678, 1, 1))
        self.assertEqual(TimeManager.apply_year_offset(datetime(1, 1, 1)),
                         datetime(1900, 1, 1))

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTimeManager)
    unittest.TextTestRunner(verbosity=2).run(suite)
