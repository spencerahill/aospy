#!/usr/bin/env python
"""Test suite for aospy.run module."""
import datetime
import sys
import unittest

from aospy.run import Run
from aospy.data_loader import DictDataLoader, GFDLDataLoader


class RunTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestRun(RunTestCase):
    def test_init_dates_valid_input(self):
        for attr in ['default_start_date', 'default_end_date']:
            for date in [None, datetime.datetime(1, 1, 1)]:
                run_ = Run(**{attr: date})
                self.assertEqual(date, getattr(run_, attr))

    def test_init_dates_invalid_input(self):
        for attr in ['default_start_date', 'default_end_date']:
            for date in [1985, False, '1750-12-10']:
                with self.assertRaises(TypeError):
                    Run(**{attr: date})

    def test_init_default_dates(self):
        gdl = GFDLDataLoader(data_start_date=datetime.datetime(1, 1, 1),
                             data_end_date=datetime.datetime(1, 12, 31))
        run_ = Run(data_loader=gdl)
        self.assertEqual(run_.default_start_date, datetime.datetime(1, 1, 1))
        self.assertEqual(run_.default_end_date, datetime.datetime(1, 12, 31))

        ddl = DictDataLoader({'monthly': '/a/'})
        run_ = Run(data_loader=ddl)
        self.assertEqual(run_.default_start_date, None)
        self.assertEqual(run_.default_end_date, None)


if __name__ == '__main__':
    sys.exit(unittest.main())
