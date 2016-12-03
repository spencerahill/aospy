#!/usr/bin/env python
"""Test suite for aospy.io module."""
import sys
import unittest

import aospy.utils.io as io


class AospyIOTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestIO(AospyIOTestCase):
    def test_dmget(self):
        # For now this just test to make sure that dmget
        # doesn't raise an exception if the command does not exist
        # on the system.
        io.dmget(['/home/Spencer.Clark/archive/imr_skc/control/'
                  'gfdl.ncrc3-default-repro/1/history/'
                  '00010101.atmos_month.nc'])

if __name__ == '__main__':
    sys.exit(unittest.main())
