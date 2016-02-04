#!/usr/bin/env python
"""
02-04-2016 SKC: Set of end-to-end tests of aospy.  Currently set up to test for
runtime errors, but could be extended to test the numerical results.
"""
import unittest

import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime

from aospy.calc import Calc, CalcInterface
from test_objs.projects import aospy_test
from test_objs.models import am2
from test_objs.runs import test_am2
from test_objs.variables import olr


class AospyTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestAospy(AospyTestCase):
    def test_am2_annual_mean(self):
        calc_int = CalcInterface(proj=aospy_test,
                                 model=am2,
                                 run=test_am2,
                                 var=olr,
                                 date_range=('0021-01-01', '0080-12-31'),
                                 intvl_in='monthly',
                                 intvl_out='ann',
                                 dtype_in_time='ts',
                                 dtype_in_vert='pressure',
                                 dtype_out_time='av',
                                 dtype_out_vert=False,
                                 level=False)
        calc = Calc(calc_int)
        calc.compute()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAospy)
    unittest.TextTestRunner(verbosity=2).run(suite)
