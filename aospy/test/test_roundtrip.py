#!/usr/bin/env python
"""
Set of end-to-end tests of aospy.  Currently set up to test for
runtime errors, but could be extended to test the numerical results.
"""
import unittest
from os.path import isfile

from aospy.calc import Calc, CalcInterface
from test_objs.projects import aospy_test
from test_objs.models import am2
from test_objs.runs import test_am2
from test_objs.variables import olr, temp
from test_objs.regions import nh, sahel, nh_ocean


class AospyTestCase(unittest.TestCase):
    def setUp(self):
        self.am2_olr_test_params = {'proj': aospy_test,
                                    'model': am2,
                                    'run': test_am2,
                                    'var': olr,
                                    'date_range': ('0021-01-01', '0080-12-31'),
                                    'intvl_in': 'monthly',
                                    'dtype_in_time': 'ts',
                                    'dtype_in_vert': 'pressure',
                                    'dtype_out_vert': False,
                                    'level': False}
        self.am2_temp_test_params = {'proj': aospy_test,
                                     'model': am2,
                                     'run': test_am2,
                                     'var': temp,
                                     'date_range': ('0021-01-01',
                                                    '0080-12-31'),
                                     'intvl_in': 'monthly',
                                     'dtype_in_time': 'ts',
                                     'level': False}

    def tearDown(self):
        pass


class TestAospy(AospyTestCase):
    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_annual_mean(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_annual_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='ts',
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_seasonal_mean(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='av',
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_seasonal_ts(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='ts',
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_monthly_mean(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='av',
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_monthly_ts(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='ts',
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_simple_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'nh': nh},
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_simple_reg_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.ts',
                                 region={'nh': nh},
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_complex_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'sahel': sahel},
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_ocean_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'nh_ocean': nh_ocean},
                                 **self.am2_olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_vert_int_pressure(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 dtype_out_vert='vert_int',
                                 dtype_in_vert='pressure',
                                 **self.am2_temp_test_params)
        calc = Calc(calc_int)
        calc.compute()

    @unittest.skipIf(not all([isfile(grid_file)
                              for grid_file in am2.grid_file_paths]),
                     'Model grid files cannot be located; note this '
                     'test can only be completed on the GFDL '
                     'filesystems.')
    def test_am2_vert_int_sigma(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 dtype_out_vert='vert_int',
                                 dtype_in_vert='sigma',
                                 **self.am2_temp_test_params)
        calc = Calc(calc_int)
        calc.compute()


if __name__ == '__main__':
    unittest.main()
