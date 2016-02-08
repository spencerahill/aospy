#!/usr/bin/env python
"""
Set of end-to-end tests of aospy.  Currently set up to test for
runtime errors, but could be extended to test the numerical results.
"""
import unittest
from os.path import isfile

from aospy.calc import Calc, CalcInterface
from test_objs.projects import aospy_test
from test_objs.models import am2, idealized_moist
from test_objs.runs import test_am2, test_idealized_moist
from test_objs.variables import olr, temp
from test_objs.regions import nh, sahel, nh_ocean

am2_files_exist = all([isfile(grid_file)
                       for grid_file in am2.grid_file_paths])
idealized_files_exist = all([isfile(grid_file)
                             for grid_file in idealized_moist.grid_file_paths])
skip_message = ('Model grid files cannot be located; note this '
                'test can only be completed on the GFDL '
                'filesystems.')


@unittest.skipIf(not am2_files_exist, skip_message)
class TestAM2(unittest.TestCase):
    def setUp(self):
        self.olr_test_params = {'proj': aospy_test,
                                'model': am2,
                                'run': test_am2,
                                'var': olr,
                                'date_range': ('0021-01-01', '0080-12-31'),
                                'intvl_in': 'monthly',
                                'dtype_in_time': 'ts',
                                'dtype_in_vert': 'pressure',
                                'dtype_out_vert': False,
                                'level': False}
        self.temp_test_params = {'proj': aospy_test,
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

    def test_annual_mean(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_annual_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='ts',
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_seasonal_mean(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='av',
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_seasonal_ts(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='ts',
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_monthly_mean(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='av',
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_monthly_ts(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='ts',
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_simple_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'nh': nh},
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_simple_reg_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.ts',
                                 region={'nh': nh},
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_complex_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'sahel': sahel},
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_ocean_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'nh_ocean': nh_ocean},
                                 **self.olr_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_vert_int_pressure(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 dtype_out_vert='vert_int',
                                 dtype_in_vert='pressure',
                                 **self.temp_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_vert_int_sigma(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 dtype_out_vert='vert_int',
                                 dtype_in_vert='sigma',
                                 **self.temp_test_params)
        calc = Calc(calc_int)
        calc.compute()


@unittest.skipIf(not idealized_files_exist, skip_message)
class TestIdealized(TestAM2):
    def setUp(self):
        self.olr_test_params = {'proj': aospy_test,
                                'model': idealized_moist,
                                'run': test_idealized_moist,
                                'var': olr,
                                'date_range': ('0001-12-27', '0002-12-22'),
                                'intvl_in': '20-day',
                                'dtype_in_time': 'ts',
                                'dtype_in_vert': 'pressure',
                                'dtype_out_vert': False,
                                'level': False}
        self.temp_test_params = {'proj': aospy_test,
                                 'model': idealized_moist,
                                 'run': test_idealized_moist,
                                 'var': temp,
                                 'date_range': ('0001-12-27',
                                                '0002-12-22'),
                                 'intvl_in': '20-day',
                                 'dtype_in_time': 'ts',
                                 'level': False}

    @unittest.skip('not valid in idealized moist model')
    def test_seasonal_mean(self):
        pass

    @unittest.skip('not valid in idealized moist model')
    def test_seasonal_ts(self):
        pass

    @unittest.skip('not valid in idealized moist model')
    def test_monthly_mean(self):
        pass

    @unittest.skip('not valid in idealized moist model')
    def test_monthly_ts(self):
        pass

    @unittest.skip('not valid in idealized moist model')
    def test_vert_int_pressure(self):
        pass

if __name__ == '__main__':
    unittest.main()
