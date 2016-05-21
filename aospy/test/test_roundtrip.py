#!/usr/bin/env python
"""
Set of end-to-end tests of aospy.  Currently set up to test for
runtime errors, but could be extended to test the numerical results.
"""
import unittest
from os.path import isfile

from aospy.calc import Calc, CalcInterface
from test_objs.projects import aospy_test
from test_objs.models import (am2, idealized_moist, am3, hiram, cesm1_cam5,
                              idealized_moist_rad)
from test_objs.runs import (test_am2, test_idealized_moist, test_am3,
                            test_hiram, test_amip, test_idealized_moist_rad)
from test_objs.variables import temp, t_surf
from test_objs.regions import nh, sahel, nh_ocean


def model_files_exist(model):
    """Returns True if the grid files specified in the given model
    can be found on the filesystem.

    Parameters
    ----------
    model : Model
        aospy Model object to check

    Returns
    -------
    files_exist : bool
        True if grid files in specified in model exist
    """
    return all([isfile(grid_file)
                for grid_file in model.grid_file_paths])
skip_message = ('Model grid files cannot be located; note this '
                'test can only be completed on the GFDL '
                'filesystems.')


@unittest.skipIf(not model_files_exist(am2), skip_message)
class TestAM2(unittest.TestCase):
    def setUp(self):
        self.two_d_test_params = {'proj': aospy_test,
                                  'model': am2,
                                  'run': test_am2,
                                  'var': t_surf,
                                  'date_range': ('0021-01-01', '0080-12-31'),
                                  'intvl_in': 'monthly',
                                  'dtype_in_time': 'ts',
                                  'dtype_in_vert': 'pressure',
                                  'dtype_out_vert': False,
                                  'level': False}
        self.three_d_test_params = {'proj': aospy_test,
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
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_annual_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='ts',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_seasonal_mean(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='av',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_seasonal_ts(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='ts',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_monthly_mean(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='av',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_monthly_ts(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='ts',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_simple_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'nh': nh},
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_simple_reg_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.ts',
                                 region={'nh': nh},
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_complex_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'sahel': sahel},
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_ocean_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'nh_ocean': nh_ocean},
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_vert_int_pressure(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 dtype_out_vert='vert_int',
                                 dtype_in_vert='pressure',
                                 **self.three_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_vert_int_sigma(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 dtype_out_vert='vert_int',
                                 dtype_in_vert='sigma',
                                 **self.three_d_test_params)
        calc = Calc(calc_int)
        calc.compute()

    def test_sub_daily(self):
        """Not tested for now in models other than idealized_moist.
        and idealized_moist_rad."""
        if not getattr(self, 'sub_daily_test_params', False):
            self.skipTest("Sub-daily test parameters not provided")
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 **self.sub_daily_test_params)
        calc = Calc(calc_int)
        calc.compute()


@unittest.skipIf(not model_files_exist(idealized_moist), skip_message)
class TestIdealizedMoist(TestAM2):
    def setUp(self):
        self.two_d_test_params = {'proj': aospy_test,
                                  'model': idealized_moist,
                                  'run': test_idealized_moist,
                                  'var': t_surf,
                                  'date_range': ('0001-12-27', '0002-12-22'),
                                  'intvl_in': '20-day',
                                  'dtype_in_time': 'ts',
                                  'dtype_in_vert': 'pressure',
                                  'dtype_out_vert': False,
                                  'level': False}
        self.three_d_test_params = {'proj': aospy_test,
                                    'model': idealized_moist,
                                    'run': test_idealized_moist,
                                    'var': temp,
                                    'date_range': ('0001-12-27',
                                                   '0002-12-22'),
                                    'intvl_in': '20-day',
                                    'dtype_in_time': 'ts',
                                    'level': False}
        self.sub_daily_test_params = {'proj': aospy_test,
                                      'model': idealized_moist,
                                      'run': test_idealized_moist,
                                      'var': temp,
                                      'date_range': ('0001-12-27',
                                                     '0002-12-22'),
                                      'intvl_in': '3-hourly',
                                      'dtype_in_time': 'inst',
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


@unittest.skipIf(not model_files_exist(idealized_moist_rad), skip_message)
class TestIdealizedMoistRad(TestIdealizedMoist):
    def setUp(self):
        self.two_d_test_params = {'proj': aospy_test,
                                  'model': idealized_moist_rad,
                                  'run': test_idealized_moist_rad,
                                  'var': t_surf,
                                  'date_range': ('0003-01-01', '0006-12-31'),
                                  'intvl_in': 'monthly',
                                  'dtype_in_time': 'ts',
                                  'dtype_out_vert': False,
                                  'level': False}
        self.three_d_test_params = {'proj': aospy_test,
                                    'model': idealized_moist_rad,
                                    'run': test_idealized_moist_rad,
                                    'var': temp,
                                    'date_range': ('0003-01-01',
                                                   '0006-12-31'),
                                    'intvl_in': 'monthly',
                                    'dtype_in_time': 'ts',
                                    'level': False}
        self.sub_daily_test_params = {'proj': aospy_test,
                                      'model': idealized_moist_rad,
                                      'run': test_idealized_moist_rad,
                                      'var': temp,
                                      'date_range': ('0003-01-01',
                                                     '0006-12-31'),
                                      'intvl_in': '3-hourly',
                                      'dtype_in_time': 'inst',
                                      'level': False}


@unittest.skipIf(not model_files_exist(am3), skip_message)
class TestAM3(TestAM2):
    def setUp(self):
        self.two_d_test_params = {'proj': aospy_test,
                                  'model': am3,
                                  'run': test_am3,
                                  'var': t_surf,
                                  'date_range': ('1981-01-01', '2010-12-31'),
                                  'intvl_in': 'monthly',
                                  'dtype_in_time': 'ts',
                                  'dtype_in_vert': 'pressure',
                                  'dtype_out_vert': False,
                                  'level': False}
        self.three_d_test_params = {'proj': aospy_test,
                                    'model': am3,
                                    'run': test_am3,
                                    'var': temp,
                                    'date_range': ('1981-01-01',
                                                   '2010-12-31'),
                                    'intvl_in': 'monthly',
                                    'dtype_in_time': 'ts',
                                    'level': False}


@unittest.skipIf(not model_files_exist(hiram), skip_message)
class TestHiRAM(TestAM2):
    def setUp(self):
        self.two_d_test_params = {'proj': aospy_test,
                                  'model': hiram,
                                  'run': test_hiram,
                                  'var': t_surf,
                                  'date_range': ('1979-01-01', '1995-12-31'),
                                  'intvl_in': 'monthly',
                                  'dtype_in_time': 'ts',
                                  'dtype_in_vert': 'pressure',
                                  'dtype_out_vert': False,
                                  'level': False}
        self.three_d_test_params = {'proj': aospy_test,
                                    'model': hiram,
                                    'run': test_hiram,
                                    'var': temp,
                                    'date_range': ('1979-01-01',
                                                   '1995-12-31'),
                                    'intvl_in': 'monthly',
                                    'dtype_in_time': 'ts',
                                    'level': False}

    @unittest.skip('no HiRAM output readily at hand exists on sigma levels')
    def test_vert_int_sigma(self):
        pass


# cesm1_cam5 has no grid files so just use hiram's to test to skip
@unittest.skipIf(not model_files_exist(hiram), skip_message)
@unittest.expectedFailure
class TestCMIP5(TestAM2):
    def setUp(self):
        self.two_d_test_params = {'proj': aospy_test,
                                  'model': cesm1_cam5,
                                  'run': test_amip,
                                  'var': t_surf,
                                  'date_range': ('1979-01-01', '2008-12-31'),
                                  'intvl_in': 'monthly',
                                  'dtype_in_time': 'ts',
                                  'dtype_in_vert': 'pressure',
                                  'dtype_out_vert': False,
                                  'level': False}
        self.three_d_test_params = {'proj': aospy_test,
                                    'model': cesm1_cam5,
                                    'run': test_amip,
                                    'var': temp,
                                    'date_range': ('1979-01-01',
                                                   '2008-12-31'),
                                    'intvl_in': 'monthly',
                                    'dtype_in_time': 'ts',
                                    'level': False}

    @unittest.skip(('This repo does not contain any output'
                    ' on sigma levels.'))
    def test_vert_int_sigma(self):
        pass


if __name__ == '__main__':
    unittest.main()
