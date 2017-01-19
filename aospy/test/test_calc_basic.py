#!/usr/bin/env python
"""Basic test of the Calc module on 2D data."""
import datetime
from os.path import isfile
import shutil
import unittest

from aospy.calc import Calc, CalcInterface
from .data.objects.examples import (
    example_proj, example_model, example_run, condensation_rain,
    precip, globe, sahel
)


class TestCalcBasic(unittest.TestCase):
    def setUp(self):
        self.two_d_test_params = {
            'proj': example_proj,
            'model': example_model,
            'run': example_run,
            'var': condensation_rain,
            'date_range': (datetime.datetime(4, 1, 1),
                           datetime.datetime(6, 12, 31)),
            'intvl_in': 'monthly',
            'dtype_in_time': 'ts'
        }

    def tearDown(self):
        for direc in [example_proj.direc_out, example_proj.tar_direc_out]:
            try:
                shutil.rmtree(direc)
            except FileNotFoundError:
                pass

    def test_annual_mean(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='av',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['av'])
        assert isfile(calc.path_tar_out)

    def test_annual_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='ts',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['ts'])
        assert isfile(calc.path_tar_out)

    def test_seasonal_mean(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='av',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['av'])
        assert isfile(calc.path_tar_out)

    def test_seasonal_ts(self):
        calc_int = CalcInterface(intvl_out='djf',
                                 dtype_out_time='ts',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['ts'])
        assert isfile(calc.path_tar_out)

    def test_monthly_mean(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='av',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['av'])
        assert isfile(calc.path_tar_out)

    def test_monthly_ts(self):
        calc_int = CalcInterface(intvl_out=1,
                                 dtype_out_time='ts',
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['ts'])
        assert isfile(calc.path_tar_out)

    def test_simple_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'globe': globe},
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['reg.av'])
        assert isfile(calc.path_tar_out)

    def test_simple_reg_ts(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.ts',
                                 region={'globe': globe},
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['reg.ts'])
        assert isfile(calc.path_tar_out)

    def test_complex_reg_av(self):
        calc_int = CalcInterface(intvl_out='ann',
                                 dtype_out_time='reg.av',
                                 region={'sahel': sahel},
                                 **self.two_d_test_params)
        calc = Calc(calc_int)
        calc.compute()
        assert isfile(calc.path_out['reg.av'])
        assert isfile(calc.path_tar_out)


class TestCalcComposite(TestCalcBasic):
    def setUp(self):
        self.two_d_test_params = {
            'proj': example_proj,
            'model': example_model,
            'run': example_run,
            'var': precip,
            'date_range': (datetime.datetime(4, 1, 1),
                           datetime.datetime(6, 12, 31)),
            'intvl_in': 'monthly',
            'dtype_in_time': 'ts'
        }


if __name__ == '__main__':
    unittest.main()
