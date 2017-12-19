#!/usr/bin/env python
"""Basic test of the Calc module on 2D data."""
import datetime
from os.path import isfile
import shutil
import unittest
import pytest
import itertools

import xarray as xr

from aospy.calc import Calc, _add_metadata_as_attrs
from .data.objects.examples import (
    example_proj, example_model, example_run, var_not_time_defined,
    condensation_rain, precip, sphum, globe, sahel
)


def _test_output_attrs(calc, dtype_out):
    with xr.open_dataset(calc.path_out[dtype_out]) as data:
        expected_units = calc.var.units
        if calc.dtype_out_vert == 'vert_int':
            if expected_units != '':
                expected_units = ("(vertical integral of {0}):"
                                  " {0} m)").format(expected_units)
            else:
                expected_units = ("(vertical integral of quantity"
                                  " with unspecified units)")
        expected_description = calc.var.description
        for name, arr in data.data_vars.items():
            assert expected_units == arr.attrs['units']
            assert expected_description == arr.attrs['description']


def _test_files_and_attrs(calc, dtype_out):
    assert isfile(calc.path_out[dtype_out])
    assert isfile(calc.path_tar_out)
    _test_output_attrs(calc, dtype_out)


class TestCalcBasic(unittest.TestCase):
    def setUp(self):
        self.test_params = {
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
            shutil.rmtree(direc)

    def test_annual_mean(self):
        calc = Calc(intvl_out='ann', dtype_out_time='av', **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'av')

    def test_annual_ts(self):
        calc = Calc(intvl_out='ann', dtype_out_time='ts', **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'ts')

    def test_seasonal_mean(self):
        calc = Calc(intvl_out='djf', dtype_out_time='av', **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'av')

    def test_seasonal_ts(self):
        calc = Calc(intvl_out='djf', dtype_out_time='ts', **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'ts')

    def test_monthly_mean(self):
        calc = Calc(intvl_out=1, dtype_out_time='av', **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'av')

    def test_monthly_ts(self):
        calc = Calc(intvl_out=1, dtype_out_time='ts', **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'ts')

    def test_simple_reg_av(self):
        calc = Calc(intvl_out='ann', dtype_out_time='reg.av', region=[globe],
                    **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'reg.av')

    def test_simple_reg_ts(self):
        calc = Calc(intvl_out='ann', dtype_out_time='reg.ts', region=[globe],
                    **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'reg.ts')

    def test_complex_reg_av(self):
        calc = Calc(intvl_out='ann', dtype_out_time='reg.av', region=[sahel],
                    **self.test_params)
        calc.compute()
        _test_files_and_attrs(calc, 'reg.av')


class TestCalcComposite(TestCalcBasic):
    def setUp(self):
        self.test_params = {
            'proj': example_proj,
            'model': example_model,
            'run': example_run,
            'var': precip,
            'date_range': (datetime.datetime(4, 1, 1),
                           datetime.datetime(6, 12, 31)),
            'intvl_in': 'monthly',
            'dtype_in_time': 'ts'
        }


class TestCalc3D(TestCalcBasic):
    def setUp(self):
        self.test_params = {
            'proj': example_proj,
            'model': example_model,
            'run': example_run,
            'var': sphum,
            'date_range': (datetime.datetime(6, 1, 1),
                           datetime.datetime(6, 1, 31)),
            'intvl_in': 'monthly',
            'dtype_in_time': 'ts',
            'dtype_in_vert': 'sigma',
            'dtype_out_vert': 'vert_int'
        }


test_params = {
    'proj': example_proj,
    'model': example_model,
    'run': example_run,
    'var': var_not_time_defined,
    'date_range': 'default',
    'intvl_in': 'monthly',
    'dtype_in_time': 'av',
    'intvl_out': 1,
}


@pytest.mark.parametrize('dtype_out_time', [None, []])
def test_calc_object_no_time_options(dtype_out_time):
    test_params['dtype_out_time'] = dtype_out_time
    calc = Calc(**test_params)
    if isinstance(dtype_out_time, list):
        assert calc.dtype_out_time == tuple(dtype_out_time)
    else:
        assert calc.dtype_out_time == tuple([dtype_out_time])


@pytest.mark.parametrize(
    'dtype_out_time',
    ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts'])
def test_calc_object_string_time_options(dtype_out_time):
    test_params['dtype_out_time'] = dtype_out_time
    with pytest.raises(ValueError):
        Calc(**test_params)


def test_calc_object_time_options():
    time_options = ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts']
    for i in range(1, len(time_options) + 1):
        for time_option in list(itertools.permutations(time_options, i)):
            if time_option != ('None',):
                test_params['dtype_out_time'] = time_option
                with pytest.raises(ValueError):
                    Calc(**test_params)


@pytest.mark.parametrize(
    ('units', 'description', 'dtype_out_vert', 'expected_units',
     'expected_description'),
    [('', '', None, '', ''),
     ('m', '', None, 'm', ''),
     ('', 'rain', None, '', 'rain'),
     ('m', 'rain', None, 'm', 'rain'),
     ('', '', 'vert_av', '', ''),
     ('m', '', 'vert_av', 'm', ''),
     ('', 'rain', 'vert_av', '', 'rain'),
     ('m', 'rain', 'vert_av', 'm', 'rain'),
     ('', '', 'vert_int',
      '(vertical integral of quantity with unspecified units)', ''),
     ('m', '', 'vert_int',
      '(vertical integral of m): m kg m^-2)', ''),
     ('', 'rain', 'vert_int',
      '(vertical integral of quantity with unspecified units)', 'rain'),
     ('m', 'rain', 'vert_int',
      '(vertical integral of m): m kg m^-2)', 'rain')])
def test_attrs(units, description, dtype_out_vert, expected_units,
               expected_description):
    da = xr.DataArray(None)
    ds = xr.Dataset({'bar': 'foo', 'boo': 'baz'})
    da = _add_metadata_as_attrs(da, units, description, dtype_out_vert)
    ds = _add_metadata_as_attrs(ds, units, description, dtype_out_vert)
    assert expected_units == da.attrs['units']
    assert expected_description == da.attrs['description']
    for name, arr in ds.data_vars.items():
        assert expected_units == arr.attrs['units']
        assert expected_description == arr.attrs['description']


if __name__ == '__main__':
    unittest.main()
