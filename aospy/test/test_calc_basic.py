#!/usr/bin/env python
"""Basic test of the Calc module on 2D data."""
import datetime
from os.path import isfile
import shutil
import unittest
import pytest
import itertools

import cftime
import numpy as np
import xarray as xr

from aospy import Var
from aospy.calc import Calc, _add_metadata_as_attrs, _replace_pressure
from aospy.internal_names import ETA_STR
from aospy.utils.vertcoord import p_eta, dp_eta, p_level, dp_level
from .data.objects.examples import (
    example_proj, example_model, example_run, var_not_time_defined,
    condensation_rain, convection_rain, precip, sphum, globe, sahel, p, dp
)


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def _test_output_attrs(calc, dtype_out):
    with xr.set_options(enable_cftimeindex=True):
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


def _clean_test_direcs():
    for direc in [example_proj.direc_out, example_proj.tar_direc_out]:
        try:
            shutil.rmtree(direc)
        except OSError:
            pass


def _test_files_and_attrs(calc, dtype_out):
    assert isfile(calc.path_out[dtype_out])
    assert isfile(calc.path_tar_out)
    _test_output_attrs(calc, dtype_out)


_2D_DATE_RANGES = {
    'datetime': (datetime.datetime(4, 1, 1), datetime.datetime(6, 12, 31)),
    'datetime64': (np.datetime64('0004-01-01'), np.datetime64('0006-12-31')),
    'cftime': (cftime.DatetimeNoLeap(4, 1, 1),
               cftime.DatetimeNoLeap(6, 12, 31)),
    'str': ('0004', '0006')
}
_3D_DATE_RANGES = {
    'datetime': (datetime.datetime(6, 1, 1), datetime.datetime(6, 1, 31)),
    'datetime64': (np.datetime64('0006-01-01'), np.datetime64('0006-01-31')),
    'cftime': (cftime.DatetimeNoLeap(6, 1, 1),
               cftime.DatetimeNoLeap(6, 1, 31)),
    'str': ('0006', '0006')
}
_2D_VARS = {'basic': condensation_rain, 'composite': precip}
_2D_DTYPE_OUT_VERT = {'None': None}
_2D_DTYPE_IN_VERT = {'None': None}
_3D_VARS = {'3D': sphum}
_3D_DTYPE_OUT_VERT = {'vert_int': 'vert_int',
                      'vert_av': 'vert_av'}
_3D_DTYPE_IN_VERT = {'sigma': 'sigma'}
_CASES = (
    list(itertools.product(_2D_DATE_RANGES.items(), _2D_VARS.items(),
                           _2D_DTYPE_IN_VERT.items(),
                           _2D_DTYPE_OUT_VERT.items())) +
    list(itertools.product(_3D_DATE_RANGES.items(), _3D_VARS.items(),
                           _3D_DTYPE_IN_VERT.items(),
                           _3D_DTYPE_OUT_VERT.items()))
)
_CALC_TESTS = {}
for ((date_type, date_range), (test_type, var),
     (vert_in_label, vert_in), (vert_out_label, vert_out)) in _CASES:
    _CALC_TESTS['{}-{}-{}-{}'.format(
        date_type, test_type, vert_in_label, vert_out_label)] = (
            date_range, var, vert_in, vert_out)


@pytest.fixture(params=_CALC_TESTS.values(), ids=list(_CALC_TESTS.keys()))
def test_params(request):
    date_range, var, vert_in, vert_out = request.param
    yield {
        'proj': example_proj,
        'model': example_model,
        'run': example_run,
        'var': var,
        'date_range': date_range,
        'intvl_in': 'monthly',
        'dtype_in_time': 'ts',
        'dtype_in_vert': vert_in,
        'dtype_out_vert': vert_out
    }
    _clean_test_direcs()


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_annual_mean(test_params):
    calc = Calc(intvl_out='ann', dtype_out_time='av', **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'av')


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_annual_ts(test_params):
    calc = Calc(intvl_out='ann', dtype_out_time='ts', **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'ts')


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_seasonal_mean(test_params):
    calc = Calc(intvl_out='djf', dtype_out_time='av', **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'av')


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_seasonal_ts(test_params):
    calc = Calc(intvl_out='djf', dtype_out_time='ts', **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'ts')


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_monthly_mean(test_params):
    calc = Calc(intvl_out=1, dtype_out_time='av', **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'av')


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_monthly_ts(test_params):
    calc = Calc(intvl_out=1, dtype_out_time='ts', **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'ts')


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_simple_reg_av(test_params):
    calc = Calc(intvl_out='ann', dtype_out_time='reg.av', region=[globe],
                **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'reg.av')


@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_simple_reg_ts(test_params):
    calc = Calc(intvl_out='ann', dtype_out_time='reg.ts', region=[globe],
                **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'reg.ts')


@pytest.mark.filterwarnings('ignore:Mean of empty slice')
@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_complex_reg_av(test_params):
    calc = Calc(intvl_out='ann', dtype_out_time='reg.av', region=[sahel],
                **test_params)
    calc.compute()
    _test_files_and_attrs(calc, 'reg.av')


test_params_not_time_defined = {
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
    test_params_not_time_defined['dtype_out_time'] = dtype_out_time
    calc = Calc(**test_params_not_time_defined)
    if isinstance(dtype_out_time, list):
        assert calc.dtype_out_time == tuple(dtype_out_time)
    else:
        assert calc.dtype_out_time == tuple([dtype_out_time])


@pytest.mark.parametrize(
    'dtype_out_time',
    ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts'])
def test_calc_object_string_time_options(dtype_out_time):
    test_params_not_time_defined['dtype_out_time'] = dtype_out_time
    with pytest.raises(ValueError):
        Calc(**test_params_not_time_defined)


def test_calc_object_time_options():
    time_options = ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts']
    for i in range(1, len(time_options) + 1):
        for time_option in list(itertools.permutations(time_options, i)):
            if time_option != ('None',):
                test_params_not_time_defined['dtype_out_time'] = time_option
                with pytest.raises(ValueError):
                    Calc(**test_params_not_time_defined)


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


@pytest.fixture()
def recursive_test_params():
    basic_params = {
        'proj': example_proj,
        'model': example_model,
        'run': example_run,
        'var': condensation_rain,
        'date_range': (datetime.datetime(4, 1, 1),
                       datetime.datetime(6, 12, 31)),
        'intvl_in': 'monthly',
        'dtype_in_time': 'ts'
    }
    recursive_params = basic_params.copy()

    recursive_condensation_rain = Var(
        name='recursive_condensation_rain',
        variables=(precip, convection_rain), func=lambda x, y: x - y,
        def_time=True)
    recursive_params['var'] = recursive_condensation_rain

    yield (basic_params, recursive_params)

    _clean_test_direcs()


# Remove autoclose argument in xarray version 0.12
@pytest.mark.filterwarnings('ignore:The autoclose')
@pytest.mark.filterwarnings('ignore:The enable_cftimeindex')
def test_recursive_calculation(recursive_test_params):
    basic_params, recursive_params = recursive_test_params

    calc = Calc(intvl_out='ann', dtype_out_time='av', **basic_params)
    calc = calc.compute()
    with xr.set_options(enable_cftimeindex=True):
        expected = xr.open_dataset(
            calc.path_out['av'], autoclose=True)['condensation_rain']
    _test_files_and_attrs(calc, 'av')

    calc = Calc(intvl_out='ann', dtype_out_time='av', **recursive_params)
    calc = calc.compute()
    with xr.set_options(enable_cftimeindex=True):
        result = xr.open_dataset(
            calc.path_out['av'], autoclose=True)['recursive_condensation_rain']
    _test_files_and_attrs(calc, 'av')

    xr.testing.assert_equal(expected, result)


def test_compute_pressure():
    calc = Calc(
        intvl_out='ann',
        dtype_out_time='av',
        var=p,
        proj=example_proj,
        model=example_model,
        run=example_run,
        date_range=('0006', '0006'),
        intvl_in='monthly',
        dtype_in_time='ts',
        dtype_in_vert='sigma',
        dtype_out_vert=None
    )
    calc.compute()
    _test_files_and_attrs(calc, 'av')
    _clean_test_direcs()


def test_compute_pressure_thicknesses():
    calc = Calc(
        intvl_out='ann',
        dtype_out_time='av',
        var=dp,
        proj=example_proj,
        model=example_model,
        run=example_run,
        date_range=('0006', '0006'),
        intvl_in='monthly',
        dtype_in_time='ts',
        dtype_in_vert='sigma',
        dtype_out_vert=None
    )
    calc.compute()
    _test_files_and_attrs(calc, 'av')
    _clean_test_direcs()


@pytest.mark.parametrize(
    ['dtype_in_vert', 'expected'],
    [(ETA_STR, [p_eta, dp_eta, condensation_rain, 5]),
     ('pressure', [p_level, dp_level, condensation_rain, 5])])
def test_replace_pressure(dtype_in_vert, expected):
    arguments = [p, dp, condensation_rain, 5]
    p_in, dp_in, cond, num = arguments
    p_expected, dp_expected, cond_expected, num_expected = expected
    assert p_in.func != p_expected.func
    assert dp_in.func != dp_expected.func
    result = _replace_pressure(arguments, dtype_in_vert)
    p_out, dp_out, cond_out, num_out = result
    assert p_out.func == p_expected.func
    assert dp_out.func == dp_expected.func
    assert cond_out.func == cond_expected.func
    assert num_out == num_expected


if __name__ == '__main__':
    unittest.main()
