#!/usr/bin/env python
"""Test suite for aospy.data_loader module."""
import datetime
import os
import unittest
import warnings

import numpy as np
import pytest
import xarray as xr

from cftime import DatetimeNoLeap

from aospy import Var
from aospy.data_loader import (DataLoader, DictDataLoader, GFDLDataLoader,
                               NestedDictDataLoader, grid_attrs_to_aospy_names,
                               set_grid_attrs_as_coords, _sel_var,
                               _prep_time_data,
                               _preprocess_and_rename_grid_attrs,
                               _maybe_cast_to_float64)
from aospy.internal_names import (LAT_STR, LON_STR, TIME_STR, TIME_BOUNDS_STR,
                                  BOUNDS_STR, SFC_AREA_STR, ETA_STR, PHALF_STR,
                                  TIME_WEIGHTS_STR, GRID_ATTRS, ZSURF_STR)
from aospy.utils import io
from .data.objects.examples import (condensation_rain, convection_rain, precip,
                                    file_map, ROOT_PATH, example_model, bk)


def _open_ds_catch_warnings(path):
    with warnings.catch_warnings(record=True):
        return xr.open_dataset(path)


@pytest.mark.parametrize(
    ('input_dtype', 'expected_dtype'),
    [(np.float32, np.float64),
     (np.int, np.int)])
def test_maybe_cast_to_float64(input_dtype, expected_dtype):
    da = xr.DataArray(np.ones(3, dtype=input_dtype))
    result = _maybe_cast_to_float64(da).dtype
    assert result == expected_dtype


_DATE_RANGES = {
    'datetime': (datetime.datetime(2000, 1, 1),
                 datetime.datetime(2002, 12, 31)),
    'datetime64': (np.datetime64('2000-01-01'),
                   np.datetime64('2002-12-31')),
    'cftime': (DatetimeNoLeap(2000, 1, 1),
               DatetimeNoLeap(2002, 12, 31)),
    'str': ('2000', '2002')
}


@pytest.fixture(params=_DATE_RANGES.values(), ids=list(_DATE_RANGES.keys()))
def generate_file_set_args(request):
    start_date, end_date = request.param
    return dict(
        var=condensation_rain, start_date=start_date, end_date=end_date,
        domain='atmos', intvl_in='monthly', dtype_in_vert='sigma',
        dtype_in_time='ts', intvl_out=None)


@pytest.fixture()
def alt_lat_str():
    return 'LATITUDE'


@pytest.fixture()
def var_name():
    return 'a'


@pytest.fixture()
def ds(alt_lat_str, var_name):
    time_bounds = np.array([[0, 31], [31, 59], [59, 90]])
    bounds = np.array([0, 1])
    time = np.array([15, 46, 74])
    data = np.zeros((3, 1, 1))
    lat = [0]
    lon = [0]
    ds = xr.DataArray(data,
                      coords=[time, lat, lon],
                      dims=[TIME_STR, alt_lat_str, LON_STR],
                      name=var_name).to_dataset()
    ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                       coords=[time, bounds],
                                       dims=[TIME_STR, BOUNDS_STR],
                                       name=TIME_BOUNDS_STR)
    units_str = 'days since 2000-01-01 00:00:00'
    ds[TIME_STR].attrs['units'] = units_str
    ds[TIME_BOUNDS_STR].attrs['units'] = units_str
    return ds


@pytest.fixture()
def inst_ds(ds):
    inst_time = np.array([3, 6, 9])
    inst_units_str = 'hours since 2000-01-01 00:00:00'
    inst_ds = ds.copy()
    inst_ds.drop(TIME_BOUNDS_STR)
    inst_ds[TIME_STR].values = inst_time
    inst_ds[TIME_STR].attrs['units'] = inst_units_str
    return inst_ds


def _gfdl_data_loader_kwargs(data_start_date, data_end_date):
    return dict(data_direc=os.path.join('.', 'test'),
                data_dur=6,
                data_start_date=data_start_date,
                data_end_date=data_end_date,
                upcast_float32=False,
                data_vars='minimal',
                coords='minimal')


_DATA_LOADER_KWARGS = {
    'DataLoader': (DataLoader, {}),
    'DictDataLoader': (DictDataLoader, dict(file_map={'monthly': ['a.nc']})),
    'NestedDictDataLoader': (
        NestedDictDataLoader,
        dict(file_map={'monthly': {'condensation_rain': ['a.nc']}})),
    'GFDLDataLoader-datetime': (
        GFDLDataLoader, _gfdl_data_loader_kwargs(
            datetime.datetime(2000, 1, 1), datetime.datetime(2012, 12, 31))),
    'GFDLDataLoader-datetime64': (
        GFDLDataLoader, _gfdl_data_loader_kwargs(
            np.datetime64('2000-01-01'), np.datetime64('2012-12-31'))),
    'GFDLDataLoader-cftime': (
        GFDLDataLoader, _gfdl_data_loader_kwargs(
            DatetimeNoLeap(2000, 1, 1), DatetimeNoLeap(2012, 12, 31))),
    'GFDLDataLoader-str': (
        GFDLDataLoader, _gfdl_data_loader_kwargs('2000', '2012'))
}


@pytest.fixture(params=_DATA_LOADER_KWARGS.values(),
                ids=list(_DATA_LOADER_KWARGS.keys()))
def data_loader(request):
    data_loader_type, kwargs = request.param
    return data_loader_type(**kwargs)


_GFDL_DATA_LOADER_KWARGS = {key: _DATA_LOADER_KWARGS[key] for key in
                            _DATA_LOADER_KWARGS if 'GFDL' in key}


@pytest.fixture(params=_GFDL_DATA_LOADER_KWARGS.values(),
                ids=list(_GFDL_DATA_LOADER_KWARGS.keys()))
def gfdl_data_loader(request):
    data_loader_type, kwargs = request.param
    return data_loader_type(**kwargs)


def test_rename_grid_attrs_ds(ds, alt_lat_str):
    assert LAT_STR not in ds
    assert alt_lat_str in ds
    ds = grid_attrs_to_aospy_names(ds)
    assert LAT_STR in ds


def test_rename_grid_attrs_dim_no_coord(ds, var_name):
    bounds_dim = 'nv'
    assert bounds_dim not in ds
    assert bounds_dim in GRID_ATTRS[BOUNDS_STR]
    # Create DataArray with all dims lacking coords
    values = ds[var_name].values
    arr = xr.DataArray(values, name='dummy')
    # Insert name to be replaced (its physical meaning doesn't matter here)
    ds = arr.rename({'dim_0': bounds_dim}).to_dataset()
    assert not ds[bounds_dim].coords
    result = grid_attrs_to_aospy_names(ds)
    assert not result[BOUNDS_STR].coords


def test_rename_grid_attrs_skip_scalar_dim(ds):
    phalf_dim = 'phalf'
    assert phalf_dim not in ds
    assert phalf_dim in GRID_ATTRS[PHALF_STR]
    ds_copy = ds.copy()
    ds_copy[phalf_dim] = 4
    ds_copy = ds_copy.set_coords(phalf_dim)
    result = grid_attrs_to_aospy_names(ds_copy)
    xr.testing.assert_identical(result[phalf_dim], ds_copy[phalf_dim])


def test_rename_grid_attrs_copy_attrs(ds, alt_lat_str):
    orig_attrs = {'dummy_key': 'dummy_val'}
    ds_orig = ds.copy()
    ds_orig[alt_lat_str].attrs = orig_attrs
    ds = grid_attrs_to_aospy_names(ds_orig)
    assert ds[LAT_STR].attrs == orig_attrs


def test_rename_grid_attrs_custom(ds, alt_lat_str):
    assert LAT_STR not in ds
    ds = ds.rename({alt_lat_str: 'custom_lat_name'})
    ds = grid_attrs_to_aospy_names(ds, {LAT_STR: 'custom_lat_name'})
    assert LAT_STR in ds
    assert 'custom_lat_name' not in ds


def test_rename_grid_attrs_custom_error(ds, alt_lat_str):
    assert LAT_STR not in ds
    ds = ds.rename({alt_lat_str: 'custom_lat_name'})
    with pytest.raises(ValueError):
        ds = grid_attrs_to_aospy_names(ds, {alt_lat_str: 'custom_lat_name'})


def test_set_grid_attrs_as_coords(ds, var_name):
    ds = grid_attrs_to_aospy_names(ds)
    sfc_area = ds[var_name].isel(**{TIME_STR: 0}).drop(TIME_STR)
    ds[SFC_AREA_STR] = sfc_area

    assert SFC_AREA_STR not in ds.coords

    ds = set_grid_attrs_as_coords(ds)
    assert SFC_AREA_STR in ds.coords
    assert TIME_BOUNDS_STR in ds.coords


def test_sel_var():
    time = np.array([0, 31, 59]) + 15
    data = np.zeros((3))
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name=convection_rain.name).to_dataset()
    condensation_rain_alt_name, = condensation_rain.alt_names
    ds[condensation_rain_alt_name] = xr.DataArray(data, coords=[ds.time])
    result = _sel_var(ds, convection_rain)
    assert result.name == convection_rain.name

    result = _sel_var(ds, condensation_rain)
    assert result.name == condensation_rain.name

    with pytest.raises(LookupError):
        _sel_var(ds, precip)


def test_maybe_apply_time_shift(data_loader, ds, inst_ds, var_name,
                                generate_file_set_args):
    ds = xr.decode_cf(ds)
    da = ds[var_name]

    result = data_loader._maybe_apply_time_shift(
        da.copy(), **generate_file_set_args)[TIME_STR]
    assert result.identical(da[TIME_STR])

    offset = data_loader._maybe_apply_time_shift(
        da.copy(), {'days': 1}, **generate_file_set_args)
    result = offset[TIME_STR]

    expected = da[TIME_STR] + np.timedelta64(1, 'D')
    expected[TIME_STR] = expected

    assert result.identical(expected)


def test_maybe_apply_time_shift_ts(gfdl_data_loader, ds, var_name,
                                   generate_file_set_args):
    ds = xr.decode_cf(ds)
    da = ds[var_name]
    result = gfdl_data_loader._maybe_apply_time_shift(
        da.copy(), **generate_file_set_args)[TIME_STR]
    assert result.identical(da[TIME_STR])


def test_maybe_apply_time_shift_inst(gfdl_data_loader, inst_ds, var_name,
                                     generate_file_set_args):
    inst_ds = xr.decode_cf(inst_ds)
    generate_file_set_args['dtype_in_time'] = 'inst'
    generate_file_set_args['intvl_in'] = '3hr'
    da = inst_ds[var_name]
    result = gfdl_data_loader._maybe_apply_time_shift(
        da.copy(), **generate_file_set_args)[TIME_STR]

    expected = da[TIME_STR] + np.timedelta64(-3, 'h')
    expected[TIME_STR] = expected
    assert result.identical(expected)

    generate_file_set_args['intvl_in'] = 'daily'
    da = inst_ds[var_name]
    result = gfdl_data_loader._maybe_apply_time_shift(
        da.copy(), **generate_file_set_args)[TIME_STR]

    expected = da[TIME_STR]
    expected[TIME_STR] = expected
    assert result.identical(expected)


def test_prep_time_data(inst_ds):
    assert (TIME_WEIGHTS_STR not in inst_ds)
    ds = _prep_time_data(inst_ds)
    assert (TIME_WEIGHTS_STR in ds)


def test_preprocess_and_rename_grid_attrs(ds, alt_lat_str):
    def preprocess_func(ds, **kwargs):
        # Corrupt a grid attribute name so that we test
        # that grid_attrs_to_aospy_names is still called
        # after
        ds = ds.rename({LON_STR: 'LONGITUDE'})
        ds.attrs['a'] = 'b'
        return ds

    assert LAT_STR not in ds
    assert alt_lat_str in ds
    assert LON_STR in ds

    expected = ds.rename({alt_lat_str: LAT_STR})
    expected = expected.set_coords(TIME_BOUNDS_STR)
    expected.attrs['a'] = 'b'
    result = _preprocess_and_rename_grid_attrs(preprocess_func)(ds)
    xr.testing.assert_identical(result, expected)


def test_generate_file_set(data_loader, generate_file_set_args):
    if type(data_loader) is DataLoader:
        with pytest.raises(NotImplementedError):
            data_loader._generate_file_set()

    elif isinstance(data_loader, DictDataLoader):
        result = data_loader._generate_file_set(
            **generate_file_set_args)
        expected = ['a.nc']
        result == expected

        with pytest.raises(KeyError):
            generate_file_set_args['intvl_in'] = 'daily'
            result = data_loader._generate_file_set(
                **generate_file_set_args)

    elif isinstance(data_loader, NestedDictDataLoader):
        result = data_loader._generate_file_set(
            **generate_file_set_args)
        expected = ['a.nc']
        assert result == expected

        with pytest.raises(KeyError):
            generate_file_set_args['var'] = convection_rain
            result = data_loader._generate_file_set(
                **generate_file_set_args)

    else:
        with pytest.raises(IOError):
            data_loader._generate_file_set(**generate_file_set_args)


def test_overriding_constructor(gfdl_data_loader, ds):
    new = GFDLDataLoader(gfdl_data_loader,
                         data_direc=os.path.join('.', 'a'))
    assert new.data_direc == os.path.join('.', 'a')
    assert new.data_dur == gfdl_data_loader.data_dur
    assert new.data_start_date == gfdl_data_loader.data_start_date
    assert new.data_end_date == gfdl_data_loader.data_end_date
    assert new.preprocess_func == gfdl_data_loader.preprocess_func
    assert new.upcast_float32 == gfdl_data_loader.upcast_float32

    new = GFDLDataLoader(gfdl_data_loader, data_dur=8)
    assert new.data_dur == 8

    new = GFDLDataLoader(gfdl_data_loader,
                         data_start_date=datetime.datetime(2001, 1, 1))
    assert new.data_start_date == datetime.datetime(2001, 1, 1)

    new = GFDLDataLoader(gfdl_data_loader,
                         data_end_date=datetime.datetime(2003, 12, 31))
    assert new.data_end_date == datetime.datetime(2003, 12, 31)

    new = GFDLDataLoader(gfdl_data_loader, preprocess_func=lambda ds: ds)
    xr.testing.assert_identical(new.preprocess_func(ds), ds)

    new = GFDLDataLoader(gfdl_data_loader, upcast_float32=True)
    assert new.upcast_float32

    new = GFDLDataLoader(gfdl_data_loader, data_vars='all')
    assert new.data_vars == 'all'

    new = GFDLDataLoader(gfdl_data_loader, coords='all')
    assert new.coords == 'all'


_GFDL_DATE_RANGES = {
    'datetime': (datetime.datetime(2010, 1, 1),
                 datetime.datetime(2010, 12, 31)),
    'datetime64': (np.datetime64('2010-01-01'),
                   np.datetime64('2010-12-31')),
    'cftime': (DatetimeNoLeap(2010, 1, 1),
               DatetimeNoLeap(2010, 12, 31)),
    'str': ('2010', '2010')
}


@pytest.mark.parametrize(['start_date', 'end_date'],
                         _GFDL_DATE_RANGES.values(),
                         ids=list(_GFDL_DATE_RANGES.keys()))
def test_input_data_paths_gfdl(gfdl_data_loader, start_date, end_date):
    expected = [os.path.join('.', 'test', 'atmos', 'ts', 'monthly', '6yr',
                             'atmos.200601-201112.temp.nc')]
    result = gfdl_data_loader._input_data_paths_gfdl(
        'temp', start_date, end_date, 'atmos',
        'monthly', 'pressure', 'ts', None)
    assert result == expected

    expected = [os.path.join('.', 'test', 'atmos_daily', 'ts', 'daily',
                             '6yr',
                             'atmos_daily.20060101-20111231.temp.nc')]
    result = gfdl_data_loader._input_data_paths_gfdl(
        'temp', start_date, end_date, 'atmos',
        'daily', 'pressure', 'ts', None)
    assert result == expected

    expected = [os.path.join('.', 'test', 'atmos_daily', 'ts', 'daily',
                             '6yr',
                             'atmos_daily.20060101-20111231.temp.nc')]
    result = gfdl_data_loader._input_data_paths_gfdl(
        'temp', start_date, end_date, 'atmos',
        'daily', 'pressure', 'ts', None)
    assert result == expected

    expected = [os.path.join('.', 'test', 'atmos_level', 'ts', 'monthly',
                             '6yr', 'atmos_level.200601-201112.temp.nc')]
    result = gfdl_data_loader._input_data_paths_gfdl(
        'temp', start_date, end_date, 'atmos',
        'monthly', ETA_STR, 'ts', None)
    assert result == expected

    expected = [os.path.join('.', 'test', 'atmos', 'ts', 'monthly',
                             '6yr', 'atmos.200601-201112.ps.nc')]
    result = gfdl_data_loader._input_data_paths_gfdl(
        'ps', start_date, end_date, 'atmos',
        'monthly', ETA_STR, 'ts', None)
    assert result == expected

    expected = [os.path.join('.', 'test', 'atmos_inst', 'ts', 'monthly',
                             '6yr', 'atmos_inst.200601-201112.temp.nc')]
    result = gfdl_data_loader._input_data_paths_gfdl(
        'temp', start_date, end_date, 'atmos',
        'monthly', 'pressure', 'inst', None)
    assert result == expected

    expected = [os.path.join('.', 'test', 'atmos', 'av', 'monthly_6yr',
                             'atmos.2006-2011.jja.nc')]
    result = gfdl_data_loader._input_data_paths_gfdl(
        'temp', start_date, end_date, 'atmos',
        'monthly', 'pressure', 'av', 'jja')
    assert result == expected


# TODO: Parametrize these tests
def test_data_name_gfdl_annual():
    for data_type in ['ts', 'inst']:
        expected = 'atmos.2010.temp.nc'
        result = io.data_name_gfdl('temp', 'atmos', data_type,
                                   'annual', 2010, None, 2000, 1)
        assert result == expected

        expected = 'atmos.2006-2011.temp.nc'
        result = io.data_name_gfdl('temp', 'atmos', data_type,
                                   'annual', 2010, None, 2000, 6)
        assert result == expected

    for intvl_type in ['annual', 'ann']:
        expected = 'atmos.2010.ann.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av',
                                   intvl_type, 2010, None, 2000, 1)
        assert result == expected

        expected = 'atmos.2006-2011.ann.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av',
                                   intvl_type, 2010, None, 2000, 6)
        assert result == expected

    expected = 'atmos.2006-2011.01-12.nc'
    result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                               'annual', 2010, None, 2000, 6)
    assert result == expected


def test_data_name_gfdl_monthly():
    for data_type in ['ts', 'inst']:
        expected = 'atmos.200601-201112.temp.nc'
        result = io.data_name_gfdl('temp', 'atmos', data_type,
                                   'monthly', 2010, 'jja', 2000, 6)
        assert result == expected

    for intvl_type in ['monthly', 'mon']:
        expected = 'atmos.2010.jja.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av',
                                   intvl_type, 2010, 'jja', 2000, 1)
        assert result == expected

        expected = 'atmos.2006-2011.jja.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av',
                                   intvl_type, 2010, 'jja', 2000, 6)
        assert result == expected

    expected = 'atmos.2006-2011.01-12.nc'
    result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                               'monthly', 2010, 'jja', 2000, 6)
    assert result == expected


def test_data_name_gfdl_daily():
    for data_type in ['ts', 'inst']:
        expected = 'atmos.20060101-20111231.temp.nc'
        result = io.data_name_gfdl('temp', 'atmos', data_type,
                                   'daily', 2010, None, 2000, 6)
        assert result == expected

    with pytest.raises(NameError):
        io.data_name_gfdl('temp', 'atmos', 'av',
                          'daily', 2010, None, 2000, 6)

    expected = 'atmos.2006-2011.01-12.nc'
    result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                               'daily', 2010, None, 2000, 6)
    assert result == expected


def test_data_name_gfdl_hr():
    for data_type in ['ts', 'inst']:
        expected = 'atmos.2006010100-2011123123.temp.nc'
        result = io.data_name_gfdl('temp', 'atmos', data_type,
                                   '3hr', 2010, None, 2000, 6)
        assert result == expected

    with pytest.raises(NameError):
        io.data_name_gfdl('temp', 'atmos', 'av',
                          '3hr', 2010, None, 2000, 6)

    expected = 'atmos.2006-2011.01-12.nc'
    result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                               '3hr', 2010, None, 2000, 6)
    assert result == expected


def test_data_name_gfdl_seasonal():
    for data_type in ['ts', 'inst']:
        with pytest.raises(NameError):
            io.data_name_gfdl('temp', 'atmos', data_type,
                              'seasonal', 2010, None, 2000, 6)

    for intvl_type in ['seasonal', 'seas']:
        expected = 'atmos.2010.JJA.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av',
                                   intvl_type, 2010, 'jja', 2000, 1)
        assert result == expected

        expected = 'atmos.2006-2011.JJA.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av',
                                   intvl_type, 2010, 'jja', 2000, 6)
        assert result == expected

    expected = 'atmos.2006-2011.01-12.nc'
    result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                               'seasonal', 2010, None, 2000, 6)
    assert result == expected


@pytest.fixture()
def load_variable_data_loader():
    return NestedDictDataLoader(file_map, upcast_float32=False)


_LOAD_VAR_DATE_RANGES = {
    'datetime': (datetime.datetime(5, 1, 1),
                 datetime.datetime(5, 12, 31)),
    'datetime64': (np.datetime64('0005-01-01'),
                   np.datetime64('0005-12-31')),
    'cftime': (DatetimeNoLeap(5, 1, 1), DatetimeNoLeap(5, 12, 31)),
    'str': ('0005', '0005')
}


@pytest.mark.parametrize(['start_date', 'end_date'],
                         _LOAD_VAR_DATE_RANGES.values(),
                         ids=list(_LOAD_VAR_DATE_RANGES.keys()))
def test_load_variable(load_variable_data_loader, start_date, end_date):
    result = load_variable_data_loader.load_variable(
        condensation_rain, start_date, end_date, intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00050101.precip_monthly.nc')
    expected = _open_ds_catch_warnings(filepath)['condensation_rain']
    np.testing.assert_array_equal(result.values, expected.values)


@pytest.mark.filterwarnings('ignore:CFTimeIndex.data is deprecated')
@pytest.mark.parametrize(['start_date', 'end_date'],
                         _LOAD_VAR_DATE_RANGES.values(),
                         ids=list(_LOAD_VAR_DATE_RANGES.keys()))
def test_load_variable_does_not_warn(load_variable_data_loader,
                                     start_date, end_date):
    with warnings.catch_warnings(record=True) as warnlog:
        load_variable_data_loader.load_variable(
            condensation_rain,
            start_date, end_date,
            intvl_in='monthly')
    assert len(warnlog) == 0


@pytest.mark.parametrize(['start_date', 'end_date'],
                         _LOAD_VAR_DATE_RANGES.values(),
                         ids=list(_LOAD_VAR_DATE_RANGES.keys()))
def test_load_variable_float32_to_float64(load_variable_data_loader,
                                          start_date, end_date):
    def preprocess(ds, **kwargs):
        # This function converts testing data to the float32 datatype
        return ds.astype(np.float32)
    load_variable_data_loader.upcast_float32 = True
    load_variable_data_loader.preprocess_func = preprocess
    result = load_variable_data_loader.load_variable(
        condensation_rain, start_date,
        end_date,
        intvl_in='monthly').dtype
    expected = np.float64
    assert result == expected


@pytest.mark.parametrize(['start_date', 'end_date'],
                         _LOAD_VAR_DATE_RANGES.values(),
                         ids=list(_LOAD_VAR_DATE_RANGES.keys()))
def test_load_variable_maintain_float32(load_variable_data_loader,
                                        start_date, end_date):
    def preprocess(ds, **kwargs):
        # This function converts testing data to the float32 datatype
        return ds.astype(np.float32)
    load_variable_data_loader.preprocess_func = preprocess
    load_variable_data_loader.upcast_float32 = False
    result = load_variable_data_loader.load_variable(
        condensation_rain, start_date,
        end_date,
        intvl_in='monthly').dtype
    expected = np.float32
    assert result == expected


_LOAD_VAR_MULTI_YEAR_RANGES = {
    'datetime': (datetime.datetime(4, 1, 1),
                 datetime.datetime(5, 12, 31)),
    'datetime64': (np.datetime64('0004-01-01'),
                   np.datetime64('0005-12-31')),
    'cftime': (DatetimeNoLeap(4, 1, 1), DatetimeNoLeap(5, 12, 31)),
    'str': ('0004', '0005')
}


@pytest.mark.parametrize(['start_date', 'end_date'],
                         _LOAD_VAR_MULTI_YEAR_RANGES.values(),
                         ids=list(_LOAD_VAR_MULTI_YEAR_RANGES.keys()))
def test_load_variable_data_vars_all(load_variable_data_loader,
                                     start_date, end_date):
    def preprocess(ds, **kwargs):
        # This function drops the time coordinate from condensation_rain
        temp = ds[condensation_rain.name]
        temp = temp.isel(time=0, drop=True)
        ds = ds.drop(condensation_rain.name)
        ds[condensation_rain.name] = temp
        assert TIME_STR not in ds[condensation_rain.name].coords
        return ds

    load_variable_data_loader.data_vars = 'all'
    load_variable_data_loader.preprocess_func = preprocess
    data = load_variable_data_loader.load_variable(
        condensation_rain, start_date, end_date,
        intvl_in='monthly')
    assert TIME_STR in data.coords


@pytest.mark.parametrize(['start_date', 'end_date'],
                         _LOAD_VAR_MULTI_YEAR_RANGES.values(),
                         ids=list(_LOAD_VAR_MULTI_YEAR_RANGES.keys()))
def test_load_variable_data_vars_default(load_variable_data_loader, start_date,
                                         end_date):
    data = load_variable_data_loader.load_variable(
        condensation_rain, start_date, end_date,
        intvl_in='monthly')
    assert TIME_STR in data.coords


@pytest.mark.parametrize(['start_date', 'end_date'],
                         _LOAD_VAR_MULTI_YEAR_RANGES.values(),
                         ids=list(_LOAD_VAR_MULTI_YEAR_RANGES.keys()))
def test_load_variable_coords_all(load_variable_data_loader, start_date,
                                  end_date):
    load_variable_data_loader.coords = 'all'
    data = load_variable_data_loader.load_variable(
        condensation_rain, start_date, end_date,
        intvl_in='monthly')
    assert TIME_STR in data[ZSURF_STR].coords


@pytest.mark.parametrize('year', [4, 5, 6])
def test_load_variable_non_0001_refdate(load_variable_data_loader, year):
    def preprocess(ds, **kwargs):
        # This function converts our testing data (encoded with a units
        # attribute with a reference data of 0001-01-01) to one
        # with a reference data of 0004-01-01 (to do so we also need
        # to offset the raw time values by three years).
        three_yrs = 1095.
        ds['time'] = ds['time'] - three_yrs
        ds['time'].attrs['units'] = 'days since 0004-01-01 00:00:00'
        ds['time'].attrs['calendar'] = 'noleap'
        ds['time_bounds'] = ds['time_bounds'] - three_yrs
        ds['time_bounds'].attrs['units'] = 'days since 0004-01-01 00:00:00'
        ds['time_bounds'].attrs['calendar'] = 'noleap'
        return ds

    load_variable_data_loader.preprocess_func = preprocess

    result = load_variable_data_loader.load_variable(
        condensation_rain, DatetimeNoLeap(year, 1, 1),
        DatetimeNoLeap(year, 12, 31),
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '000{}0101.precip_monthly.nc'.format(year))
    expected = _open_ds_catch_warnings(filepath)['condensation_rain']
    np.testing.assert_allclose(result.values, expected.values)


def test_load_variable_preprocess(load_variable_data_loader):
    def preprocess(ds, **kwargs):
        if kwargs['start_date'] == DatetimeNoLeap(5, 1, 1):
            ds['condensation_rain'] = 10. * ds['condensation_rain']
        return ds

    load_variable_data_loader.preprocess_func = preprocess

    result = load_variable_data_loader.load_variable(
        condensation_rain, DatetimeNoLeap(5, 1, 1),
        DatetimeNoLeap(5, 12, 31),
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00050101.precip_monthly.nc')
    expected = 10. * _open_ds_catch_warnings(filepath)['condensation_rain']
    np.testing.assert_allclose(result.values, expected.values)

    result = load_variable_data_loader.load_variable(
        condensation_rain, DatetimeNoLeap(4, 1, 1),
        DatetimeNoLeap(4, 12, 31),
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00040101.precip_monthly.nc')
    expected = _open_ds_catch_warnings(filepath)['condensation_rain']
    np.testing.assert_allclose(result.values, expected.values)


def test_load_variable_mask_and_scale(load_variable_data_loader):
    def convert_all_to_missing_val(ds, **kwargs):
        ds['condensation_rain'] = 0. * ds['condensation_rain'] + 1.0e20
        ds['condensation_rain'].attrs['_FillValue'] = 1.0e20
        return ds

    load_variable_data_loader.preprocess_func = convert_all_to_missing_val

    data = load_variable_data_loader.load_variable(
        condensation_rain, DatetimeNoLeap(5, 1, 1),
        DatetimeNoLeap(5, 12, 31),
        intvl_in='monthly')

    num_non_missing = np.isfinite(data).sum().item()
    expected_num_non_missing = 0
    assert num_non_missing == expected_num_non_missing


def test_recursively_compute_variable_native(load_variable_data_loader):
    result = load_variable_data_loader.recursively_compute_variable(
        condensation_rain, DatetimeNoLeap(5, 1, 1),
        DatetimeNoLeap(5, 12, 31),
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00050101.precip_monthly.nc')
    expected = _open_ds_catch_warnings(filepath)['condensation_rain']
    np.testing.assert_array_equal(result.values, expected.values)


def test_recursively_compute_variable_one_level(load_variable_data_loader):
    one_level = Var(
        name='one_level', variables=(condensation_rain, condensation_rain),
        func=lambda x, y: x + y)
    result = load_variable_data_loader.recursively_compute_variable(
        one_level, DatetimeNoLeap(5, 1, 1), DatetimeNoLeap(5, 12, 31),
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00050101.precip_monthly.nc')
    expected = 2. * _open_ds_catch_warnings(filepath)['condensation_rain']
    np.testing.assert_array_equal(result.values, expected.values)


def test_recursively_compute_variable_multi_level(load_variable_data_loader):
    one_level = Var(
        name='one_level', variables=(condensation_rain, condensation_rain),
        func=lambda x, y: x + y)
    multi_level = Var(
        name='multi_level', variables=(one_level, condensation_rain),
        func=lambda x, y: x + y)
    result = load_variable_data_loader.recursively_compute_variable(
        multi_level, DatetimeNoLeap(5, 1, 1), DatetimeNoLeap(5, 12, 31),
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00050101.precip_monthly.nc')
    expected = 3. * _open_ds_catch_warnings(filepath)['condensation_rain']
    np.testing.assert_array_equal(result.values, expected.values)


def test_recursively_compute_grid_attr(load_variable_data_loader):
    result = load_variable_data_loader.recursively_compute_variable(
        bk, DatetimeNoLeap(5, 1, 1),
        DatetimeNoLeap(5, 12, 31), model=example_model,
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00060101.sphum_monthly.nc')
    expected = _open_ds_catch_warnings(filepath)['bk']
    np.testing.assert_array_equal(result.values, expected.values)


def test_recursively_compute_grid_attr_multi_level(load_variable_data_loader):
    one_level = Var(
        name='one_level', variables=(bk, ),
        func=lambda x: 2 * x)
    multi_level = Var(
        name='multi_level', variables=(one_level, bk),
        func=lambda x, y: x + y)
    result = load_variable_data_loader.recursively_compute_variable(
        multi_level, DatetimeNoLeap(5, 1, 1),
        DatetimeNoLeap(5, 12, 31), model=example_model,
        intvl_in='monthly')
    filepath = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '00060101.sphum_monthly.nc')
    expected = 3 * _open_ds_catch_warnings(filepath)['bk']
    np.testing.assert_array_equal(result.values, expected.values)


def test_recursively_compute_grid_attr_error(load_variable_data_loader):
    # Should fail because zsurf is not provided to the example_model object
    zsurf = Var(name=ZSURF_STR, def_time=False, def_vert=False,
                def_lon=True, def_lat=True)
    with pytest.raises(AttributeError):
        load_variable_data_loader.recursively_compute_variable(
            zsurf, DatetimeNoLeap(5, 1, 1),
            DatetimeNoLeap(5, 12, 31), model=example_model,
            intvl_in='monthly')


if __name__ == '__main__':
    unittest.main()
