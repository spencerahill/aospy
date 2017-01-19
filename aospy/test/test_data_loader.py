#!/usr/bin/env python
"""Test suite for aospy.data_loader module."""
from datetime import datetime
import os
import unittest

import numpy as np
import xarray as xr

from aospy.data_loader import (DataLoader, DictDataLoader, GFDLDataLoader,
                               NestedDictDataLoader, rename_grid_attrs,
                               set_grid_attrs_as_coords, _sel_var)
from aospy.internal_names import (LAT_STR, LON_STR, TIME_STR, TIME_BOUNDS_STR,
                                  NV_STR, SFC_AREA_STR, ETA_STR)
from aospy.utils import io
from .data.objects.examples import condensation_rain, convection_rain, precip


class AospyDataLoaderTestCase(unittest.TestCase):
    def setUp(self):
        self.DataLoader = DataLoader()
        self.generate_file_set_args = dict(
            var=condensation_rain, start_date=datetime(2000, 1, 1),
            end_date=datetime(2002, 12, 31), domain='atmos',
            intvl_in='monthly', dtype_in_vert='sigma', dtype_in_time='ts',
            intvl_out=None)
        time_bounds = np.array([[0, 31], [31, 59], [59, 90]])
        nv = np.array([0, 1])
        time = np.array([15, 46, 74])
        data = np.zeros((3, 1, 1))
        lat = [0]
        lon = [0]
        self.ALT_LAT_STR = 'LATITUDE'
        self.var_name = 'a'
        ds = xr.DataArray(data,
                          coords=[time, lat, lon],
                          dims=[TIME_STR, self.ALT_LAT_STR, LON_STR],
                          name=self.var_name).to_dataset()
        ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                           coords=[time, nv],
                                           dims=[TIME_STR, NV_STR],
                                           name=TIME_BOUNDS_STR)
        units_str = 'days since 2000-01-01 00:00:00'
        ds[TIME_STR].attrs['units'] = units_str
        ds[TIME_BOUNDS_STR].attrs['units'] = units_str
        self.ds = ds

        inst_time = np.array([3, 6, 9])
        inst_units_str = 'hours since 2000-01-01 00:00:00'
        inst_ds = ds.copy()
        inst_ds.drop(TIME_BOUNDS_STR)
        inst_ds[TIME_STR].values = inst_time
        inst_ds[TIME_STR].attrs['units'] = inst_units_str
        self.inst_ds = inst_ds

    def tearDown(self):
        pass


class TestDataLoader(AospyDataLoaderTestCase):
    def test_rename_grid_attrs_ds(self):
        assert LAT_STR not in self.ds
        assert self.ALT_LAT_STR in self.ds
        ds = rename_grid_attrs(self.ds)
        assert LAT_STR in ds

    def test_set_grid_attrs_as_coords(self):
        ds = rename_grid_attrs(self.ds)
        sfc_area = ds[self.var_name].isel(**{TIME_STR: 0}).drop(TIME_STR)
        ds[SFC_AREA_STR] = sfc_area

        assert SFC_AREA_STR not in ds[self.var_name]

        ds = set_grid_attrs_as_coords(ds)
        assert SFC_AREA_STR in ds[self.var_name]

    def test_sel_var(self):
        time = np.array([0, 31, 59]) + 15
        data = np.zeros((3))
        ds = xr.DataArray(data,
                          coords=[time],
                          dims=[TIME_STR],
                          name=convection_rain.name).to_dataset()
        condensation_rain_alt_name, = condensation_rain.alt_names
        ds[condensation_rain_alt_name] = xr.DataArray(data, coords=[ds.time])
        result = _sel_var(ds, convection_rain)
        self.assertEqual(result.name, convection_rain.name)

        result = _sel_var(ds, condensation_rain)
        self.assertEqual(result.name, condensation_rain.name)

        with self.assertRaises(KeyError):
            _sel_var(ds, precip)

    def test_maybe_apply_time_shift(self):
        ds = xr.decode_cf(self.ds)
        da = ds[self.var_name]

        result = self.DataLoader._maybe_apply_time_shift(
            da.copy(), **self.generate_file_set_args)[TIME_STR]
        assert result.identical(da[TIME_STR])

        offset = self.DataLoader._maybe_apply_time_shift(
            da.copy(), {'days': 1}, **self.generate_file_set_args)
        result = offset[TIME_STR]

        expected = da[TIME_STR] + np.timedelta64(1, 'D')
        expected[TIME_STR] = expected

        assert result.identical(expected)

    def test_generate_file_set(self):
        with self.assertRaises(NotImplementedError):
            self.DataLoader._generate_file_set()


class TestDictDataLoader(TestDataLoader):
    def setUp(self):
        super(TestDictDataLoader, self).setUp()
        file_map = {'monthly': ['a.nc']}
        self.DataLoader = DictDataLoader(file_map)

    def test_generate_file_set(self):
        result = self.DataLoader._generate_file_set(
            **self.generate_file_set_args)
        expected = ['a.nc']
        self.assertEqual(result, expected)

        with self.assertRaises(KeyError):
            self.generate_file_set_args['intvl_in'] = 'daily'
            result = self.DataLoader._generate_file_set(
                **self.generate_file_set_args)


class TestNestedDictDataLoader(TestDataLoader):
    def setUp(self):
        super(TestNestedDictDataLoader, self).setUp()
        file_map = {'monthly': {'condensation_rain': ['a.nc']}}
        self.DataLoader = NestedDictDataLoader(file_map)

    def test_generate_file_set(self):
        result = self.DataLoader._generate_file_set(
            **self.generate_file_set_args)
        expected = ['a.nc']
        self.assertEqual(result, expected)

        with self.assertRaises(KeyError):
            self.generate_file_set_args['var'] = convection_rain
            result = self.DataLoader._generate_file_set(
                **self.generate_file_set_args)


class TestGFDLDataLoader(TestDataLoader):
    def setUp(self):
        super(TestGFDLDataLoader, self).setUp()
        self.DataLoader = GFDLDataLoader(
            data_direc=os.path.join('.', 'test'),
            data_dur=6,
            data_start_date=datetime(2000, 1, 1),
            data_end_date=datetime(2012, 12, 31)
        )

    def test_overriding_constructor(self):
        new = GFDLDataLoader(self.DataLoader,
                             data_direc=os.path.join('.', 'a'))
        self.assertEqual(new.data_direc, os.path.join('.', 'a'))
        self.assertEqual(new.data_dur, self.DataLoader.data_dur)
        self.assertEqual(new.data_start_date, self.DataLoader.data_start_date)
        self.assertEqual(new.data_end_date, self.DataLoader.data_end_date)

        new = GFDLDataLoader(self.DataLoader, data_dur=8)
        self.assertEqual(new.data_dur, 8)

        new = GFDLDataLoader(self.DataLoader,
                             data_start_date=datetime(2001, 1, 1))
        self.assertEqual(new.data_start_date, datetime(2001, 1, 1))

        new = GFDLDataLoader(self.DataLoader,
                             data_end_date=datetime(2003, 12, 31))
        self.assertEqual(new.data_end_date, datetime(2003, 12, 31))

    def test_maybe_apply_time_offset_inst(self):
        inst_ds = xr.decode_cf(self.inst_ds)
        self.generate_file_set_args['dtype_in_time'] = 'inst'
        self.generate_file_set_args['intvl_in'] = '3hr'
        da = inst_ds[self.var_name]
        result = self.DataLoader._maybe_apply_time_shift(
            da.copy(), **self.generate_file_set_args)[TIME_STR]

        expected = da[TIME_STR] + np.timedelta64(-3, 'h')
        expected[TIME_STR] = expected
        assert result.identical(expected)

        self.generate_file_set_args['intvl_in'] = 'daily'
        da = inst_ds[self.var_name]
        result = self.DataLoader._maybe_apply_time_shift(
            da.copy(), **self.generate_file_set_args)[TIME_STR]

        expected = da[TIME_STR]
        expected[TIME_STR] = expected
        assert result.identical(expected)

    def test_maybe_apply_time_offset_ts(self):
        ds = xr.decode_cf(self.ds)
        da = ds[self.var_name]

        result = self.DataLoader._maybe_apply_time_shift(
            da.copy(), **self.generate_file_set_args)[TIME_STR]
        assert result.identical(da[TIME_STR])

    def test_generate_file_set(self):
        with self.assertRaises(IOError):
            self.DataLoader._generate_file_set(**self.generate_file_set_args)

    def test_input_data_paths_gfdl(self):
        expected = [os.path.join('.', 'test', 'atmos', 'ts', 'monthly', '6yr',
                                 'atmos.200601-201112.temp.nc')]
        result = self.DataLoader._input_data_paths_gfdl(
            'temp', datetime(2010, 1, 1), datetime(2010, 12, 31), 'atmos',
            'monthly', 'pressure', 'ts', None)
        self.assertEqual(result, expected)

        expected = [os.path.join('.', 'test', 'atmos_daily', 'ts', 'daily',
                                 '6yr',
                                 'atmos_daily.20060101-20111231.temp.nc')]
        result = self.DataLoader._input_data_paths_gfdl(
            'temp', datetime(2010, 1, 1), datetime(2010, 12, 31), 'atmos',
            'daily', 'pressure', 'ts', None)
        self.assertEqual(result, expected)

        expected = [os.path.join('.', 'test', 'atmos_level', 'ts', 'monthly',
                                 '6yr', 'atmos_level.200601-201112.temp.nc')]
        result = self.DataLoader._input_data_paths_gfdl(
            'temp', datetime(2010, 1, 1), datetime(2010, 12, 31), 'atmos',
            'monthly', ETA_STR, 'ts', None)
        self.assertEqual(result, expected)

        expected = [os.path.join('.', 'test', 'atmos', 'ts', 'monthly',
                                 '6yr', 'atmos.200601-201112.ps.nc')]
        result = self.DataLoader._input_data_paths_gfdl(
            'ps', datetime(2010, 1, 1), datetime(2010, 12, 31), 'atmos',
            'monthly', ETA_STR, 'ts', None)
        self.assertEqual(result, expected)

        expected = [os.path.join('.', 'test', 'atmos_inst', 'ts', 'monthly',
                                 '6yr', 'atmos_inst.200601-201112.temp.nc')]
        result = self.DataLoader._input_data_paths_gfdl(
            'temp', datetime(2010, 1, 1), datetime(2010, 12, 31), 'atmos',
            'monthly', 'pressure', 'inst', None)
        self.assertEqual(result, expected)

        expected = [os.path.join('.', 'test', 'atmos', 'av', 'monthly_6yr',
                                 'atmos.2006-2011.jja.nc')]
        result = self.DataLoader._input_data_paths_gfdl(
            'temp', datetime(2010, 1, 1), datetime(2010, 12, 31), 'atmos',
            'monthly', 'pressure', 'av', 'jja')
        self.assertEqual(result, expected)

    def test_data_name_gfdl_annual(self):
        for data_type in ['ts', 'inst']:
            expected = 'atmos.2010.temp.nc'
            result = io.data_name_gfdl('temp', 'atmos', data_type,
                                       'annual', 2010, None, 2000, 1)
            self.assertEqual(result, expected)

            expected = 'atmos.2006-2011.temp.nc'
            result = io.data_name_gfdl('temp', 'atmos', data_type,
                                       'annual', 2010, None, 2000, 6)
            self.assertEqual(result, expected)

        for intvl_type in ['annual', 'ann']:
            expected = 'atmos.2010.ann.nc'
            result = io.data_name_gfdl('temp', 'atmos', 'av',
                                       intvl_type, 2010, None, 2000, 1)
            self.assertEqual(result, expected)

            expected = 'atmos.2006-2011.ann.nc'
            result = io.data_name_gfdl('temp', 'atmos', 'av',
                                       intvl_type, 2010, None, 2000, 6)
            self.assertEqual(result, expected)

        expected = 'atmos.2006-2011.01-12.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                                   'annual', 2010, None, 2000, 6)
        self.assertEqual(result, expected)

    def test_data_name_gfdl_monthly(self):
        for data_type in ['ts', 'inst']:
            expected = 'atmos.200601-201112.temp.nc'
            result = io.data_name_gfdl('temp', 'atmos', data_type,
                                       'monthly', 2010, 'jja', 2000, 6)
            self.assertEqual(result, expected)

        for intvl_type in ['monthly', 'mon']:
            expected = 'atmos.2010.jja.nc'
            result = io.data_name_gfdl('temp', 'atmos', 'av',
                                       intvl_type, 2010, 'jja', 2000, 1)
            self.assertEqual(result, expected)

            expected = 'atmos.2006-2011.jja.nc'
            result = io.data_name_gfdl('temp', 'atmos', 'av',
                                       intvl_type, 2010, 'jja', 2000, 6)
            self.assertEqual(result, expected)

        expected = 'atmos.2006-2011.01-12.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                                   'monthly', 2010, 'jja', 2000, 6)
        self.assertEqual(result, expected)

    def test_data_name_gfdl_daily(self):
        for data_type in ['ts', 'inst']:
            expected = 'atmos.20060101-20111231.temp.nc'
            result = io.data_name_gfdl('temp', 'atmos', data_type,
                                       'daily', 2010, None, 2000, 6)
            self.assertEqual(result, expected)

        with self.assertRaises(NameError):
            io.data_name_gfdl('temp', 'atmos', 'av',
                              'daily', 2010, None, 2000, 6)

        expected = 'atmos.2006-2011.01-12.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                                   'daily', 2010, None, 2000, 6)
        self.assertEqual(result, expected)

    def test_data_name_gfdl_hr(self):
        for data_type in ['ts', 'inst']:
            expected = 'atmos.2006010100-2011123123.temp.nc'
            result = io.data_name_gfdl('temp', 'atmos', data_type,
                                       '3hr', 2010, None, 2000, 6)
            self.assertEqual(result, expected)

        with self.assertRaises(NameError):
            io.data_name_gfdl('temp', 'atmos', 'av',
                              '3hr', 2010, None, 2000, 6)

        expected = 'atmos.2006-2011.01-12.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                                   '3hr', 2010, None, 2000, 6)
        self.assertEqual(result, expected)

    def test_data_name_gfdl_seasonal(self):
        for data_type in ['ts', 'inst']:
            with self.assertRaises(NameError):
                io.data_name_gfdl('temp', 'atmos', data_type,
                                  'seasonal', 2010, None, 2000, 6)

        for intvl_type in ['seasonal', 'seas']:
            expected = 'atmos.2010.JJA.nc'
            result = io.data_name_gfdl('temp', 'atmos', 'av',
                                       intvl_type, 2010, 'jja', 2000, 1)
            self.assertEqual(result, expected)

            expected = 'atmos.2006-2011.JJA.nc'
            result = io.data_name_gfdl('temp', 'atmos', 'av',
                                       intvl_type, 2010, 'jja', 2000, 6)
            self.assertEqual(result, expected)

        expected = 'atmos.2006-2011.01-12.nc'
        result = io.data_name_gfdl('temp', 'atmos', 'av_ts',
                                   'seasonal', 2010, None, 2000, 6)
        self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
