"""pytest conftest.py file for sharing fixtures across modules."""
import datetime

from cftime import DatetimeNoLeap
import numpy as np
import pytest
import xarray as xr

from aospy.internal_names import (
    LON_STR,
    TIME_STR,
    TIME_BOUNDS_STR,
    BOUNDS_STR,
)


_DATE_RANGES = {
    'datetime': (datetime.datetime(2000, 1, 1),
                 datetime.datetime(2002, 12, 31)),
    'datetime64': (np.datetime64('2000-01-01'),
                   np.datetime64('2002-12-31')),
    'cftime': (DatetimeNoLeap(2000, 1, 1),
               DatetimeNoLeap(2002, 12, 31)),
    'str': ('2000', '2002')
}


@pytest.fixture()
def alt_lat_str():
    return 'LATITUDE'


@pytest.fixture()
def var_name():
    return 'a'


@pytest.fixture
def ds_time_encoded_cf():
    time_bounds = np.array([[0, 31], [31, 59], [59, 90]])
    bounds = np.array([0, 1])
    time = np.array([15, 46, 74])
    data = np.zeros((3))
    ds = xr.DataArray(data,
                      coords=[time],
                      dims=[TIME_STR],
                      name='a').to_dataset()
    ds[TIME_BOUNDS_STR] = xr.DataArray(time_bounds,
                                       coords=[time, bounds],
                                       dims=[TIME_STR, BOUNDS_STR],
                                       name=TIME_BOUNDS_STR)
    units_str = 'days since 2000-01-01 00:00:00'
    cal_str = 'noleap'
    ds[TIME_STR].attrs['units'] = units_str
    ds[TIME_STR].attrs['calendar'] = cal_str
    return ds


@pytest.fixture()
def ds_with_time_bounds(ds_time_encoded_cf, alt_lat_str, var_name):
    time = ds_time_encoded_cf[TIME_STR]
    data = np.zeros((3, 1, 1))
    lat = [0]
    lon = [0]

    ds = xr.DataArray(
        data,
        coords=[time, lat, lon],
        dims=[TIME_STR, alt_lat_str, LON_STR],
        name=var_name,
    ).to_dataset()
    ds[TIME_BOUNDS_STR] = ds_time_encoded_cf[TIME_BOUNDS_STR]
    return ds


@pytest.fixture()
def ds_inst(ds_with_time_bounds):
    inst_time = np.array([3, 6, 9])
    inst_units_str = 'hours since 2000-01-01 00:00:00'
    ds = ds_with_time_bounds.drop_vars([BOUNDS_STR, TIME_BOUNDS_STR])
    ds[TIME_STR].values = inst_time
    ds[TIME_STR].attrs['units'] = inst_units_str
    ds[TIME_STR].attrs['calendar'] = 'noleap'
    return ds


@pytest.fixture()
def ds_no_time(ds_with_time_bounds):
    return ds_with_time_bounds.drop_vars(TIME_STR)
