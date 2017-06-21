import numpy as np
import pytest
import xarray as xr

from aospy import Region
from aospy.region import _get_land_mask
from aospy.internal_names import (LAT_STR, LON_STR,
                                  SFC_AREA_STR, LAND_MASK_STR)


@pytest.fixture()
def data_for_reg_calcs():
    lat = [-10., 1., 10., 20.]
    lon = [1., 10.]
    sfc_area = [0.5, 1., 0.5, 0.25]
    land_mask = [1., 1., 0., 1.]

    lat = xr.DataArray(lat, dims=[LAT_STR], coords=[lat])
    lon = xr.DataArray(lon, dims=[LON_STR], coords=[lon])
    sfc_area = xr.DataArray(sfc_area, dims=[LAT_STR], coords=[lat])
    land_mask = xr.DataArray(land_mask, dims=[LAT_STR], coords=[lat])

    sfc_area, _ = xr.broadcast(sfc_area, lon)
    land_mask, _ = xr.broadcast(land_mask, lon)

    da = xr.DataArray([[2., 2.],
                       [np.nan, 5.],
                       [3., 3.],
                       [4., 4.]], coords=[lat, lon])
    da[SFC_AREA_STR] = sfc_area
    da[LAND_MASK_STR] = land_mask
    return da


region_no_land_mask = Region(
    name='test',
    description='Test region with no land mask',
    lat_bounds=(0., 90.),
    lon_bounds=(0., 5.),
    do_land_mask=False
)


region_land_mask = Region(
    name='test',
    description='Test region with no land mask',
    lat_bounds=(0., 90.),
    lon_bounds=(0., 5.),
    do_land_mask=True
)


@pytest.mark.parametrize(
    'region',
    [region_no_land_mask, region_land_mask],
    ids=['no-land-mask', 'land-mask'])
def test_mask_var(data_for_reg_calcs, region):
    # Test region masks first row and second column
    # of test data.  Note that first element of
    # second row is np.nan in initial dataset.
    expected_data = [[np.nan, np.nan],
                     [np.nan, np.nan],
                     [3., np.nan],
                     [4., np.nan]]
    expected = data_for_reg_calcs.copy(deep=True)
    expected.values = expected_data

    result = region.mask_var(data_for_reg_calcs)
    xr.testing.assert_identical(result, expected)


def test_get_land_mask_without_land_mask(data_for_reg_calcs):
    result = _get_land_mask(data_for_reg_calcs,
                            region_no_land_mask.do_land_mask)
    expected = 1
    assert result == expected


def test_get_land_mask_with_land_mask(data_for_reg_calcs):
    result = _get_land_mask(data_for_reg_calcs, region_land_mask.do_land_mask)
    expected = data_for_reg_calcs[LAND_MASK_STR]
    xr.testing.assert_identical(result, expected)


def test_ts_no_land_mask(data_for_reg_calcs):
    result = region_no_land_mask.ts(data_for_reg_calcs)

    data = data_for_reg_calcs.values
    sfc_area = data_for_reg_calcs.sfc_area.values
    exp_numerator = data[2, 0] * sfc_area[2, 0] + data[3, 0] * sfc_area[3, 0]
    exp_denominator = sfc_area[2, 0] + sfc_area[3, 0]
    expected = xr.DataArray(exp_numerator / exp_denominator)
    xr.testing.assert_identical(result, expected)


def test_ts_land_mask(data_for_reg_calcs):
    result = region_land_mask.ts(data_for_reg_calcs)
    expected = xr.DataArray(data_for_reg_calcs.values[3, 0])
    xr.testing.assert_identical(result, expected)
