import numpy as np
import pytest
import xarray as xr

from aospy import Region
from aospy.region import (
    _get_land_mask,
    BoundsRect,
)
from aospy.internal_names import (
    LAT_STR,
    LON_STR,
    SFC_AREA_STR,
    LAND_MASK_STR
)
from aospy.utils import Longitude


@pytest.fixture()
def values_for_reg_arr():
    return np.array([[-2., 1.],
                     [np.nan, 5.],
                     [3., 3.],
                     [4., 4.2]])


@pytest.fixture()
def data_for_reg_calcs(values_for_reg_arr):
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

    da = xr.DataArray(values_for_reg_arr, coords=[lat, lon])
    da.coords[SFC_AREA_STR] = sfc_area
    da.coords[LAND_MASK_STR] = land_mask
    return da


_alt_names = {LON_STR: 'LONS', LAT_STR: 'LATS', LAND_MASK_STR: 'lm',
              SFC_AREA_STR: 'AREA'}


@pytest.fixture()
def data_reg_alt_names(data_for_reg_calcs):
    return data_for_reg_calcs.rename(_alt_names)


region_no_land_mask = Region(
    name='test',
    description='Test region with no land mask',
    west_bound=0.,
    east_bound=5,
    south_bound=0,
    north_bound=90.,
    do_land_mask=False
)


region_land_mask = Region(
    name='test',
    description='Test region with land mask',
    west_bound=0.,
    east_bound=5,
    south_bound=0,
    north_bound=90.,
    do_land_mask=True
)


_expected_mask = [[False, False],
                  [True, False],
                  [True, False],
                  [True, False]]


def test_get_land_mask_without_land_mask(data_for_reg_calcs):
    result = _get_land_mask(data_for_reg_calcs,
                            region_no_land_mask.do_land_mask)
    expected = 1
    assert result == expected


def test_get_land_mask_with_land_mask(data_for_reg_calcs):
    result = _get_land_mask(data_for_reg_calcs, region_land_mask.do_land_mask)
    expected = data_for_reg_calcs[LAND_MASK_STR]
    xr.testing.assert_identical(result, expected)


def test_get_land_mask_non_aospy_name(data_reg_alt_names):
    result = _get_land_mask(data_reg_alt_names, region_land_mask.do_land_mask,
                            land_mask_str=_alt_names[LAND_MASK_STR])
    expected = data_reg_alt_names[_alt_names[LAND_MASK_STR]]
    xr.testing.assert_identical(result, expected)


def test_region_init():
    region = Region(
        name='test',
        description='region description',
        west_bound=0.,
        east_bound=5,
        south_bound=0,
        north_bound=90.,
        do_land_mask=True
    )
    assert region.name == 'test'
    assert region.description == 'region description'
    assert isinstance(region.mask_bounds, tuple)
    assert len(region.mask_bounds) == 1
    assert isinstance(region.mask_bounds[0], BoundsRect)
    assert np.all(region.mask_bounds[0] ==
                  (Longitude(0.), Longitude(5), 0, 90.))
    assert region.do_land_mask is True


def test_region_init_mult_rect():
    bounds_in = [[1, 2, 3, 4], (-12, -30, 2.3, 9)]
    region = Region(name='test', mask_bounds=bounds_in)
    assert isinstance(region.mask_bounds, tuple)
    assert len(region.mask_bounds) == 2
    for (w, e, s, n), bounds in zip(bounds_in, region.mask_bounds):
        assert isinstance(bounds, tuple)
        assert np.all(bounds == (Longitude(w), Longitude(e), s, n))


def test_region_init_bad_bounds():
    with pytest.raises(ValueError):
        Region(mask_bounds=[(1, 2, 3)])
        Region(mask_bounds=[(1, 2, 3, 4),
                            (1, 2, 3)])


def test_make_mask_single_rect(data_for_reg_calcs):
    result = region_land_mask._make_mask(data_for_reg_calcs)
    expected = xr.DataArray(_expected_mask, dims=[LAT_STR, LON_STR],
                            coords={LAT_STR: data_for_reg_calcs[LAT_STR],
                                    LON_STR: data_for_reg_calcs[LON_STR]})
    xr.testing.assert_equal(result.transpose(), expected)


def test_make_mask_mult_rect(data_for_reg_calcs):
    mask_bounds = (region_land_mask.mask_bounds[0], [0, 360, -20, -5])
    region = Region(name='mult_rect', mask_bounds=mask_bounds)
    result = region._make_mask(data_for_reg_calcs)
    expected_values = [[True, True],
                       [True, False],
                       [True, False],
                       [True, False]]
    expected = xr.DataArray(expected_values, dims=[LAT_STR, LON_STR],
                            coords={LAT_STR: data_for_reg_calcs[LAT_STR],
                                    LON_STR: data_for_reg_calcs[LON_STR]})
    xr.testing.assert_equal(result.transpose(), expected)


@pytest.mark.parametrize(
    'region',
    [region_no_land_mask, region_land_mask],
    ids=['no-land-mask', 'land-mask'])
def test_mask_var(data_for_reg_calcs, region):
    # Test region masks first row and second column of test data.  Note that
    # first element of second row is np.nan in initial dataset.
    expected_data = [[np.nan, np.nan],
                     [np.nan, np.nan],
                     [3., np.nan],
                     [4., np.nan]]
    expected = data_for_reg_calcs.copy(deep=True)
    expected.values = expected_data
    result = region.mask_var(data_for_reg_calcs)
    xr.testing.assert_identical(result, expected)


@pytest.mark.parametrize(
    'region',
    [region_no_land_mask, region_land_mask],
    ids=['no-land-mask', 'land-mask'])
def test_mask_var_non_aospy_names(data_reg_alt_names, region):
    # Test region masks first row and second column of test data.  Note that
    # first element of second row is np.nan in initial dataset.
    expected_data = [[np.nan, np.nan],
                     [np.nan, np.nan],
                     [3., np.nan],
                     [4., np.nan]]
    expected = data_reg_alt_names.copy(deep=True)
    expected.values = expected_data
    result = region.mask_var(data_reg_alt_names, lon_str=_alt_names[LON_STR],
                             lat_str=_alt_names[LAT_STR])
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


_map_to_alt_names = {'lon_str': _alt_names[LON_STR],
                     'lat_str': _alt_names[LAT_STR],
                     'land_mask_str': _alt_names[LAND_MASK_STR],
                     'sfc_area_str': _alt_names[SFC_AREA_STR]}


def test_ts_non_aospy_names(data_reg_alt_names):
    result = region_land_mask.ts(data_reg_alt_names, **_map_to_alt_names)
    expected = xr.DataArray(data_reg_alt_names.values[3, 0])
    xr.testing.assert_identical(result, expected)
