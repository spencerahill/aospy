#!/usr/bin/env python
import numpy as np
import pytest
import xarray as xr

from aospy.utils.longitude import Longitude, _maybe_cast_to_lon


_good_init_vals_attrs_objs = {
    -10: [10, 'W', Longitude('10W')],
    190.2: [169.8, 'W', Longitude('169.8W')],
    25: [25, 'E', Longitude('25E')],
    365: [5, 'E', Longitude('5E')],
    '45.5e': [45.5, 'E', Longitude('45.5E')],
    '22.2w': [22.2, 'W', Longitude('22.2W')],
    '0': [0, 'E', Longitude('0E')],
    }


_bad_init_vals = ['10ee', '-20e', '190w', None, 'abc', {'a': 1}]


@pytest.mark.parametrize(('val', 'attrs_and_obj'),
                         zip(_good_init_vals_attrs_objs.keys(),
                             _good_init_vals_attrs_objs.values()))
def test_longitude_init_good(val, attrs_and_obj):
    obj = Longitude(val)
    expected_lon = attrs_and_obj[0]
    expected_hem = attrs_and_obj[1]
    expected_obj = attrs_and_obj[2]
    assert obj.longitude == expected_lon
    assert obj.hemisphere == expected_hem
    assert Longitude(val) == expected_obj


def test_longitude_properties():
    lon = Longitude(5)
    with pytest.raises(ValueError):
        lon.longitude = 10
        lon.hemisphere = 'W'


@pytest.mark.parametrize(
    ('obj', 'expected_val'),
    [(Longitude('10w'), "Longitude('10.0W')"),
     (Longitude(0), "Longitude('0.0E')"),
     (Longitude(180), "Longitude('180.0W')")])
def test_longitude_repr(obj, expected_val):
    assert obj.__repr__() == expected_val


@pytest.mark.parametrize('bad_val', _bad_init_vals)
def test_longitude_init_bad(bad_val):
    with pytest.raises(ValueError):
        Longitude(bad_val)


@pytest.mark.parametrize('val', _good_init_vals_attrs_objs.keys())
def test_maybe_cast_to_lon_good(val):
    assert isinstance(_maybe_cast_to_lon(val), Longitude)


@pytest.mark.parametrize('bad_val', _bad_init_vals)
def test_maybe_cast_to_lon_bad(bad_val):
    assert isinstance(_maybe_cast_to_lon(bad_val), type(bad_val))
    with pytest.raises((ValueError, TypeError)):
        _maybe_cast_to_lon(bad_val, strict=True)


@pytest.mark.parametrize(
    ('obj1', 'obj2'),
    [(Longitude('100W'), Longitude('100W')),
     (Longitude('90E'), Longitude('90E')),
     (Longitude(0), Longitude(0)),
     (Longitude('0E'), 0),
     (Longitude('0E'), 720),
     (Longitude('0E'), [0, 720]),
     (Longitude('0E'), xr.DataArray([0, 720]))])
def test_lon_eq(obj1, obj2):
    assert np.all(obj1 == obj2)
    assert np.all(obj2 == obj1)


@pytest.mark.parametrize(
    ('obj1', 'obj2'),
    [(Longitude('100W'), Longitude('90W')),
     (Longitude('90E'), Longitude('100E')),
     (Longitude('10W'), Longitude('0E')),
     (Longitude('0E'), 10),
     (Longitude('0E'), [5, 10]),
     (Longitude('0E'), xr.DataArray([5, 10]))])
def test_lon_lt(obj1, obj2):
    assert np.all(obj1 < obj2)
    assert np.all(obj2 > obj1)


@pytest.mark.parametrize(
    ('obj1', 'obj2'),
    [(Longitude('90W'), Longitude('100W')),
     (Longitude('100E'), Longitude('90E')),
     (Longitude('0E'), Longitude('10W')),
     (Longitude('0E'), -10),
     (Longitude('0E'), [-10, -5]),
     (Longitude('0E'), xr.DataArray([-10, -5]))])
def test_lon_gt(obj1, obj2):
    assert np.all(obj1 > obj2)
    assert np.all(obj2 < obj1)


@pytest.mark.parametrize(
    ('obj1', 'obj2'),
    [(Longitude('100W'), Longitude('100W')),
     (Longitude('90E'), Longitude('90E')),
     (Longitude(0), Longitude(0)),
     (Longitude('100W'), Longitude('90W')),
     (Longitude('90E'), Longitude('100E')),
     (Longitude('10W'), Longitude('0E')),
     (Longitude('0E'), 10),
     (Longitude('0E'), [10, 0]),
     (Longitude('0E'), xr.DataArray([10, 0]))])
def test_lon_leq(obj1, obj2):
    assert np.all(obj1 <= obj2)
    assert np.all(obj2 >= obj2)


@pytest.mark.parametrize(
    ('obj1', 'obj2'),
    [(Longitude('100W'), Longitude('100W')),
     (Longitude('90E'), Longitude('90E')),
     (Longitude(0), Longitude(0)),
     (Longitude('90W'), Longitude('100W')),
     (Longitude('100E'), Longitude('90E')),
     (Longitude('0E'), Longitude('10W')),
     (Longitude('0E'), -10),
     (Longitude('0E'), [0, -10]),
     (Longitude('0E'), xr.DataArray([0, -10]))])
def test_lon_geq(obj1, obj2):
    assert np.all(obj1 >= obj2)
    assert np.all(obj2 <= obj1)


@pytest.mark.parametrize(
    ('obj', 'expected_val'),
    [(Longitude('100W'), 260),
     (Longitude(0), 0),
     (Longitude('20E'), 20)])
def test_to_0360(obj, expected_val):
    assert obj.to_0360() == expected_val


@pytest.mark.parametrize(
    ('obj', 'expected_val'),
    [(Longitude('100W'), -100),
     (Longitude(0), 0),
     (Longitude('20E'), 20)])
def test_to_pm180(obj, expected_val):
    assert obj.to_pm180() == expected_val


@pytest.mark.parametrize(
    ('obj1', 'obj2', 'expected_val'),
    [(Longitude(1), Longitude(1), Longitude(2)),
     (Longitude(175), Longitude(10), Longitude('175W'))])
def test_lon_add(obj1, obj2, expected_val):
    assert obj1 + obj2 == expected_val


@pytest.mark.parametrize(
    ('obj1', 'obj2', 'expected_val'),
    [(Longitude(1), Longitude(1), Longitude(0)),
     (Longitude(185), Longitude(10), Longitude('175E')),
     (Longitude(370), Longitude(20), Longitude('10W'))])
def test_lon_sub(obj1, obj2, expected_val):
    assert obj1 - obj2 == expected_val


if __name__ == '__main__':
    pass
