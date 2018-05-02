#!/usr/bin/env python
"""Functionality relating to parsing and comparing longitudes."""

import numpy as np
import xarray as xr


def lon_to_0360(lon):
    """Convert longitude(s) to be within [0, 360).

    The Eastern hemisphere corresponds to 0 <= lon + (n*360) < 180, and the
    Western Hemisphere corresponds to 180 <= lon + (n*360) < 360, where 'n' is
    any integer (positive, negative, or zero).

    Parameters
    ----------
    lon : scalar or sequence of scalars
        One or more longitude values to be converted to lie in the [0, 360)
        range

    Returns
    -------
        If ``lon`` is a scalar, then a scalar of the same type in the range [0,
        360).  If ``lon`` is array-like, then an array-like of the same type
        with each element a scalar in the range [0, 360).

    """
    quotient = lon // 360
    return lon - quotient*360


def _lon_in_west_hem(lon):
    if lon_to_0360(lon) >= 180:
        return True
    else:
        return False


def lon_to_pm180(lon):
    """Convert longitude(s) to be within [-180, 180).

    The Eastern hemisphere corresponds to 0 <= lon + (n*360) < 180, and the
    Western Hemisphere corresponds to 180 <= lon + (n*360) < 360, where 'n' is
    any integer (positive, negative, or zero).

    Parameters
    ----------
    lon : scalar or sequence of scalars
        One or more longitude values to be converted to lie in the [-180, 180)
        range

    Returns
    -------
        If ``lon`` is a scalar, then a scalar of the same type in the range
        [-180, 180).  If ``lon`` is array-like, then an array-like of the same
        type with each element a scalar in the range [-180, 180).

    """
    lon0360 = lon_to_0360(lon)
    if _lon_in_west_hem(lon0360):
        return lon0360 - 360
    else:
        return lon0360


def _maybe_cast_to_lon(obj, strict=False):
    if isinstance(obj, Longitude):
        return obj
    try:
        return Longitude(obj)
    except (ValueError, TypeError) as e:
        if strict:
            raise type(e)(str(e))
        else:
            return obj


def _other_to_lon(func):
    """Wrapper for casting Longitude operator arguments to Longitude"""
    def func_other_to_lon(obj, other):
        return func(obj, _maybe_cast_to_lon(other))
    return func_other_to_lon


class Longitude(object):
    """Geographic longitude.

    Enables unambiguous comparison of longitudes using the standard comparison
    operators, regardless of they were initially represented with a 0 to 360
    convention, -180 to 180 convention, or anything else, and even if the
    original convention differs between them.

    Specifically, the ``<`` operator assesses if the first object is to the
    west of the second object, with the standard convention that longitudes in
    the Western Hemisphere are always to the west of longitudes in the Eastern
    Hemisphere.  The ``>`` operator is defined analogously.  ``==``, ``>=``,
    and ``<=`` are also all defined.

    In addition to other Longitude objects, the operators can be used to
    compare a Longitude object to anything that can be casted to a Longitude
    object, or to any sequence (e.g. a list or xarray.DataArray) whose elements
    can be casted to Longitude objects.

    """
    def __init__(self, value):
        """
        Parameters
        ----------
        value : {scalar, str}
            Scalars get converted to longitudes using the convention that 0-180
            corresponds to the Eastern Hemisphere, 180-360 corresponds to the
            Western Hemisphere, 360-540 the Eastern Hemisphere, and so on,
            including for negative numbers.

            Strings must be castable to a float or be a positive number in the
            range 0-180 followed by a single letter 'e' or 'w' (case
            insensitive).  For example, ``Longitude('10w')`` would yield a
            ``Longitude`` object corresponding to 10 degrees west longitude.
        """
        try:
            val_as_float = float(value)
        except (ValueError, TypeError):
            if not isinstance(value, str):
                raise ValueError('value must be a scalar or a string')
            if value[-1].lower() not in ('w', 'e'):
                raise ValueError("string inputs must end in 'e' or 'w'")
            try:
                lon_value = float(value[:-1])
            except ValueError:
                raise ValueError('improperly formatted string')
            if (lon_value < 0) or (lon_value > 180):
                raise ValueError('Value given as strings with hemisphere '
                                 'identifier must have numerical values '
                                 'within 0 and +180.  Value given: '
                                 '{}'.format(lon_value))
            self._longitude = lon_value
            self._hemisphere = value[-1].upper()
        else:
            lon_pm180 = lon_to_pm180(val_as_float)
            if _lon_in_west_hem(val_as_float):
                self._longitude = abs(lon_pm180)
                self._hemisphere = 'W'
            else:
                self._longitude = lon_pm180
                self._hemisphere = 'E'

    @property
    def longitude(self):
        """(scalar) The unsigned numerical value of the longitude.

        Always in the range 0 to 180.  Must be combined with the ``hemisphere``
        attribute to specify the exact latitude.

        """
        return self._longitude

    @longitude.setter
    def longitude(self, value):
        raise ValueError("'longitude' property cannot be modified after "
                         "Longitude object has been created.")

    @property
    def hemisphere(self):
        """{'W', 'E'} The longitude's hemisphere, either western or eastern."""
        return self._hemisphere

    @hemisphere.setter
    def hemisphere(self, value):
        raise ValueError("'hemisphere' property cannot be modified after "
                         "Longitude object has been created.")

    def __repr__(self):
        return "Longitude('{0}{1}')".format(self.longitude, self.hemisphere)

    @_other_to_lon
    def __eq__(self, other):
        if isinstance(other, Longitude):
            return (self.hemisphere == other.hemisphere and
                    self.longitude == other.longitude)
        else:
            return xr.apply_ufunc(np.equal, other, self)

    @_other_to_lon
    def __lt__(self, other):
        if isinstance(other, Longitude):
            if self.hemisphere == 'W':
                if other.hemisphere == 'E':
                    return True
                else:
                    return self.longitude > other.longitude
            else:
                if other.hemisphere == 'W':
                    return False
                else:
                    return self.longitude < other.longitude
        else:
            return xr.apply_ufunc(np.greater, other, self)

    @_other_to_lon
    def __gt__(self, other):
        if isinstance(other, Longitude):
            if self.hemisphere == 'W':
                if other.hemisphere == 'E':
                    return False
                else:
                    return self.longitude < other.longitude
            else:
                if other.hemisphere == 'W':
                    return True
                else:
                    return self.longitude > other.longitude
        else:
            return xr.apply_ufunc(np.less, other, self)

    @_other_to_lon
    def __le__(self, other):
        if isinstance(other, Longitude):
            return self < other or self == other
        else:
            return xr.apply_ufunc(np.greater_equal, other, self)

    @_other_to_lon
    def __ge__(self, other):
        if isinstance(other, Longitude):
            return self > other or self == other
        else:
            return xr.apply_ufunc(np.less_equal, other, self)

    def to_0360(self):
        """Convert longitude to its numerical value within [0, 360)."""
        if self.hemisphere == 'W':
            return -1*self.longitude + 360
        else:
            return self.longitude

    def to_pm180(self):
        """Convert longitude to its numerical value within [-180, 180)."""
        if self.hemisphere == 'W':
            return -1*self.longitude
        else:
            return self.longitude

    @_other_to_lon
    def __add__(self, other):
        return Longitude(self.to_0360() + other.to_0360())

    @_other_to_lon
    def __sub__(self, other):
        return Longitude(self.to_0360() - other.to_0360())


if __name__ == '__main__':
    pass
