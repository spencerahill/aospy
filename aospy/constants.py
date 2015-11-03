"""Classes and objects pertaining to physical constants."""
from __future__ import division
import numpy as np


class Constant(object):
    """Physical constants used in atmospheric and oceanic sciences."""
    def __init__(self, value, units, description=''):
        if isinstance(value, int):
            self.value = np.int_(value)
        elif isinstance(value, float):
            self.value = np.float_(value)
        else:
            raise NotImplementedError
        self.units = units
        self.description = description

    def __add__(self, other):
        if isinstance(other, Constant):
            if self.units != other.units:
                raise TypeError(
                    "Can't add constants (%s, %s) with different units: "
                    "(%s, %s)" % (self, other, self.units, other.units)
                )
            else:
                other_value = other.value
        else:
            other_value = other
        return np.ma.add(self.value, other_value)

    def __sub__(self, other):
        if isinstance(other, Constant):
            if self.units != other.units:
                raise TypeError(
                    "Can't subtract constants (%s, %s) with different units: "
                    "(%s, %s)" % (self, other, self.units, other.units)
                )
            else:
                other_value = other.value
        else:
            other_value = other
        return np.ma.subtract(self.value, other_value)

    def __mul__(self, other):
        if isinstance(other, Constant):
            other_value = other.value
        else:
            other_value = other
        return np.ma.multiply(self.value, other_value)

    def __truediv__(self, other):
        if isinstance(other, Constant):
            other_value = other.value
        else:
            other_value = other
        return np.ma.divide(self.value, other_value)

    def __floordiv__(self, other):
        if isinstance(other, Constant):
            other_value = other.value
        else:
            other_value = other
        return np.ma.divide(self.value, other_value)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return np.ma.subtract(other, self.value)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rtruediv__(self, other):
        return np.ma.divide(other, self.value)

    def __rfloordiv__(self, other):
        return np.ma.divide(other, self.value)

    def __pow__(self, other):
        if isinstance(other, Constant):
            other_value = other.value
        else:
            other_value = other
        return np.ma.power(self.value, other_value)

    def __rpow__(self, other):
        return np.ma.power(other, self.value)

r_e = Constant(
    6370997.,
    'm',
    description='Mean radius of Earth'
)
c_p = Constant(
    1003.5,
    'J/K/kg.',
    description='Specific heat capacity of dry air at constant pressure'
)
c_v = Constant(
    717.,
    'J/K/kg.',
    description='Specific heat capacity of dry air at constant volume'
)
L_v = Constant(
    2.5e6,
    'J/kg',
    description='Latent heat of vaporization of water'
)
L_f = Constant(
    3.34e5,
    'J/kg',
    description='Latent heat of fusion of water'
)
grav = Constant(
    9.81,
    r'm s$^{-2}$',
    description='Acceleration due to gravity'
)
R_d = Constant(
    287.04,
    'J/K/kg',
    description='Dry air gas constant'
)
R_v = Constant(
    461.5,
    'J/K/kg',
    description='Water vapor gas constant'
)
pi = Constant(
    np.pi,
    'dimensionless',
    description='Mathematical pi: 3.14159...'
)
seconds_in_day = Constant(
    24.*3600.,
    'seconds',
    description='Number of seconds in a day'
)
Omega = Constant(
    2.*pi / seconds_in_day.value,
    's^{-1}',
    description='Rotation rate of Earth'
)
epsilon = Constant(
    R_d.value / R_v.value,
    'dimensionless',
    description='Ratio of gas constants of dry air and water vapor'
)
kappa = Constant(
    R_d.value / c_p.value,
    'dimensionless',
    description=('Ratio of gas constant and specific heat at constant '
                 'pressure of dry air')
)

c_va = Constant(
    719.,
    'J/kg/K',
    description='Specific heat capacity of dry air at constant volume'
)
c_vv = Constant(
    1418.,
    'J/kg/K',
    description='Specific heat capacity of water vapor at constant volume'
)
c_vl = Constant(
    4216.,
    'J/kg/K',
    description='Specific heat capacity of liquid water at constant volume'
)
c_vs = Constant(
    2106.,
    'J/kg/K',
    description='Specific heat capacity of solid water at constant volume'
)
R_a = Constant(
    287.04,
    'J/kg/K',
    description='Dry air gas constant'
)
R_v = Constant(
    461.4,
    'J/kg/K',
    description='Water vapor gas constant'
)
p_trip = Constant(
    611.65,
    'hPa',
    description='Pressure of water triple point'
)
T_trip = Constant(
    273.16,
    'K',
    description='Temperature of water triple point'
)
E_0v = Constant(
    2.374e6,
    'J/kg',
    description=('Difference in specific internal energy between water vapor '
                 'and liquid water at the triple-point temperature')
)
E_0s = Constant(
    3.337e5,
    'J/kg',
    description=('Difference in specific internal energy between liquid water '
                 'and solid water at the triple-point temperature')
)
s_0v = Constant(
    E_0v.value / T_trip.value + R_v.value,
    'J/kg/K',
    description='Specific entropy of water vapor at the triple point.'
)
s_0s = Constant(
    E_0s.value / T_trip.value,
    'J/kg/K',
    description='Specific entropy of solid water at the triple point.'
)
