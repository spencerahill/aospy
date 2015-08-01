"""Classes and objects pertaining to physical constants."""
import numpy as np


class Constant(object):
    """Physical constants used in atmospheric and oceanic sciences."""
    def __init__(self, value, units, description=''):
        self.value = value
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

    def __div__(self, other):
        if isinstance(other, Constant):
            other_value = other.value
        else:
            other_value = other
        return np.ma.divide(self.value, other_value)

    def __radd__(self, other):
        return self.__add__(other)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __rdiv__(self, other):
        return self.__div__(other)


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
    2.*pi / seconds_in_day,
    's^{-1}',
    description='Rotation rate of Earth'
)
epsilon = Constant(
    R_d / R_v,
    'dimensionless',
    description='Ratio of gas constants of dry air and water vapor'
)
kappa = Constant(
    R_d / c_p,
    'dimensionless',
    description=('Ratio of gas constant and specific heat at constant '
                 'pressure of dry air')
)

c_va = 719
c_vv = 1418
c_vl = 4216
c_vs = 2106
R_a = 287.04
