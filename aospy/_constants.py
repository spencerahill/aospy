"""Physical constants used in aospy"""
import xarray as xr


RADIUS_EARTH = xr.DataArray(
    6370997.,
    attrs={'units': 'm',
           'description': 'Mean radius of Earth'}
)
GRAV_EARTH = xr.DataArray(
    9.81,
    attrs={'units': r'm s$^{-2}$',
           'description': 'Acceleration due to gravity'}
)
