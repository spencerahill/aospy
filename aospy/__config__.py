import os

user_path = os.path.join(os.getenv('HOME'), 'aospy_user', 'aospy_user')

ETA_STR = 'sigma'
LON_STR = 'lon'
LAT_STR = 'lat'
LON_BOUNDS_STR = 'lon_bounds'
LAT_BOUNDS_STR = 'lat_bounds'
PHALF_STR = 'phalf'
PFULL_STR = 'pfull'
PLEVEL_STR = 'level'
TIME_STR = 'time'
TIME_STR_IDEALIZED = 'year'
TIME_BOUNDS_STR = 'time_bounds'
YEAR_STR = 'year'
BOUNDS_STR = 'bnds'

del os
