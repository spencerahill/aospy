import os

default_colormap = 'RdBu'
user_path = os.path.join(os.getenv('HOME'), 'aospy_user', 'aospy_user')

LON_STR = 'lon'
LAT_STR = 'lat'
PHALF_STR = 'phalf'
PFULL_STR = 'pfull'
PLEVEL_STR = 'level'
TIME_STR = 'time'
TIME_STR_IDEALIZED = 'year'
YEAR_STR = 'year'

del os
