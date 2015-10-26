import os as _os

default_colormap = 'RdBu'
user_path = '/'.join([_os.getenv('HOME'),
                      'aospy_user/aospy_user']).replace('//', '/')

LON_STR = 'lon'
LAT_STR = 'lat'
PHALF_STR = 'phalf'
PFULL_STR = 'pfull'
PLEVEL_STR = 'level'
TIME_STR = 'time'
TIME_STR_IDEALIZED = 'year'

del _os
