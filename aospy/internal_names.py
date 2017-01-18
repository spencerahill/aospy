"""Internal names used throughout aospy

Depending on its source, data corresponding to the same physical
variable can have different internal names. To ensure that all
variables can be easily accessed, upon loading in data ``aospy``
overwrites these names with names specified either in ``aospy.Var``
objects (in the case of variables like 'precip'), or names specified
in this module (in the case of frequently encountered coordinate variables,
like 'latitude').

The ``GRID_ATTRS`` dictionary maps the internal names for coordinates used
in ``aospy`` to lists of potential names based on all previously encountered
names for each coordinate in different datasets.  ``aospy`` uses this
dictionary to identify given physical coordinates, and rename them with
their proper internal names.
"""
from collections import OrderedDict

# Horizontal coordinates
LON_STR = 'lon'
LAT_STR = 'lat'
LON_BOUNDS_STR = 'lon_bounds'
LAT_BOUNDS_STR = 'lat_bounds'
SFC_AREA_STR = 'sfc_area'
LAND_MASK_STR = 'land_mask'

# Vertical coordinates
ETA_STR = 'sigma'
PHALF_STR = 'phalf'
PFULL_STR = 'pfull'
PLEVEL_STR = 'level'
PK_STR = 'pk'
BK_STR = 'bk'
ZSURF_STR = 'zsurf'

# Time coordinates
TIME_STR = 'time'
TIME_BOUNDS_STR = 'time_bounds'
YEAR_STR = 'year'
BOUNDS_STR = 'bnds'
AVERAGE_T1_STR = 'average_T1'
AVERAGE_T2_STR = 'average_T2'
AVERAGE_DT_STR = 'average_DT'
NV_STR = 'nv'
AVG_START_DATE_STR = 'avg_start_date'
AVG_END_DATE_STR = 'avg_end_date'

# AVERAGE_DT_STR does not enter here, because it is not a 'time since'
# quantity.
TIME_VAR_STRS = [TIME_STR, TIME_BOUNDS_STR]

GRID_ATTRS = OrderedDict(
    [(LAT_STR, ('lat', 'latitude', 'LATITUDE', 'y', 'yto')),
     (LAT_BOUNDS_STR, ('latb', 'lat_bnds', 'lat_bounds')),
     (LON_STR, ('lon', 'longitude', 'LONGITUDE', 'x', 'xto')),
     (LON_BOUNDS_STR, ('lonb', 'lon_bnds', 'lon_bounds')),
     (ZSURF_STR, ('zsurf',)),
     (SFC_AREA_STR, ('area', 'sfc_area')),
     (LAND_MASK_STR, ('land_mask',)),
     (PK_STR, ('pk',)),
     (BK_STR, ('bk',)),
     (PHALF_STR, ('phalf',)),
     (PFULL_STR, ('pfull',)),
     (PLEVEL_STR, ('level', 'lev', 'plev')),
     (TIME_STR, ('time',)),
     (AVERAGE_DT_STR, ('average_DT',)),
     (TIME_BOUNDS_STR, ('time_bounds', 'time_bnds')),
     (NV_STR, ('nv',)),
     (AVG_START_DATE_STR, ('avg_start_date',)),
     (AVG_END_DATE_STR, ('avg_end_date',))]
)
