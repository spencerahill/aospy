from . import user_path

def _load_user_data(name):
    """Load user data from aospy_path for given module name.

    File must be located in the `aospy_path` directory and be the same name
    as the desired aospy module subpackage, namely one of `regions`, `calcs`,
    `variables`, and `projects`.
    """
    import imp
    return imp.load_source(
        name, '/'.join([user_path, name, '__init__.py']).replace('//','/')
    )

def _set_attr_from_init_kwargs(obj, kwargs):
    """
    Given a dict of keyword arguments and their values as input for an
    __init__() call, set each attribute of the object with the name and value
    of the inputted arguments.
    """
    for kwarg_key, kwarg_value in kwargs.iteritems():
        setattr(obj, kwarg_key, kwarg_value)

def _get_parent_attr(obj, attr):
    """
    Check if the object has the given attribute and it is non-empty.  If not,
    check each parent object for the attribute and use the first one found.
    """
    # Get the attribute from the object itself, if it exists.
    try:
        return getattr(obj, attr)
    # Otherwise, loop through the parent objects until finding the attr.
    except AttributeError:
        for parent_obj in ('var', 'run', 'model', 'proj'):
            try:
                return getattr(getattr(obj, parent_obj), attr)
            except AttributeError:
                pass

def _set_named_attr_dict(obj, attr, attr_list):
    """Set attribute lists that are dicts of {name: value} type."""
    setattr(obj, attr, {item.name: item for item in attr_list})

def _set_parent_object_of_child(parent_obj, parent_attr, child):
    setattr(child, parent_attr, parent_obj)

def _get_nc_grid_objs(obj):
    """Get the nc_grid_objs of an aospy object."""
    from netCDF4 import Dataset, MFDataset
    nc_grid_paths = _get_parent_attr(obj, 'nc_grid_paths')
    nc_objs = []
    for path in nc_grid_paths:
        try:
            nc_obj = Dataset(path)
        except TypeError:
            nc_obj = MFDataset(path)
        finally:
            nc_objs.append(nc_obj)
    return tuple(nc_objs)

def _set_attr_from_nc_grid(nc_grid_obj, obj, attr, attr_name):
    """Set attribute that comes from an nc_grid file."""
    for nc in nc_grid_obj:
        try:
            val = nc.variables[attr_name][:]
        except KeyError:
            pass
        else:
            setattr(obj, attr, val)
            break
    # If not in any nc_grid file, set to None
    if not hasattr(obj, attr):
        setattr(obj, attr, None)

def _set_mult_nc_grid_attr(obj):
    """
    Set multiple attrs from grid file given their names in the grid file.
    """
    nc_grid_obj = _get_nc_grid_objs(obj)
    grid_attrs = {
        'lat':         ('lat', 'latitude', 'LATITUDE', 'y'),
        'lat_bounds':  ('latb', 'lat_bnds', 'lat_bounds'),
        'lon':         ('lon', 'longitude', 'LONGITUDE', 'x'),
        'lon_bounds':  ('lonb', 'lon_bnds', 'lon_bounds'),
        'level':       ('level', 'lev', 'plev'),
        'time':        ('time', 'TIME'),
        'time_st':     ('average_T1',),
        'time_end':    ('average_T2',),
        'time_dur':    ('average_DT',),
        'time_bounds': ('time_bounds', 'time_bnds'),
        'sfc_area':    ('area', 'sfc_area'),
        'zsurf':       ('zsurf',),
        'land_mask':   ('land_mask',),
        'pk':          ('pk',),
        'bk':          ('bk',),
        'phalf':       ('phalf',),
        'pfull':       ('pfull',)
    }
    for attr, attr_names in grid_attrs.iteritems():
        for name in attr_names:
            _set_attr_from_nc_grid(nc_grid_obj, obj, attr, name)
    # Close file objects.
    for nc in nc_grid_obj:
        nc.close()

def calc_grid_sfc_area(lonb, latb, gridcenter=False):
    """Calculate surface area of each grid cell in a lon-lat grid."""
    import numpy as np
    from aospy.constants import r_e

    def diff_latlon_bnds(array):
        import numpy as np
        if array.ndim == 1:
            return array[1:] - array[:-1]
        else:
            return array[:,1] - array[:,0]

    dlon = diff_latlon_bnds(lonb)
    # Must compute gridcell edges if given values at cell centers.
    if gridcenter:
        # Grid must be evenly spaced for algorithm to work.
        assert np.allclose(dlon[0], dlon)
        dlon = np.append(dlon, lonb[0] - lonb[-1] + 360.)
        dlat = diff_latlon_bnds(latb)
        assert np.allclose(dlat[0], dlat)
        # First and last array values done separately.
        lat = latb
        latb = np.append(lat[0] - 0.5*dlat[0], lat[:-1] + 0.5*dlat)
        latb = np.append(latb, lat[-1]+0.5*dlat[-1])

    dsinlat = diff_latlon_bnds(np.sin(np.deg2rad(latb)))
    X, Y = np.meshgrid(np.abs(dlon), np.abs(dsinlat))

    return X*Y*(r_e**2)*(np.pi/180.)

def calc_levs_thick(levs):
    """
    Calculates the thickness, in Pa, of each pressure level.

    Assumes that the pressure values given are at the center of that model
    level, except for the lowest value (typically 1000 hPa), which is the
    bottom boundary. The uppermost level extends to 0 hPa.

    """
    import numpy as np
    from aospy.calcs import to_pascal
    # Bottom level extends from levs[0] to halfway betwen levs[0]
    # and levs[1].
    levs = to_pascal(levs)
    lev_thick = [0.5*(levs[0] - levs[1])]
    # Middle levels extend from halfway between [k-1], [k] and [k], [k+1].
    for k in range(1, levs.size-1):
        lev_thick.append(0.5*(levs[k-1] - levs[k+1]))
    # Top level extends from halfway between top two levels to 0 hPa.
    lev_thick.append(0.5*(levs[-2] + levs[-1]))
    # Convert to numpy array and from hectopascals (hPa) to Pascals (Pa).
    return np.array(lev_thick)
