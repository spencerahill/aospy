"""aospy.utils: utility functions for the aospy module."""
import numpy as np
import pandas as pd
import xray

from . import user_path
from .constants import grav

TIME_STR = 'time'


def coord_to_new_dataarray(arr, dim):
    """Create a DataArray comprising the coord for the specified dim.

    Useful, for example, when wanting to resample in time, because at least
    for xray 0.6.0 and prior, the `resample` method doesn't work when applied
    to coords.  The DataArray returned by this method lacks that limitation.
    """
    return xray.DataArray(arr[dim].values, coords=[arr[dim].values],
                          dims=[dim])


def apply_time_offset(time, months=0, days=0, hours=0):
    """Apply the given offset to the given time array.

    This is useful for GFDL model output of instantaneous values.  For example,
    3 hourly data postprocessed to netCDF files spanning 1 year each will
    actually have time values that are offset by 3 hours, such that the first
    value is for 1 Jan 03:00 and the last value is 1 Jan 00:00 of the
    subsequent year.  This causes problems in xray, e.g. when trying to group
    by month.  It is resolved by manually subtracting off those three hours,
    such that the dates span from 1 Jan 00:00 to 31 Dec 21:00 as desired.
    """
    return (pd.to_datetime(time.values) +
            pd.tseries.offsets.DateOffset(months=months, days=days,
                                          hours=hours))


def monthly_mean_ts(arr):
    """Convert a sub-monthly time-series into one of monthly means."""
    return arr.resample('1M', TIME_STR, how='mean')


def monthly_mean_at_each_ind(arr_mon, arr_sub):
    """Copy monthly mean over each time index in that month."""
    time = arr_mon[TIME_STR]
    start = time.indexes[TIME_STR][0].replace(day=1, hour=0)
    end = time.indexes[TIME_STR][-1]
    new_indices = pd.DatetimeIndex(start=start, end=end, freq='MS')
    arr_new = arr_mon.reindex(time=new_indices, method='backfill')
    return arr_new.reindex_like(arr_sub, method='pad')


def load_user_data(name):
    """Load user data from aospy_path for given module name.

    File must be located in the `aospy_path` directory and be the same name
    as the desired aospy module subpackage, namely one of `regions`, `calcs`,
    `variables`, and `projects`.
    """
    import imp
    return imp.load_source(
        name, '/'.join([user_path, name, '__init__.py']).replace('//', '/')
    )


def robust_bool(obj):
    try:
        return bool(obj)
    except ValueError:
        return obj.any()


def get_parent_attr(obj, attr, strict=False):
    """
    Check if the object has the given attribute and it is non-empty.  If not,
    check each parent object for the attribute and use the first one found.
    """
    attr_val = getattr(obj, attr, False)
    if robust_bool(attr_val):
        return attr_val

    else:
        for parent in ('parent', 'var', 'run', 'model', 'proj'):
            parent_obj = getattr(obj, parent, False)
            if parent_obj:
                return get_parent_attr(parent_obj, attr, strict=strict)

        if strict:
            raise AttributeError('Attribute %s not found in parent of %s'
                                 % (attr, obj))
        else:
            return None


def dict_name_keys(objs):
    """Create dict whose keys are the 'name' attr of the objects."""
    assert isinstance(objs, (tuple, list, dict))
    if isinstance(objs, (tuple, list)):
        try:
            return {obj.name: obj for obj in objs}
        except AttributeError:
            raise AttributeError
    else:
        return objs


def to_radians(field):
    if np.max(np.abs(field)) > 2*np.pi:
        return np.deg2rad(field)
    else:
        return field


def to_pascal(field):
    # For dp fields, this won't work if the input data is already Pascals and
    # the largest level thickness is < 1200 Pa, i.e. 12 hPa.  This will almost
    # never come up in practice for data interpolated to pressure levels, but
    # could come up in sigma data if model has sufficiently high vertical
    # resolution.
    if np.max(np.abs(field)) < 1200.:
        field *= 100.
    return field


def to_hpa(field):
    """Convert pressure array from Pa to hPa (if needed)."""
    if np.max(np.abs(field)) > 1200.:
        field /= 100.
    return field


def level_thickness(p):
    """
    Calculates the thickness, in Pa, of each pressure level.

    Assumes that the pressure values given are at the center of that model
    level, except for the lowest value (typically 1000 hPa), which is the
    bottom boundary. The uppermost level extends to 0 hPa.

    """
    # Bottom level extends from p[0] to halfway betwen p[0]
    # and p[1].
    p = to_pascal(p)
    dp = [0.5*(p[0] - p[1])]
    # Middle levels extend from halfway between [k-1], [k] and [k], [k+1].
    for k in range(1, p.size-1):
        dp.append(0.5*(p[k-1] - p[k+1]))
    # Top level extends from halfway between top two levels to 0 hPa.
    dp.append(0.5*(p[-2] + p[-1]))
    # Convert to numpy array and from hectopascals (hPa) to Pascals (Pa).
    return xray.DataArray(dp, coords=[p/100.0], dims=['level'])


def phalf_from_sigma(bk, pk, ps):
    """
    This should work.
    """
    return (ps*bk + pk)


def pfull_from_phalf(phalf, pfull_coord):
    """
    Compute data at full sigma levels from the values at the half levels.
    """
    # We will need to be smart in how we set the coordinates so that we can
    # add things gracefully within xray.
    phalf_top = phalf.isel(phalf=slice(1,None))
    phalf_top = phalf_top.rename({'phalf' : 'pfull'})
    phalf_top['pfull'] = pfull_coord

    phalf_bot = phalf.isel(phalf=slice(None,-1))
    phalf_bot = phalf_bot.rename({'phalf' : 'pfull'})
    phalf_bot['pfull'] = pfull_coord

    return 0.5*(phalf_bot + phalf_top)


def phalf_from_pfull(pfull, val_toa=0, val_sfc=0):
    """
    Compute data at half sigma levels from the values at full levels, given the
    specified top and bottom boundary conditions.

    Could be the pressure array itself, but it could also be any other data
    defined at pressure levels.
    """
    phalf = np.empty((pfull.shape[0] + 1, pfull.shape[1], pfull.shape[2]))
    phalf[0] = val_toa
    phalf[-1] = val_sfc
    phalf[1:-1] = 0.5*(pfull[:-1] + pfull[1:])
    return phalf


def pfull_from_sigma(bk, pk, ps, pfull_coord):
    return pfull_from_phalf(phalf_from_sigma(bk, pk, ps), pfull_coord)


def dp_from_phalf(phalf, pfull_coord):
    # We need to make sure dp is on a pfull coord.
    dp = phalf.diff(dim='phalf', n=1)
    dp = dp.rename({'phalf': 'pfull'})
    dp['pfull'] = pfull_coord
    return dp


def dp_from_sigma(bk, pk, ps, pfull_coord):
    return dp_from_phalf(phalf_from_sigma(bk, pk, ps), pfull_coord)


def weight_by_delta(integrand, delta):
    """Multiply an xray.DataArray by some weights."""
    return integrand*delta


def integrate(integrand, delta, dim):
    """Integrate along the given dimension."""
    prod = weight_by_delta(integrand, delta)
    return prod.sum(dim=dim)


def int_dp_g(integrand, dp):
    """
    Mass weighted integral.
    """
    return integrate(integrand, dp, vert_coord_name(dp)) * (1. / grav.value)


def vert_coord_name(dp):
    for name in ['level', 'pfull']:
        if name in dp.coords:
            return name
    return None


def dp_from_p(p, ps):
    """Get level thickness of pressure data, incorporating surface pressure."""
    # Top layer goes to 0 hPa; bottom layer goes to 1100 hPa.
    p = to_pascal(p)[np.newaxis,:,np.newaxis,np.newaxis]
    p_top = np.array([0])[np.newaxis,:,np.newaxis,np.newaxis]
    p_bot = np.array([1.1e5])[np.newaxis,:,np.newaxis,np.newaxis]

    # Layer edges are halfway between the given pressure levels.
    p_edges_interior = 0.5*(p[:,:-1] + p[:,1:])
    p_edges = np.concatenate((p_bot, p_edges_interior, p_top), axis=1)
    p_edge_above = p_edges[:, 1:]
    p_edge_below = p_edges[:, :-1]
    dp_interior = p_edge_below - p_edge_above

    ps = to_pascal(ps)[:,np.newaxis,:,:]
    # If ps < p_edge_below, then ps becomes the layer's bottom boundary.
    dp_adj_sfc = ps - p_edge_above[np.newaxis,:,np.newaxis,np.newaxis]
    dp = np.where(np.sign(ps - p_edge_below) > 0, dp_interior, dp_adj_sfc)
    # Mask where ps is less than the p.
    return np.ma.masked_where(ps < p, dp)
