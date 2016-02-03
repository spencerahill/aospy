"""aospy.utils: utility functions for the aospy module."""
import logging
import warnings

from infinite_diff import FiniteDiff
import numpy as np
import pandas as pd
import xarray as xr

from .__config__ import (PHALF_STR, PFULL_STR, PLEVEL_STR, TIME_STR,
                         LAT_STR, LON_STR, user_path)
from .constants import grav, Constant


def coord_to_new_dataarray(arr, dim):
    """Create a DataArray comprising the coord for the specified dim.

    Useful, for example, when wanting to resample in time, because at least
    for xarray 0.6.0 and prior, the `resample` method doesn't work when applied
    to coords.  The DataArray returned by this method lacks that limitation.
    """
    return xr.DataArray(arr[dim].values, coords=[arr[dim].values], dims=[dim])


# 2015-11-16 S. Hill: This time-related function should be moved to time module
def apply_time_offset(time, months=0, days=0, hours=0):
    """Apply the given offset to the given time array.

    This is useful for GFDL model output of instantaneous values.  For example,
    3 hourly data postprocessed to netCDF files spanning 1 year each will
    actually have time values that are offset by 3 hours, such that the first
    value is for 1 Jan 03:00 and the last value is 1 Jan 00:00 of the
    subsequent year.  This causes problems in xarray, e.g. when trying to group
    by month.  It is resolved by manually subtracting off those three hours,
    such that the dates span from 1 Jan 00:00 to 31 Dec 21:00 as desired.
    """
    return (pd.to_datetime(time.values) +
            pd.tseries.offsets.DateOffset(months=months, days=days,
                                          hours=hours))


# 2015-11-16 S. Hill: This time-related function should be moved to time module
def monthly_mean_ts(arr):
    """Convert a sub-monthly time-series into one of monthly means."""
    if isinstance(arr, (float, complex, int, Constant)):
        return arr
    if isinstance(arr, (list, tuple)):
        arr_type = type(arr)
        return arr_type([monthly_mean_ts(a) for a in arr])
    try:
        return arr.resample('1M', TIME_STR, how='mean').dropna(TIME_STR)
    except KeyError:
        raise KeyError("`{}` lacks time dimension with "
                       "label `{}`.".format(arr, TIME_STR))


# 2015-11-16 S. Hill: This time-related function should be moved to time module
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
    assert isinstance(objs, (tuple, list, dict, set))
    if isinstance(objs, (tuple, list, set)):
        try:
            return {obj.name: obj for obj in objs}
        except AttributeError as e:
            raise AttributeError(e)
    return objs


def to_radians(arr, is_delta=False):
    """Force data with units either degrees or radians to be radians."""
    # Infer the units from embedded metadata, if it's there.
    try:
        units = arr.units
    except AttributeError:
        pass
    else:
        if units.lower().startswith('degrees'):
            warn_msg = ("Conversion applied: degrees -> radians to array: "
                        "{}".format(arr))
            warnings.warn(warn_msg, UserWarning)
            return np.deg2rad(arr)
    # Otherwise, assume degrees if the values are sufficiently large.
    threshold = 0.1*np.pi if is_delta else 4*np.pi
    if np.max(np.abs(arr)) > threshold:
        warn_msg = ("Conversion applied: degrees -> radians to array: "
                    "{}".format(arr))
        warnings.warn(warn_msg, UserWarning)
        return np.deg2rad(arr)
    return arr


def to_pascal(arr, is_dp=False):
    """Force data with units either hPa or Pa to be in Pa."""
    threshold = 400 if is_dp else 1200
    if np.max(np.abs(arr)) < threshold:
        warn_msg = "Conversion applied: hPa -> Pa to array: {}".format(arr)
        warnings.warn(warn_msg, UserWarning)
        return arr*100.
    return arr


def to_hpa(arr):
    """Convert pressure array from Pa to hPa (if needed)."""
    if np.max(np.abs(arr)) > 1200.:
        warn_msg = "Conversion applied: Pa -> hPa to array: {}".format(arr)
        warnings.warn(warn_msg, UserWarning)
        return arr / 100.
    return arr


def phalf_from_ps(bk, pk, ps):
    """Compute pressure of half levels of hybrid sigma-pressure coordinates."""
    return (ps*bk + pk)


def replace_coord(arr, old_dim, new_dim, new_coord):
    """Replace a coordinate with new one; new and old must have same shape."""
    new_arr = arr.rename({old_dim: new_dim})
    # new_arr[new_dim] = new_coord
    # 2016-01-25 S. Hill: Temporary workaround to deal with xarray 0.7.0 bug.
    # new_arr[new_dim].values = new_coord.values
    ds = new_arr.to_dataset(name='new_arr')
    ds[new_dim] = new_coord
    return ds['new_arr']


def to_pfull_from_phalf(arr, pfull_coord):
    """Compute data at full pressure levels from values at half levels."""
    phalf_top = arr.isel(**{PHALF_STR: slice(1, None)})
    phalf_top = replace_coord(phalf_top, PHALF_STR, PFULL_STR, pfull_coord)

    phalf_bot = arr.isel(**{PHALF_STR: slice(None, -1)})
    phalf_bot = replace_coord(phalf_bot, PHALF_STR, PFULL_STR, pfull_coord)
    return 0.5*(phalf_bot + phalf_top)


def to_phalf_from_pfull(arr, val_toa=0, val_sfc=0):
    """Compute data at half pressure levels from values at full levels.

    Could be the pressure array itself, but it could also be any other data
    defined at pressure levels.  Requires specification of values at surface
    and top of atmosphere.
    """
    phalf = np.zeros((arr.shape[0] + 1, arr.shape[1], arr.shape[2]))
    phalf[0] = val_toa
    phalf[-1] = val_sfc
    phalf[1:-1] = 0.5*(arr[:-1] + arr[1:])
    return phalf


def pfull_from_ps(bk, pk, ps, pfull_coord):
    """Compute pressure at full levels from surface pressure."""
    return to_pfull_from_phalf(phalf_from_ps(bk, pk, ps), pfull_coord)


def d_deta_from_phalf(arr, pfull_coord):
    """Compute pressure level thickness from half level pressures."""
    d_deta = arr.diff(dim=PHALF_STR, n=1)
    return replace_coord(d_deta, PHALF_STR, PFULL_STR, pfull_coord)


def d_deta_from_pfull(arr):
    """Compute $\partial/\partial\eta$ of the array on full hybrid levels.

    $\eta$ is the model vertical coordinate, and its value is assumed to simply
    increment by 1 from 0 at the surface upwards.  The data to be differenced
    is assumed to be defined at full pressure levels.
    """
    deriv = FiniteDiff.cen_diff(arr, PFULL_STR, do_edges_one_sided=True) / 2.
    # Edges use 1-sided differencing, so only spanning one level, not two.
    deriv[{PFULL_STR: 0}] = deriv[{PFULL_STR: 0}] * 2.
    deriv[{PFULL_STR: -1}] = deriv[{PFULL_STR: -1}] * 2.
    return deriv


def dp_from_ps(bk, pk, ps, pfull_coord):
    """Compute pressure level thickness from surface pressure"""
    return d_deta_from_phalf(phalf_from_ps(bk, pk, ps), pfull_coord)


def integrate(arr, ddim, dim=False, is_pressure=False):
    """Integrate along the given dimension."""
    if is_pressure:
        dim = vert_coord_name(ddim)
    return (arr*ddim).sum(dim=dim)


def vert_coord_name(dp):
    for name in [PLEVEL_STR, PFULL_STR]:
        if name in dp.coords:
            return name
    return


def int_dp_g(arr, dp):
    """Mass weighted integral."""
    return integrate(arr, to_pascal(dp, is_dp=True),
                     vert_coord_name(dp)) / grav.value


def dp_from_p(p, ps, p_top=0., p_bot=1.1e5):
    """Get level thickness of pressure data, incorporating surface pressure.

    Level edges are defined as halfway between the levels, as well as the user-
    specified uppermost and lowermost values.  The dp of levels whose bottom
    pressure is less than the surface pressure is not changed by ps, since they
    don't intersect the surface.  If ps is in between a level's top and bottom
    pressures, then its dp becomes the pressure difference between its top and
    ps.  If ps is less than a level's top and bottom pressures, then that level
    is underground and its values are masked.

    Note that postprocessing routines (e.g. at GFDL) typically mask out data
    wherever the surface pressure is less than the level's given value, not the
    level's upper edge.  This masks out more levels than the

    """
    p_vals = to_pascal(p.values)

    # Layer edges are halfway between the given pressure levels.
    p_edges_interior = 0.5*(p_vals[:-1] + p_vals[1:])
    p_edges = np.concatenate(([p_bot], p_edges_interior, [p_top]))
    p_edge_above = p_edges[1:]
    p_edge_below = p_edges[:-1]
    dp = p_edge_below - p_edge_above
    if not all(np.sign(dp)):
        raise ValueError("dp array not all > 0 : {}".format(dp))
    # Pressure difference between ps and the upper edge of each pressure level.
    p_edge_above_xarray = xr.DataArray(p_edge_above, dims=p.dims,
                                       coords=p.coords)
    dp_to_sfc = ps - p_edge_above_xarray
    # Find the level adjacent to the masked, under-ground levels.
    change = xr.DataArray(np.zeros(dp_to_sfc.shape), dims=dp_to_sfc.dims,
                            coords=dp_to_sfc.coords)
    change[{PLEVEL_STR: slice(1, None)}] = np.diff(
        np.sign(ps - to_pascal(p.copy()))
    )
    dp_combined = xr.DataArray(np.where(change, dp_to_sfc, dp),
                                 dims=dp_to_sfc.dims, coords=dp_to_sfc.coords)
    # Mask levels that are under ground.
    above_ground = ps > to_pascal(p.copy())
    above_ground[PLEVEL_STR].values = p.values
    dp_with_ps = dp_combined.where(above_ground)
    # Revert to original dim order.
    possible_dim_orders = [
        (TIME_STR, PLEVEL_STR, LAT_STR, LON_STR),
        (TIME_STR, PLEVEL_STR, LAT_STR),
        (TIME_STR, PLEVEL_STR, LON_STR),
        (TIME_STR, PLEVEL_STR),
        (PLEVEL_STR, LAT_STR, LON_STR),
        (PLEVEL_STR, LAT_STR),
        (PLEVEL_STR, LON_STR),
        (PLEVEL_STR,),
    ]
    for dim_order in possible_dim_orders:
        try:
            return dp_with_ps.transpose(*dim_order)
        except ValueError:
            logging.debug("Failed transpose to dims: {}".format(dim_order))
    else:
        logging.debug("No transpose was successful.")
        return dp_with_ps


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
    return xr.DataArray(dp, coords=[p/100.0], dims=['level'])


def does_coord_increase_w_index(arr):
    """Determine if the array values increase with the index.

    Useful, e.g., for pressure, which sometimes is indexed surface to TOA and
    sometimes the opposite.
    """
    diff = np.diff(arr)
    if not np.all(np.abs(np.sign(diff))):
        raise ValueError("Array is not monotonic: {}".format(arr))
    # Since we know its monotonic, just test the first value.
    return bool(diff[0])
