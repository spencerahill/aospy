"""Utility functions for dealing with vertical coordinates."""
import logging

import numpy as np
import xarray as xr

from .._constants import GRAV_EARTH
from ..var import Var
from .. import internal_names


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
            logging.debug(warn_msg)
            return np.deg2rad(arr)
    # Otherwise, assume degrees if the values are sufficiently large.
    threshold = 0.1*np.pi if is_delta else 4*np.pi
    if np.max(np.abs(arr)) > threshold:
        warn_msg = ("Conversion applied: degrees -> radians to array: "
                    "{}".format(arr))
        logging.debug(warn_msg)
        return np.deg2rad(arr)
    return arr


def to_pascal(arr, is_dp=False):
    """Force data with units either hPa or Pa to be in Pa."""
    threshold = 400 if is_dp else 1200
    if np.max(np.abs(arr)) < threshold:
        warn_msg = "Conversion applied: hPa -> Pa to array: {}".format(arr)
        logging.debug(warn_msg)
        return arr*100.
    return arr


def to_hpa(arr):
    """Convert pressure array from Pa to hPa (if needed)."""
    if np.max(np.abs(arr)) > 1200.:
        warn_msg = "Conversion applied: Pa -> hPa to array: {}".format(arr)
        logging.debug(warn_msg)
        return arr / 100.
    return arr


def phalf_from_ps(bk, pk, ps):
    """Compute pressure of half levels of hybrid sigma-pressure coordinates."""
    return ps*bk + pk


def replace_coord(arr, old_dim, new_dim, new_coord):
    """Replace a coordinate with new one; new and old must have same shape."""
    new_arr = arr.rename({old_dim: new_dim})
    new_arr[new_dim] = new_coord
    return new_arr


def to_pfull_from_phalf(arr, pfull_coord):
    """Compute data at full pressure levels from values at half levels."""
    phalf_top = arr.isel(**{internal_names.PHALF_STR: slice(1, None)})
    phalf_top = replace_coord(phalf_top, internal_names.PHALF_STR,
                              internal_names.PFULL_STR, pfull_coord)

    phalf_bot = arr.isel(**{internal_names.PHALF_STR: slice(None, -1)})
    phalf_bot = replace_coord(phalf_bot, internal_names.PHALF_STR,
                              internal_names.PFULL_STR, pfull_coord)
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
    d_deta = arr.diff(dim=internal_names.PHALF_STR, n=1)
    return replace_coord(d_deta, internal_names.PHALF_STR,
                         internal_names.PFULL_STR, pfull_coord)


def d_deta_from_pfull(arr):
    """Compute $\partial/\partial\eta$ of the array on full hybrid levels.

    $\eta$ is the model vertical coordinate, and its value is assumed to simply
    increment by 1 from 0 at the surface upwards.  The data to be differenced
    is assumed to be defined at full pressure levels.

    Parameters
    ----------
    arr : xarray.DataArray containing the 'pfull' dim

    Returns
    -------
    deriv : xarray.DataArray with the derivative along 'pfull' computed via
            2nd order centered differencing.
    """  # noqa: W605
    right = arr[{internal_names.PFULL_STR: slice(2, None, None)}].values
    left = arr[{internal_names.PFULL_STR: slice(0, -2, 1)}].values
    deriv = xr.DataArray(np.zeros(arr.shape), dims=arr.dims,
                         coords=arr.coords)
    deriv[{internal_names.PFULL_STR: slice(1, -1, 1)}] = (right - left) / 2.
    deriv[{internal_names.PFULL_STR: 0}] = (
        arr[{internal_names.PFULL_STR: 1}].values -
        arr[{internal_names.PFULL_STR: 0}].values)
    deriv[{internal_names.PFULL_STR: -1}] = (
        arr[{internal_names.PFULL_STR: -1}].values -
        arr[{internal_names.PFULL_STR: -2}].values)
    return deriv


def dp_from_ps(bk, pk, ps, pfull_coord):
    """Compute pressure level thickness from surface pressure"""
    return d_deta_from_phalf(phalf_from_ps(bk, pk, ps), pfull_coord)


def integrate(arr, ddim, dim=False, is_pressure=False):
    """Integrate along the given dimension."""
    if is_pressure:
        dim = vert_coord_name(ddim)
    return (arr*ddim).sum(dim=dim)


def get_dim_name(arr, names):
    """Determine if an object has an attribute name matching a given list."""
    for name in names:
        # TODO: raise warning/exception when multiple names arr attrs.
        if hasattr(arr, name):
            return name
    raise AttributeError("No attributes of the object `{0}` match the "
                         "specified names of `{1}`".format(arr, names))


def vert_coord_name(arr):
    return get_dim_name(arr, [internal_names.PLEVEL_STR,
                              internal_names.PFULL_STR])


def int_dp_g(arr, dp):
    """Mass weighted integral."""
    return integrate(arr, to_pascal(dp, is_dp=True),
                     vert_coord_name(dp)) / GRAV_EARTH


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
    p_str = get_dim_name(p, (internal_names.PLEVEL_STR, 'plev'))
    p_vals = to_pascal(p.values.copy())

    # Layer edges are halfway between the given pressure levels.
    p_edges_interior = 0.5*(p_vals[:-1] + p_vals[1:])
    p_edges = np.concatenate(([p_bot], p_edges_interior, [p_top]))
    p_edge_above = p_edges[1:]
    p_edge_below = p_edges[:-1]
    dp = p_edge_below - p_edge_above
    if not all(np.sign(dp)):
        raise ValueError("dp array not all > 0 : {}".format(dp))
    # Pressure difference between ps and the upper edge of each pressure level.
    p_edge_above_xr = xr.DataArray(p_edge_above, dims=p.dims, coords=p.coords)
    dp_to_sfc = ps - p_edge_above_xr
    # Find the level adjacent to the masked, under-ground levels.
    change = xr.DataArray(np.zeros(dp_to_sfc.shape), dims=dp_to_sfc.dims,
                          coords=dp_to_sfc.coords)
    change[{p_str: slice(1, None)}] = np.diff(
        np.sign(ps - to_pascal(p.copy()))
    )
    dp_combined = xr.DataArray(np.where(change, dp_to_sfc, dp),
                               dims=dp_to_sfc.dims, coords=dp_to_sfc.coords)
    # Mask levels that are under ground.
    above_ground = ps > to_pascal(p.copy())
    above_ground[p_str] = p[p_str]
    dp_with_ps = dp_combined.where(above_ground)
    # Revert to original dim order.
    possible_dim_orders = [
        (internal_names.TIME_STR, p_str, internal_names.LAT_STR,
         internal_names.LON_STR),
        (internal_names.TIME_STR, p_str, internal_names.LAT_STR),
        (internal_names.TIME_STR, p_str, internal_names.LON_STR),
        (internal_names.TIME_STR, p_str),
        (p_str, internal_names.LAT_STR, internal_names.LON_STR),
        (p_str, internal_names.LAT_STR),
        (p_str, internal_names.LON_STR),
        (p_str,),
    ]
    for dim_order in possible_dim_orders:
        try:
            return dp_with_ps.transpose(*dim_order)
        except ValueError:
            logging.debug("Failed transpose to dims: {}".format(dim_order))
    else:
        logging.debug("No transpose was successful.")
        return dp_with_ps


def level_thickness(p, p_top=0., p_bot=1.01325e5):
    """
    Calculates the thickness, in Pa, of each pressure level.

    Assumes that the pressure values given are at the center of that model
    level, except for the lowest value (typically 1000 hPa), which is the
    bottom boundary. The uppermost level extends to 0 hPa.

    Unlike `dp_from_p`, this does not incorporate the surface pressure.

    """
    p_vals = to_pascal(p.values.copy())
    dp_vals = np.empty_like(p_vals)
    # Bottom level extends from p[0] to halfway betwen p[0] and p[1].
    dp_vals[0] = p_bot - 0.5*(p_vals[0] + p_vals[1])
    # Middle levels extend from halfway between [k-1], [k] and [k], [k+1].
    dp_vals[1:-1] = 0.5*(p_vals[0:-2] - p_vals[2:])
    # Top level extends from halfway between top two levels to 0 hPa.
    dp_vals[-1] = 0.5*(p_vals[-2] + p_vals[-1]) - p_top
    dp = p.copy()
    dp.values = dp_vals
    return dp


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


bk = Var(
    name=internal_names.BK_STR,
    alt_names=internal_names.GRID_ATTRS[internal_names.BK_STR],
    def_vert=True,
    def_time=False,
    def_lon=False,
    def_lat=False
)


pk = Var(
    name=internal_names.PK_STR,
    alt_names=internal_names.GRID_ATTRS[internal_names.PK_STR],
    def_vert=True,
    def_time=False,
    def_lon=False,
    def_lat=False
)


pfull_coord = Var(
    name=internal_names.PFULL_STR,
    alt_names=internal_names.GRID_ATTRS[internal_names.PFULL_STR],
    def_vert=True,
    def_time=False,
    def_lon=False,
    def_lat=False
)


ps = Var(
    name='ps',
    domain='atmos',
    description='Surface pressure',
    units='Pa',
    def_vert=False,
    def_time=True,
    def_lon=True,
    def_lat=True
)


p_level = Var(
    name='p',
    alt_names=internal_names.GRID_ATTRS[internal_names.PLEVEL_STR],
    domain='atmos',
    description='Pressure on interpolated levels',
    units='Pa',
    def_vert=True,
    def_time=False,
    def_lon=False,
    def_lat=False
)


dp_level = Var(
    name='dp',
    description='Pressure thickness of model levels',
    units='Pa',
    def_vert=True,
    def_time=True,
    def_lon=True,
    def_lat=True,
    func=dp_from_p,
    variables=(p_level, ps)
)


p_eta = Var(
    name='p',
    description='Pressure at model-native level midpoints',
    units='Pa',
    def_vert=True,
    def_time=True,
    def_lon=True,
    def_lat=True,
    func=pfull_from_ps,
    variables=(bk, pk, ps, pfull_coord)
)


dp_eta = Var(
    name='dp',
    description='Pressure thickness of model levels',
    units='Pa',
    def_vert=True,
    def_time=True,
    def_lon=True,
    def_lat=True,
    func=dp_from_ps,
    variables=(bk, pk, ps, pfull_coord)
)
