"""Calcs submodule of aospy module.
"""
from hadley_itcz import *
from vert_fluxes import *
from trop_height import trop_height

### Utility functions for sigma coordinates. ###

def phalf_from_sigma(bk, pk, ps):
    """
    Compute pressure at sigma half levels from the sigma coordinate arrays and
    the surface pressure.

    Assume pk, bk, and ps are in Pa, unitless, and Pa, respectively.  Assume
    pk[-1] and bk[-1] are at surface where pressure equals ps and pk[0] and
    bk[0] are at top of atmosphere where pressure equals zero.  pk and bk are
    1-d arrays; ps has last two dimensions (lat, lon) but may also have time as
    first dimension.  Assume pk and bk include both endpoints, i.e. the surface
    and TOA, such that the number of sigma layers is one less than the length
    of either pk or bk.
    """
    import numpy as np
    # 3D ps array assumed to be (time, lat, lon).
    if ps.ndim in (3, 4):
        bk = bk[np.newaxis,:,np.newaxis,np.newaxis]
        pk = pk[np.newaxis,:,np.newaxis,np.newaxis]
        if ps.ndim == 3:
            ps = ps[:,np.newaxis,:,:]
    # 2D ps array assumed to be (lat, lon).
    elif ps.ndim == 2:
        bk = bk[:,np.newaxis,np.newaxis]
        pk = pk[:,np.newaxis,np.newaxis]
        ps = ps[np.newaxis,:,:]
    return np.squeeze(pk + ps*bk)

def pfull_from_phalf(phalf):
    """
    Compute data at full sigma levels from the values at half levels.

    Could be the pressure array itself, but it could also be any other data 
    defined at half levels.
    """
    # 4D array assumed to be (time, p, lat, lon).
    if phalf.ndim == 4:
        return 0.5*(phalf[:,1:] + phalf[:,:-1])
    # Anything else assumed to have p as first dimension.
    else:
        return 0.5*(phalf[1:] + phalf[:-1])

def phalf_from_pfull(pfull, val_toa=0, val_sfc=0):
    """
    Compute data at half sigma levels from the values at full levels, given the 
    specified top and bottom boundary conditions.

    Could be the pressure array itself, but it could also be any other data 
    defined at pressure levels.
    """
    import numpy as np
    phalf = np.empty((pfull.shape[0] + 1, pfull.shape[1], pfull.shape[2]))
    phalf[0] = val_toa
    phalf[-1] = val_sfc
    phalf[1:-1] = 0.5*(pfull[:-1] + pfull[1:])
    return phalf

def pfull_from_sigma(bk, pk, ps):
    """
    Compute pressure at full sigma levels from the sigma coordinate arrays and
    surface pressure.
    """
    phalf = phalf_from_sigma(bk, pk, ps)
    return pfull_from_phalf(phalf)

def dp_from_phalf(phalf):
    """Compute pressure-depth of vertical levels from level edge pressures."""
    # If 4D, assume dimensions (time, p, lat, lon).
    if phalf.ndim == 4:
        return phalf[:,1:] - phalf[:,:-1]
    # Otherwise assume first dimension is p.
    else:
        return phalf[1:] - phalf[:-1]

def dp_from_sigma(bk, pk, ps):
    """Compute sigma layer pressure thickness."""
    import numpy as np
    phalf = phalf_from_sigma(bk, pk, ps)
    return dp_from_phalf(phalf)

### x, dx from longitude, y, dy from latitude

# def x_from_latlon(lat, lon, radius):
#     """Create array of x-distances based on arrays of lons and lats."""
#     import numpy as np
#     rad = radius*np.abs(np.cos(np.deg2rad(lat)))[np.newaxis,:]
#     return np.squeeze(rad*np.deg2rad(lon[:,np.newaxis]))

# def y_from_latlon(lat, lon, radius):
#     """Create an array of y-distances from an array of lats."""
#     import numpy as np
#     return radius*np.tile(np.sin(np.deg2rad(lat))[:,np.newaxis], len(lon)).T

# def dx_from_latlon(lat, lon, radius):
#     """Compute zonal spacing of grid cell centers."""
#     import numpy as np
#     x = x_from_latlon(lon, lat, radius)
#     dx = np.empty(x.shape)
#     dx[:,:-1] = np.diff(x, n=1, axis=-1)
#     # Assume lon has smallest value in 1st index.
#     assert np.argmin(x) == 0
#     dx[:,-1] = (x[:,0] + radius*2*np.pi) - x[:,-1]
#     return dx

# def dy_from_latlon(lat, lon, radius):
#     """Compute meridional spacing of grid cell centers."""
#     y = y_from_latlon(lat, lon, radius)
#     return y[1:] - y[:-1]

# def dlon_from_latlon(lat, lon):
    # """Compute longitude spacing of grid cell centers."""
    # import numpy as np
    # return 

# def dlat_from_latlon(lat, lon):
    # """Compute latitude spacing of grid cell centers."""
    # import numpy as np
    # return 

### General finite differencing.
def fwd_diff(f, dx):
    """1st order accurate forward differencing."""
    return (f[1:] - f[:-1]) / dx

def cen_diff2(f, dx, dx_array=False):
    """2nd order accurate centered differencing."""
    if dx_array:
        df_dx = (f[2:] - f[:-2]) / (dx[:-1] + dx[:1])
    else:
        df_dx = (f[2:] - f[:-2]) / (2.*dx)
    return df_dx

def cen_diff4(f, dx, dx_array=False):
    """4th order accurate centered differencing."""
    if dx_array:
        # Assume len(dx) == len(f) - 1.
        d1 = (f[3:-1] - f[1:-3]) / (dx[1:-2] + dx[2:-1])
        d2 = (f[4:] - f[:-4]) / (dx[:-3] + dx[1:-2] + dx[2:-1] + dx[3:])
        return (4*d1 - d2) / 3.
    else:
        return (8*(f[3:-1] - f[1:-3]) - (f[4:] - f[:-4])) / (12.*dx)
    return numer / denom

### Derivatives in x, y, and p.

def to_rad(field):
    import numpy as np
    if np.max(np.abs(field)) > 2*np.pi:
        return np.deg2rad(field)
    else:
        return field

def to_pascal(field):
    import numpy as np
    # For dp fields, this won't work if the input data is already Pascals and
    # the largest level thickness is < 1200 Pa, i.e. 12 hPa.  This will almost
    # never come up in practice for data interpolated to pressure levels, but
    # could come up in sigma data if model has sufficiently high vertical
    # resolution.
    if np.max(np.abs(field)) < 1200.:
        field *= 100.
    return field

def d_dx_from_latlon(field, lat, lon, radius):
    """Compute \partial(field)/\partial x using centered differencing."""
    import numpy as np
    # Verify longitude spacing is uniform.
    dlon = lon[1] - lon[0]
    assert np.allclose(360. - (lon[-1] - lon[0]), dlon)
    assert np.allclose(lon[2:] - lon[1:-1], lon[1:-1] - lon[:-2])
    lat, lon = to_rad(lat), to_rad(lon)
    dlon = lon[1] - lon[0]
    prefactor = 1. / (radius*np.cos(lat))[:,np.newaxis]
    # Assume latitude and longitude are last two axes.
    # Transpose the arrays: lon is 1st axis, lat is 2nd.
    f = field.T
    # Wrap around at 0/360 degrees longitude.
    f = np.ma.concatenate((f[-2:], f, f[:2]), axis=0)
    df_dx = cen_diff4(f, dlon)
    # Transpose again to regain original axis order.
    return prefactor*df_dx.T

def d_dy_from_latlon(field, lat, lon, radius, vec_field=False):
    """Compute \partial(field)/\partial y using centered differencing."""
    import numpy as np
    lat = to_rad(lat)
    dlat = lat[1:] - lat[:-1]
    # Assume latitude and longitude are last two axes.
    if field.ndim == 2:
        lat = lat[np.newaxis,:]
        dlat = dlat[np.newaxis,:]
    if field.ndim == 3:
        lat = lat[np.newaxis,:,np.newaxis]
        dlat = dlat[np.newaxis,:,np.newaxis]
    elif field.ndim == 4:
        lat = lat[np.newaxis,:,np.newaxis,np.newaxis]
        dlat = dlat[np.newaxis,:,np.newaxis,np.newaxis]
    f = field.T
    prefactor = 1. / radius
    # Del operator differs in spherical coords for scalar v. vector fields.
    if vec_field:
        f *= np.cos(lat)
        prefactor /= np.cos(lat)
        prefactor = prefactor.T
    # Roll lat to be leading axis for broadcasting.
    f = np.rollaxis(f, 1, 0)
    dlat = np.rollaxis(dlat, 1, 0)
    df_dy = np.ma.empty(f.shape)
    df_dy[2:-2] = cen_diff4(f, dlat, dx_array=True)
    df_dy[1]    = cen_diff2(f[:3], dlat[:2], dx_array=True)
    df_dy[-2]   = cen_diff2(f[-3:], dlat[-2:], dx_array=True)
    df_dy[0]    = fwd_diff(f[:2], dlat[0])
    df_dy[-1]   = fwd_diff(f[-2:], dlat[-1])
    # Roll axis and transpose again to regain original axis order.
    df_dy = np.rollaxis(df_dy, 0, 2)
    return prefactor*df_dy.T

def d_dp_from_p(field, p):
    import numpy as np
    """Derivative in pressure of a given field."""
    # Assume pressure is 3rd to last axis: ([time,] p, lat, lon)
    f = field.T
    p = to_pascal(p)
    # Amend array dimensions as necessary for broadcasting purposes.
    if p.ndim == 1:
        p = p[np.newaxis,np.newaxis,:,np.newaxis]
    elif p.ndim == 3:
        p = p[np.newaxis,:,:,:].T
    else:
        p = p.T
    # One-sided difference at TOA and surface; centered difference elsewhere.
    df_dp = np.ma.empty(f.shape)
    df_dp[:,:,1:-1] = (f[:,:,2:] - f[:,:,:-2]) / (p[:,:,2:] - p[:,:,:-2])
    df_dp[:,:,0]    = (f[:,:,1]  - f[:,:,0])   / (p[:,:,1]  - p[:,:,0])
    df_dp[:,:,-1]   = (f[:,:,-1] - f[:,:,-2])  / (p[:,:,-1] - p[:,:,-2])
    # Transpose again to regain original axis order.
    return df_dp.T

# def add_ps_value_to_field(field, p, ps, field_at_ps):
#     """Add value at surface pressure to pressure-defined data.
#
#     Result can be used to easily perform derivatives or integrals in pressure
#     that extend all the way to the surface.
#
#     Algorithm is as follows:
#     1) Determine p_bot, where p_bot is the largest p such that p < ps.
#     2) Unmask the element directly below p_bot.
#     3) Replace that element with the specified value at ps.
#     """
#     import numpy as np
#     if p.ndim == 1:
#         p = p[np.newaxis,:,np.newaaxis,np.newaxis]
#     p_bot = np.ma.where(p < ps)
#     return

### Integrals in x, y, and p.

def weight_by_delta(integrand, delta):
    """
    Weight the integrand by delta, usually for subsequent integration.

    delta array may be one dimension or three; in the latter it is assumed to 
    be of shape (vertical, lat, lon).  integrand is assumed to be 3 or 4
    dimensions, with time 1st if 4-D.  Both are assumed to be numpy arrays.
    """
    import numpy as np
    try:
        intdel = integrand*delta
    except ValueError:
        delta = delta[np.newaxis,:,np.newaxis, np.newaxis]
        intdel = integrand*delta
    return intdel

def integrate(integrand, delta, axis):
    """Integrate the array along the given axis using the given delta array."""
    import numpy as np
    prod = weight_by_delta(integrand, delta)
    return np.ma.sum(prod, axis=axis)

# def int_dlon(integrand, lat, lon, start=0., end=360.):
    # """Integrate in longitude."""
    # import numpy as np
    # # Assume pressure is 3rd to last axis.
    # inds = np.where((lon >= start) & (lon <= end))
    # dlon = dlon_from_latlon(lat, lon)
    # dlon = dlon[:,inds]
    # try:
        # integrand = integrand[inds]
    # except IndexError:
        # integrand = integrand[:,inds]
    # return integrate(integrand, dlon, -3)

# def int_dlat(integrand, lat, lon, start=-90., end=-90.):
    # """Integrate in latitude."""
    # return

def int_dp_g(integrand, dp, start=0., end=None):
    """Integrate vertically in pressure."""
    # Assume pressure is 3rd to last axis.
    import numpy as np
    from aospy.constants import grav
    dp = to_pascal(dp)
    return integrate(integrand, dp, -3) / grav

### Horizontal & vertical advection, divergence, and flux divergence functions.

def horiz_advec(field, u, v, lat, lon, radius):
    """Horizontal advection of the given field."""
    df_dx = d_dx_from_latlon(field, lat, lon, radius)
    df_dy = d_dy_from_latlon(field, lat, lon, radius, vec_field=False)
    return u*df_dx + v*df_dy

def vert_advec(field, omega, p):
    """Vertical advection of the given field."""
    return omega*d_dp_from_p(field, p)

def horiz_divg(u, v, lat, lon, radius):
    """Flow horizontal divergence."""
    du_dx = d_dx_from_latlon(u, lat, lon, radius)
    dv_dy = d_dy_from_latlon(v, lat, lon, radius, vec_field=True)
    return du_dx + dv_dy

def horiz_divg_from_windspharm(u, v):
    """Horizontal divergence computed using the windspharm package.

    The windspharm package can't handle masked data, so the procedure is as
    follows:

    1. Unmask the data, setting masked gridpoints to 0.
    2. Reshape the array, moving lat and lon to the first two dimensions
       and combining time and level into a single 3rd dimension.
    3. Perform the windspharm computations.
    4. Reshape the array back to its original shape.
    5. Re-apply the original mask and return the resulting array.

    Windspharm can handle 3-D data, but lat and lon have to be the first two
    dimensions.  So we have to reshape the array, moving time and/or level
    to be after lat and lon, and then reshape the array back.
    through the vertical and time indices.

    All steps except #3 are independent of the windspharm function being called.
    So we should have a single function that performs all steps, taking the
    function to be called in #3 as an argument.
    """
    import numpy as np
    from windspharm.standard import VectorWind
    w = VectorWind(u, v)
    return w.divergence()

def vert_divg(omega, p):
    """Flow vertical divergence."""
    return d_dp_from_p(omega, p)

def field_vert_int_bal(field, dp):
    """Impose vertical balance to the field, i.e. column integral = 0.

    Most frequently used for mass flux divergence to impose mass balance.
    """
    import numpy as np
    pos = np.ma.where(field > 0, field, 0)
    neg = np.ma.where(field < 0, field, 0)
    pos_int = int_dp_g(pos, dp)[:,np.newaxis,:,:]
    neg_int = int_dp_g(neg, dp)[:,np.newaxis,:,:]
    return pos - (pos_int/neg_int)*neg

def horiz_divg_mass_bal(u, v, lat, lon, radius, dp):
    """Horizontal divergence with column mass-balance correction applied."""
    from aospy import calc_levs_thick
    div = horiz_divg(u, v, lat, lon, radius)
    return field_vert_int_bal(div, dp)

def vert_divg_mass_bal(omega, p, dp):
    """Vertical divergence with column mass-balance correction applied."""
    from aospy import calc_levs_thick
    div = vert_divg(omega, p)
    return field_vert_int_bal(div, dp)

def divg_3d(u, v, omega, lat, lon, p):
    """Total (3-D) divergence.  Should equal 0 by continuity."""
    horiz = horiz_divg(u, v, lat, lon, radius)
    vert = vert_divg(omega, p)
    return horiz + vert

def field_times_horiz_divg(field, u, v, lat, lon, radius):
    """Scalar field times horizontal convergence."""
    return field*horiz_divg(u, v, lat, lon, radius)

def field_times_horiz_divg_mass_bal(field, u, v, lat, lon, radius, dp):
    """Scalar field times horizontal convergence."""
    return field*horiz_divg_mass_bal(u, v, lat, lon, radius, dp)

def field_times_vert_divg_mass_bal(field, omega, p, dp):
    """Scalar field times vertical convergence."""
    return field*vert_divg_mass_bal(omega, p, dp)

def field_horiz_flux_divg(field, u, v, lat, lon, radius):
    """Horizontal flux divergence of given scalar field."""
    dfu_dx = d_dx_from_latlon(u*field, lat, lon, radius)
    dfv_dy = d_dy_from_latlon(v*field, lat, lon, radius, vec_field=True)
    return dfu_dx + dfv_dy

def field_vert_flux_divg(field, omega, p):
    """Vertical flux divergence of a scalar field."""
    return d_dp_from_p(omega*field, p)

def field_horiz_advec_divg_sum(field, u, v, lat, lon, radius, dp):
    return (field_times_horiz_divg_mass_bal(field, u, v, lat, lon, radius, dp)
            + horiz_advec(field, u, v, lat, lon, radius))

def field_vert_advec_divg_sum(field, omega, p, dp):
    return (field_times_vert_divg_mass_bal(field, omega, p, dp) +
            vert_advec(field, omega, p, dp))

### Advection, divergence etc. applied to specific fields.

def mse_horiz_flux_divg(temp, hght, sphum, u, v, lat, lon, radius):
    """Horizontal flux convergence of moist static energy."""
    return field_horiz_flux_divg(mse(temp, hght, sphum),
                                 u, v, lat, lon, radius)

def mse_horiz_advec(temp, hght, sphum, u, v, lat, lon, radius):
    """Horizontal advection of moist static energy."""
    return horiz_advec(mse(temp, hght, sphum), u, v, lat, lon, radius)

def mse_times_horiz_divg(temp, hght, sphum, u, v, lat, lon, radius, dp):
    """Horizontal divergence of moist static energy."""
    return field_times_horiz_divg_mass_bal(
        mse(temp, hght, sphum), u, v, lat, lon, radius, dp
    )

def mse_horiz_advec_divg_sum(T, gz, q, u, v, lat, lon, rad, dp):
    return field_horiz_advec_divg_sum(mse(T, gz, q), u, v, lat, lon, rad, dp)

def mse_vert_flux_divg(T, gz, q, omega, p):
    """Vertical divergence times moist static energy."""
    return field_vert_flux_divg(mse(T, gz, q), omega, p)

def mse_vert_advec(temp, hght, sphum, omega, p):
    """Vertical advection of moist static energy."""
    return vert_advec(mse(temp, hght, sphum), omega, p)

def mse_times_vert_divg(T, gz, q, omega, p, dp):
    """Vertical divergence times moist static energy."""
    return field_times_vert_divg_mass_bal(mse(T, gz, q), omega, p, dp)

def dse_horiz_flux_divg(temp, hght, u, v, lat, lon, radius):
    """Horizontal flux convergence of moist static energy."""
    return field_horiz_flux_divg(dse(temp, hght), u, v, lat, lon, radius)

def dse_horiz_advec(temp, hght, u, v, lat, lon, radius):
    """Horizontal advection of moist static energy."""
    return horiz_advec(dse(temp, hght), u, v, lat, lon, radius)

def dse_times_horiz_divg(temp, hght, u, v, lat, lon, radius, dp):
    """Horizontal divergence of moist static energy."""
    return field_times_horiz_divg_mass_bal(
        dse(temp, hght), u, v, lat, lon, radius, dp
    )

def dse_horiz_advec_divg_sum(T, gz, u, v, lat, lon, rad, dp):
    return field_horiz_advec_divg_sum(dse(T, gz), u, v, lat, lon, rad, dp)

def dse_vert_advec(temp, hght, omega, p):
    """Vertical advection of moist static energy."""
    return vert_advec(dse(temp, hght), omega, p)

def q_horiz_flux_divg(q, u, v, lat, lon, radius):
    """Horizontal flux convergence of specific humidity."""
    return field_horiz_flux_divg(q, u, v, lat, lon, radius)

def q_horiz_advec(q, u, v, lat, lon, radius):
    """Horizontal advection of specific humidity."""
    return horiz_advec(q, u, v, lat, lon, radius)

def q_times_horiz_divg(q, u, v, lat, lon, radius, dp):
    """Horizontal divergence of specific humidity."""
    return field_times_horiz_divg_mass_bal(q, u, v, lat, lon, radius, dp)

def q_vert_advec(q, omega, p):
    """Vertical advection of specific humidity."""
    return vert_advec(q, omega, p)

def qu(sphum, u):
    """"Zonal moisture flux."""
    return sphum*u

def qv(sphum, v):
    """Meridional moisture flux."""
    return sphum*v

### Eddy computations.

def covariance(array1, array2, axis=None, weights=None):
    """
    Average `covariance` along the specified axis of two arrays.
    
    :param array1, array2: The two arrays to compute the covariance of.
    :type array1 array2: numpy.ndarray or numpy.ma.core.MaskedArray
    :param axis: Array axis number over which average is taken.
    :type axis: int or None
    :param weights: Weights used to perform the average.
    :type weights: int or None
    """
    from numpy import ma
    prod = ma.multiply(array1, array2)
    return ma.average(prod, axis=axis, weights=weights)

def eddy_component(array, axis=None, weights=None):
    """Compute the deviation from the average along the specified axis.

    :param array: The array on which to compute the eddy component.
    :type array: Numpy array.
    :param axis: The index of the axis on which to compute.
    :type axis: int
    :param weights: An array of weights used to compute the average.
    :type weights: Numpy array or `None`
    """
    from numpy import ma
    avg = ma.average(array, axis=axis, weights=weights)
    return ma.subtract(array, avg)

def eddy_covar_avg(array1, array2, axis=None, weights=None):
    """Compute the average eddy covariance of two fields.

    :param array1, array2: The two arrays to compute the covariance of.
    :type array1 array2: numpy.ndarray or numpy.ma.core.MaskedArray
    :param axis: Array axis number over which average is taken.
    :type axis: int or None
    :param weights: Weights used to perform the average.
    :type weights: int or None
    """
    from numpy import ma
    cov1 = eddy_component(array1, axis=axis, weights=weights)
    cov2 = eddy_component(array2, axis=axis, weights=weights)
    return covariance(cov1, cov2, axis=axis, weights=weights)

### Thermodynamic functions.

def mse(temp, hght, sphum):
    """Moist static energy, in Joules per kilogram."""
    from aospy.constants import L_v
    return dse(temp, hght) + L_v*sphum

def fmse(temp, hght, sphum, ice_wat):
    """Frozen moist static energy, in Joules per kilogram."""
    from aospy.constants import L_f
    return mse(temp, hght, sphum) - L_f*ice_wat

def dse(temp, hght):
    """Dry static energy, in Joules per kilogram."""
    from aospy.constants import c_p, grav
    try:
        ds = c_p*temp + grav*hght
    except ValueError:
        # On sigma coords, hght is at half levels; temp and sphum on full.
        try:
            ds = c_p*temp + grav*0.5*(hght[:,:-1] + hght[:,1:])
        except ValueError:
            ds = c_p*temp + grav*0.5*(hght[:,:,:-1] + hght[:,:,1:])
    return ds

def pot_temp(temp, p, p0=1000.):
    """Potential temperature."""
    import numpy as np
    from aospy.constants import kappa
    return temp*(p0/p[:,np.newaxis,np.newaxis])**kappa

def virt_pot_temp(temp, p, sphum, liq_wat, p0=1000.):
    """Virtual potential temperature, approximating the mixing ratios as
    specific humidities.

    """
    import numpy as np
    from aospy.constants import kappa
    return ((temp*(p0/p[:,np.newaxis,np.newaxis])**kappa) *
            (1. + 0.61*sphum - liq_wat))

def equiv_pot_temp(temp, p, sphum, p0=1000.):
    """Equivalent potential temperature."""
    import numpy as np
    from aospy.constants import kappa, L_v, c_p
    return (temp + L_v*sphum/c_p)*(p0/p[:,np.newaxis,np.newaxis])**kappa
    
### Gross moist stability-related quantities

# def dse_int_sigma(bk, pk, ps, hght, temp):
#     """Column integral of DSE using data on sigma-coordinates."""
#     import numpy as np
#     from aospy.constants import grav
#     dp_g = dp_from_sigma(bk, pk, ps) / grav
#     # temp is on level centers; hght is on level borders.
#     # So interpolate to level centers by simple averaging.
#     hght = 0.5*(hght[1:] + hght[:-1])
#     dry = dse(temp, hght)
#     return np.sum(dry*dp_g, axis=0)

def field_vert_int_max(field, dp):
    """Maximum magnitude of integral of a field from surface up."""
    import numpy as np
    from aospy.constants import grav
    dp = to_pascal(dp)
    # 2015-05-15: Problem: Sigma data indexing starts at TOA, while pressure
    #             data indexing starts at 1000 hPa.  So for now only do for
    #             sigma data and flip array direction to start from sfc.
    field_dp_g = weight_by_delta(field, dp)[::-1] / grav
    # Input array dimensions are assumed ([time dims,] level, lat, lon).
    pos_max = np.amax(np.cumsum(field_dp_g, axis=0), axis=-3)
    neg_max = np.amin(np.cumsum(field_dp_g, axis=0), axis=-3)
    # Flip sign because integrating from p_sfc up, i.e. with dp negative.
    return -1*np.where(pos_max > -neg_max, pos_max, neg_max)

def horiz_divg_vert_int_max(u, v, lat, lon, radius, dp):
    """Maximum magnitude of integral from surface up of horizontal divergence."""
    return field_vert_int_max(horiz_divg_mass_bal(u, v, lat, lon, radius, dp),
                              dp)

def vert_divg_vert_int_max(omega, p, dp):
    """Maximum magnitude of integral from surface up of vertical divergence."""
    return field_vert_int_max(vert_divg_mass_bal(omega, p, dp), dp)

def gms_like_ratio(weights, tracer, dp):
    """Compute ratio of integrals in the style of gross moist stability."""
    import numpy as np
    from aospy.constants import grav
    # Integrate weights over lower tropospheric layer
    dp = to_pascal(dp)
    denominator = field_vert_int_max(weights, dp)
    # Integrate tracer*weights over whole column and divide.
    numerator = np.sum(weight_by_delta(weights*tracer, dp), axis=-3) / grav
    return numerator / denominator

def gross_moist_strat(sphum, u, v, lat, lon, radius, dp):
    """Gross moisture stratification, in horizontal divergence form."""
    from aospy.constants import L_v
    divg = horiz_divg(u, v, lat, lon, radius)
    return L_v*gms_like_ratio(divg, sphum, dp)

def gross_dry_stab(temp, hght, u, v, lat, lon, radius, dp):
    """Gross dry stability, in horizontal divergence form."""
    divg = horiz_divg(u, v, lat, lon, radius)
    return -gms_like_ratio(divg, dse(temp, hght), dp)

def gross_moist_stab(temp, hght, sphum, u, v, lat, lon, radius, dp):
    """Gross moist stability, in horizontal divergence form."""
    divg = horiz_divg(u, v, lat, lon, radius)
    return -gms_like_ratio(divg, mse(temp, hght, sphum), dp)

def msf(lats, levs, v):
    """Meridional mass streamfunction."""
    import numpy as np
    # Compute half level boundaries and widths.
    p_top = 5.; p_bot = 1005.
    p_half = 0.5*(levs[1:] + levs[:-1])
    p_half = np.insert(np.append(p_half, p_top), 0, p_bot)
    dp = to_pascal(p_half[:-1] - p_half[1:])
    # Integrate from TOA down to surface.
    strmfunc = 2.*np.pi*r_e/grav*(np.cos(np.deg2rad(lats))[np.newaxis,:] *
                                  np.cumsum((v.mean(axis=-1) *
                                  dp[:,np.newaxis])[::-1], axis=0))[::-1]
    # Average the values calculated at half levels; flip sign by convention.
    strmfunc[:-1] = -0.5*(strmfunc[1:] + strmfunc[:-1])
    # Uppermost level goes to 0 hPa (so divide by 2); surface value is zero.
    strmfunc[-1]*=0.5; strmfunc[0] = 0.
    return strmfunc

def msf_max(lats, levs, v):
    """Maximum meridional mass streamfunction magnitude at each latitude."""
    import numpy as np
    from aospy.calcs import msf
    strmfunc = msf(lats, levs, v)
    pos_max = np.amax(strmfunc, axis=0)
    neg_max = np.amin(strmfunc, axis=0)
    return np.where(pos_max > -neg_max, pos_max, neg_max)
    #return strmfunc[:,4]

def aht(swdn_toa, swup_toa, olr, swup_sfc, swdn_sfc, lwup_sfc, lwdn_sfc,
        shflx, evap, snow_ls, snow_conv, sfc_area):
    """Total atmospheric northward energy flux."""
    import numpy as np
    from aospy.constants import L_f, L_v
    # Calculate energy balance at each grid point.
    loc = -1*(swdn_toa - swup_toa - olr +                   # TOA radiation
              swup_sfc - swdn_sfc + lwup_sfc - lwdn_sfc +   # sfc radiation
              shflx +                                       # sfc SH flux
              L_f*(snow_ls + snow_conv) + L_v*evap)         # column LH flux
    # Calculate meridional heat transport.
    glb = np.average(loc, weights=sfc_area)
    return np.cumsum(np.sum(sfc_area*(glb - loc), axis=-1))
    
def gms_up_low(temp, hght, sphum, level, lev_up=400., lev_dn=925.):
    """Gross moist stability. Upper minus lower level MSE."""
    import numpy as np
    from aospy.calcs import mse
    from aospy.constants import c_p
    m = mse(temp, hght, sphum)
    return (np.squeeze(m[np.where(level == lev_up)] -
                       m[np.where(level == lev_dn)])/c_p)

def gms_each_level(temp, hght, sphum, level, lev_dn=925.):
    import numpy as np
    from aospy.calcs import mse
    from aospy.constants import c_p
    m = mse(temp, hght, sphum)
    return (m - m[np.where(level == lev_dn)])/c_p

def dry_static_stab(temp, hght, level, lev_dn=925.):
    import numpy as np
    from aospy.calcs import dse
    from aospy.constants import c_p
    d = dse(temp, hght)
    return (d - d[np.where(level == lev_dn)])/c_p

def moist_static_stab(temp, p, sphum, p0=1000., lev_dn=925.):
    import numpy as np
    theta_e = equiv_pot_temp(temp, p, sphum, p0=p0)
    return (theta_e - theta_e[np.where(p == lev_dn)])

def gms_change_up_therm_low(temp, hght, sphum, level, lev_up=200., lev_dn=850.):
    import numpy as np
    from aospy.calcs import mse
    from aospy.constants import c_p, L_v
    """Gross moist stability. Upper minus lower level MSE with thermodynamic
    scaling estimate for low level MSE."""
    m = mse(temp, hght, sphum).mean(axis=-1)
    return (m[np.where(level == lev_up)] - m[np.where(level == lev_dn)])/c_p

def gms_h01(temp, hght, sphum, precip, level, lev_sfc=925.):
    """
    Gross moist stability. Near surface MSE diff b/w ITCZ and the given latitude.
    """
    import numpy as np
    from aospy.calcs import mse
    from aospy.constants import c_p
    # ITCZ defined as latitude with maximum zonal mean precip.
    itcz_ind = np.argmax(precip.mean(axis=-1))
    m = mse(np.squeeze(temp[np.where(level == lev_sfc)].mean(axis=-1)), 
            np.squeeze(hght[np.where(level == lev_sfc)].mean(axis=-1)), 
            np.squeeze(sphum[np.where(level == lev_sfc)].mean(axis=-1)))
    return (m[itcz_ind] - m)/c_p

def gms_h01est(temp, sphum, precip, level, lev_sfc=925.):
    """
    Gross moist stability. Near surface MSE diff b/w ITCZ and the given latitude
    neglecting the geopotential term.
    """
    import numpy as np
    from aospy.constants import c_p, L_v
    sphum = np.squeeze(sphum[np.where(level == lev_sfc)].mean(axis=-1))
    temp = np.squeeze(temp[np.where(level == lev_sfc)].mean(axis=-1))
    # ITCZ defined as latitude with maximum zonal mean precip.
    itcz_ind = np.argmax(np.mean(precip, axis=-1))
    # GMS is difference between surface values
    return temp[itcz_ind] - temp + L_v*(sphum[itcz_ind] - sphum)/c_p

def gms_h01est2(temp, hght, sphum, precip, level, lev_up=200., lev_sfc=925.):
    """
    Gross moist stability. MSE diff b/w ITCZ aloft and near surface at the 
    given latitude.
    """
    import numpy as np
    from aospy.calcs import mse
    from aospy.constants import c_p
    # ITCZ defined as latitude with maximum zonal mean precip.
    itcz_ind = np.argmax(precip.mean(axis=-1))
    m_up = mse(np.squeeze(temp[np.where(level == lev_up)].mean(axis=-1)), 
            np.squeeze(hght[np.where(level == lev_up)].mean(axis=-1)), 
            np.squeeze(sphum[np.where(level == lev_up)].mean(axis=-1)))
    m_sfc = mse(np.squeeze(temp[np.where(level == lev_sfc)].mean(axis=-1)), 
            np.squeeze(hght[np.where(level == lev_sfc)].mean(axis=-1)), 
            np.squeeze(sphum[np.where(level == lev_sfc)].mean(axis=-1)))
    return (m_up[itcz_ind] - m_sfc)/c_p

def gms_change_est(T_cont, T_pert, q_cont, precip, level, lev_sfc=925.):
    """
    Gross moist stability change estimated as near surface MSE difference
    between ITCZ and local latitude, neglecting geopotential term and applying
    a thermodynamic scaling for the moisture term.
    """
    import numpy as np
    from aospy.constants import c_p, L_v
    # ITCZ defined as latitude with maximum zonal mean precip.
    itcz_ind = np.argmax(precip.mean(axis=-1))
    # Need temperature change at 
    T_pert = np.squeeze(T_pert[np.where(level == lev_sfc)].mean(axis=-1))
    T_cont = np.squeeze(T_cont[np.where(level == lev_sfc)].mean(axis=-1))
    dT = T_pert - T_cont
    dT_itcz = T_pert[itcz_ind] - T_cont[itcz_ind]
    q_cont = np.squeeze(q_cont[np.where(level == lev_sfc)].mean(axis=-1))
    # GMS is difference between surface
    alpha = 0.07
    return ((c_p + L_v*alpha*q_cont[itcz_ind])*dT_itcz - 
            (c_p + L_v*alpha*q_cont)*dT)/c_p

def gms_change_est2(T_cont, T_pert, q_cont, precip, level, lat,
                    lev_sfc=925., gamma=1.):
    """
    Gross moist stability change estimated as near surface MSE difference
    between ITCZ and local latitude, neglecting geopotential term and applying
    a thermodynamic scaling for the moisture term, and multiplying the ITCZ
    terms by cos(lat) and a fixed fraction gamma to account for deviation of
    upper level MSE from the near surface ITCZ value.
    """
    import numpy as np
    from aospy.constants import c_p, L_v
    # ITCZ defined as latitude with maximum zonal mean precip.
    itcz_ind = np.argmax(precip.mean(axis=-1))
    # Need temperature change at 
    T_pert = np.squeeze(T_pert[np.where(level == lev_sfc)].mean(axis=-1))
    T_cont = np.squeeze(T_cont[np.where(level == lev_sfc)].mean(axis=-1))
    dT = T_pert - T_cont
    dT_itcz = T_pert[itcz_ind] - T_cont[itcz_ind]
    q_cont = np.squeeze(q_cont[np.where(level == lev_sfc)].mean(axis=-1))
    # GMS is difference between surface
    alpha = 0.07
    return (np.cos(np.deg2rad(lat))**2*gamma*
            (c_p + L_v*alpha*q_cont[itcz_ind])*dT_itcz - 
            (c_p + L_v*alpha*q_cont)*dT)/c_p

def prec_conv_frac(prec_conv, precip, prec_ls=False):
    """Fraction of precipitation coming from convection scheme."""
    # Mask where precip is zero to avoid dividing by zero.
    from numpy.ma import masked_where
    prec_conv = masked_where(precip == 0., prec_conv)
    precip = masked_where(precip == 0., precip)
    if prec_ls:
        return prec_conv/(precip + prec_conv)
    else:
        return prec_conv/precip

def descent_tot(omega, mc):
    """Vertical motion from both convection and large-scale."""
    from aospy.constants import grav
    return omega + grav*mc

def ascent_ls(omega):
    """Large-scale vertically upward motion."""
    # Get positive values and replace negative ones with zero.
    from numpy import where
    return where(omega > 0., omega, 0)

def vert_centroid(field, level, p_bot=850., p_top=150.):
    """
    Compute the vertical centroid of some vertically defined field.
    """
    import numpy as np
    from aospy import calc_levs_thick
    desired_levs = np.where((level <= p_bot) & (level >= p_top))
    lev_thick = calc_levs_thick(level)/100.
    # Add axes for later broadcasting and truncate to desired vertical levels.
    level = level[desired_levs]; level = level[:,np.newaxis,np.newaxis]
    lev_thick = lev_thick[desired_levs]
    lev_thick = lev_thick[:,np.newaxis,np.newaxis]
    field = field[desired_levs]
    # For 1D arrays, have to move the vertical coordinate to leftmost dim.
    if field.ndim == 1:
        field = np.atleast_3d(field).swapaxes(0,1)
    else:
        field = np.atleast_3d(field)
    return (np.sum(field*level*lev_thick, axis=0) /
            np.sum(field*lev_thick, axis=0))

###############################
# Functions below this line haven't been converted to new argument format.

def tht(variables, **kwargs):
    """Total atmospheric plus oceanic northward energy flux."""
    import numpy as np
    
    # Calculate energy balance at each grid point.
    loc = -1*(variables[0] - variables[1] - variables[2])             # TOA radiation
    # Calculate meridional heat transport.
    len_dt = variables[0].shape[0]
    sfc_area = grid_sfc_area(nc)
    glb = np.average(loc.reshape(len_dt, -1), 
                     weights=sfc_area.ravel(), axis=1)
    # AHT is meridionally integrated energy flux divergence.
    flux_div = np.sum(sfc_area*(glb[:,np.newaxis,np.newaxis] - loc), axis=-1)
    return np.cumsum(flux_div, axis=-1)
    
def oht(variables, **kwargs):
    """Total oceanic northward energy flux as residual of total minus atmospheric flux."""
    import numpy as np
    from aospy.constants import L_f, L_v
    
    # Calculate energy balance at each grid point.
    loc = (variables[0] - variables[1] + variables[2] - variables[3] +   # sfc radiation
           variables[4] +                                 # sfc SH flux
           L_f*(variables[5] + variables[6]) + L_v*variables[7])   # column LH flux
    # Calculate meridional heat transport.
    len_dt = variables[0].shape[0]
    sfc_area = grid_sfc_area(nc)
    glb = np.average(loc.reshape(len_dt, -1), 
                     weights=sfc_area.ravel(), axis=1)
    # AHT is meridionally integrated energy flux divergence.
    flux_div = np.sum(sfc_area*(glb[:,np.newaxis,np.newaxis] - loc), axis=-1)
    return np.cumsum(flux_div, axis=-1)
    #return tht(variables, **kwargs) - aht(variables, **kwargs)

def moc_flux(variables, **kwargs):
    """Mass weighted column integrated meridional flux by time and 
    zonal mean flow."""
    import numpy as np
    from aospy.av_stat import levs_thick
    from aospy.calcs import dse, mse
    # Specify upper bound of vertical integrals.
    p_top = kwargs.get('p_top', 0.)
    trop = np.where(nc.variables['level'][:] >= p_top)
    # Apply mass flux correction to zonal mean data.
    v_znl = np.squeeze(variables[-1][:,trop]).mean(axis=-1)
    v_north = np.where(v_znl > 0., v_znl, 0.)
    v_south = np.where(v_znl < 0., v_znl, 0.)
    lev_thick = np.squeeze(levs_thick(nc)[:,trop])/grav
    lev_thick = lev_thick[np.newaxis,:,np.newaxis]
    # Adjustment imposes that column integrated mass flux is zero.
    mass_adj = -((v_north*lev_thick).sum(axis=1) / 
                 (v_south*lev_thick).sum(axis=1))
    # Integrate the specified flux by the adjusted v vertically and zonally.
    flux_type = kwargs['flux_type']
    if flux_type == 'dse':
        flux = (np.squeeze(dse(variables[:2])[:,trop]).mean(axis=-1) *
                (v_north + v_south * mass_adj[:,np.newaxis,:]))
    elif flux_type == 'mse':
        flux = (np.squeeze(mse(variables[:3])[:,trop]).mean(axis=-1) *
                (v_north + v_south * mass_adj[:,np.newaxis,:]))
    elif flux_type == 'moisture':
        flux = L_v*(np.squeeze(variables[0][:,trop]).mean(axis=-1) *
                (v_north + v_south * mass_adj[:,np.newaxis,:]))
    return (2.*np.pi*r_e*np.cos(np.deg2rad(nc.variables['lat'][:])) * 
            (flux*lev_thick).sum(axis=1))
    
def moc_flux_raw(variables, **kwargs):
    """Mass weighted column integrated meridional flux by time and zonal mean flow, without applying column mass flux correction."""
    from aospy.av_stat import levs_thick
    # Specify upper bound of vertical integrals.
    p_top = kwargs.get('p_top', 0.)
    trop = np.where(nc.variables['level'][:] >= p_top)
    # Take zonal mean and calculate grid level thicknesses.
    v_znl = np.squeeze(variables[-1][:,trop]).mean(axis=-1)
    lev_thick = np.squeeze(levs_thick(nc)[:,trop])/grav
    lev_thick = lev_thick[np.newaxis,:,np.newaxis]
    # Integrate the specified flux vertically and zonally.
    flux_type = kwargs.get('flux_type', 'mse')
    if flux_type == 'dse':
        flux = np.squeeze(dse(variables[:2])[:,trop]).mean(axis=-1) * v_znl
    elif flux_type == 'mse':
        flux = np.squeeze(mse(variables[:3])[:,trop]).mean(axis=-1) * v_znl
    elif flux_type == 'moisture':
        flux = L_v*np.squeeze(variables[0][:,trop]).mean(axis=-1) * v_znl
    return (2.*np.pi*r_e*np.cos(np.deg2rad(nc.variables['lat'][:])) * 
            (flux*lev_thick).sum(axis=1))

def st_eddy_flux(variables, **kwargs):
    """Mass weighted column integrated meridional flux by stationary eddies."""
    import numpy as np
    from aospy.av_stat import levs_thick    
    p_top = kwargs.get('p_top', 0.)
    trop = np.where(nc.variables['level'][:] >= p_top)
    v = np.squeeze(variables[-1][:,trop])
    flux_type = kwargs['flux_type']
    if flux_type == 'dse':
        m = np.squeeze(dse(variables[:2])[:,trop])
    elif flux_type == 'mse':
        m = np.squeeze(mse(variables[:3])[:,trop])
    elif flux_type == 'moisture':
        m = np.squeeze(variables[0][:,trop])*L_v
    lev_thick = np.squeeze(levs_thick(nc)[:,trop])/grav
    lev_thick = lev_thick[np.newaxis,:,np.newaxis,np.newaxis]
    flux = ((m - m.mean(axis=-1)[:,:,:,np.newaxis]) * 
            (v - v.mean(axis=-1)[:,:,:,np.newaxis]))
    return (2.*np.pi*r_e*np.cos(np.deg2rad(nc.variables['lat'][:])) * 
            (flux*lev_thick).sum(axis=1).mean(axis=-1))
            
def moc_st_eddy_flux(variables, **kwargs):
    """Mass weighted column integrated flux by time mean flow."""
    from aospy.calcs import moc_flux, st_eddy_flux
    return moc_flux(variables, **kwargs) + st_eddy_flux(variables, **kwargs)

def trans_eddy_flux(variables, **kwargs):
    """Meridional flux by transient eddies."""
    from aospy.calcs import aht, moc_st_eddy_flux
    return aht(variables[4:], **kwargs) - moc_st_eddy_flux(variables[:4], **kwargs)

def eddy_flux(variables, **kwargs):
    """Meridional flux by stationary and transient eddies."""
    from aospy.calcs import aht, moc_flux
    return aht(variables[4:], **kwargs) - moc_flux(variables[:4], **kwargs)

def mse_flux(variables, **kwargs):
    """Column integrated moist static energy meridional flux."""
    flux_type = kwargs.get('flux_type', 'moc')
    if flux_type == 'moc':
        from aospy.calcs import moc_flux
        return moc_flux(variables[:4], **kwargs)
    elif flux_type == 'st_eddy':
        from aospy.calcs import st_eddy_flux
        return st_eddy_flux(variables[:4], **kwargs)
    elif flux_type == 'moc_st_eddy':
        from aospy.calcs import moc_st_eddy_flux
        return moc_st_eddy_flux(variables[:4], **kwargs)
    elif flux_type =='trans_eddy':
        from aospy.calcs import trans_eddy_flux
        return trans_eddy_fux(variables, **kwargs)
    elif flux_type == 'all':
        from aospy.calcs import aht
        return aht(variables[4:], **kwargs)
    elif flux_type == 'eddy':
        from aospy.calcs import eddy_flux
        return eddy_flux(variables, **kwargs)
    
def mass_flux(variables, **kwargs):
    """Meridional mass flux by time and zonal mean flow."""
    import numpy as np
    # Apply mass flux correction.
    lev_thick = levs_thick(nc)[np.newaxis,:,np.newaxis]/grav
    v_znl = variables[0].mean(axis=-1)
    v_north = np.where(v_znl > 0., v_znl, 0.)
    v_south = np.where(v_znl < 0., v_znl, 0.)
    mass_adj = -((v_north*lev_thick).sum(axis=1) / 
                 (v_south*lev_thick).sum(axis=1))
    # Integrate vertically and pick level where magnitude maximized.
    int_flux = ((v_north + v_south*mass_adj[:,np.newaxis,:]) * 
                lev_thick).cumsum(axis=1)
    flux_pos = np.amax(int_flux, axis=1)
    flux_neg = np.amin(int_flux, axis=1) 
    return (2.*np.pi*r_e*np.cos(np.deg2rad(nc.variables['lat'][:])) * 
            np.where(flux_pos > - flux_neg, flux_pos, flux_neg))

def gms_moc(variables, **kwargs):
    """Gross moist stability."""
    from aospy.calcs import moc_flux, msf_max
    from aospy.constants import c_p
    return -moc_flux(variables, **kwargs)/msf_max([variables[-1]], **kwargs)/c_p

def gms_msf(variables, **kwargs):
    """Gross moist stability."""
    from aospy.calcs import moc_st_eddy_flux, msf_max
    from aospy.constants import c_p
    return -(moc_st_eddy_flux(variables, **kwargs) / 
            (msf_max([variables[-1]], **kwargs)*c_p))

def total_gms(variables, **kwargs):
    """Total (mean plus eddy) gross moist stability."""
    from aospy.calcs import aht, msf_max
    from aospy.constants import c_p
    return -(aht(variables[:-1], **kwargs) /
             msf_max([variables[-1]], **kwargs))/c_p

def aht_no_snow(variables, **kwargs):
    """Total atmospheric northward energy flux."""
    import numpy as np
    from aospy.constants import L_v
    
    # Calculate energy balance at each grid point.
    loc = -1*(variables[0] - variables[1] - variables[2] +             # TOA radiation
              variables[3] - variables[4] + variables[5] - variables[6] +   # sfc radiation
              variables[7] +                                 # sfc SH flux
              L_v*variables[-1])   # column LH flux
    # Calculate meridional heat transport.
    len_dt = variables[0].shape[0]
    sfc_area = grid_sfc_area(nc)
    glb = np.average(loc.reshape(len_dt, -1), 
                     weights=sfc_area.ravel(), axis=1)
    # AHT is meridionally integrated energy flux divergence.
    flux_div = np.sum(sfc_area*(glb[:,np.newaxis,np.newaxis] - loc), axis=-1)
    return np.cumsum(flux_div, axis=-1)
