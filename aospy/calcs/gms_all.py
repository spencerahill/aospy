def _integrate_up_at_lat(integrand, lat, dz):
    """SAH 2014-11-06.  Attempt to generalize vertical integral at a given
    latitude, which is used by all meridional flux computations."""
    import numpy as np
    from aospy.constants import r_e
    return 2.*np.pi*r_e * np.sum(np.cos(lat)*integrand*dz, axis=1)

def mse(temp, hght, sphum):
    """Moist static energy, in Joules per kilogram."""
    from aospy.constants import c_p, grav, L_v
    return c_p*temp + grav*hght + L_v*sphum

def fmse(temp, hght, sphum, ice_wat):
    """Frozen moist static energy, in Joules per kilogram."""
    from aospy.constants import c_p, grav, L_v, L_f
    return c_p*temp + grav*hght + L_v*sphum - L_f*ice_wat

def dse(temp, hght):
    """Dry static energy, in Joules per kilogram."""
    from aospy.constants import c_p, grav
    return c_p*temp + grav*hght

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

def int_dp_g(integrand, level):
    """Weight the integrand by dp/g for use in mass-weighted integrals."""
    import numpy as np
    from aospy import calc_levs_thick
    from aospy.constants import grav
    dp_g = calc_levs_thick(level)[:,np.newaxis,np.newaxis]/grav
    return integrand*dp_g

def divg_int_max(level, divg):
    """Maximum magnitude of divergence integral from surface up."""
    import numpy as np
    divg_dp_g = int_dp_g(divg, level)
    # Input array dimensions are assumed (level, lat, lon).
    pos_max = np.amax(np.cumsum(divg_dp_g, axis=0), axis=0)
    neg_max = np.amin(np.cumsum(divg_dp_g, axis=0), axis=0)
    return np.where(pos_max > -neg_max, pos_max, neg_max)

def gms_like_ratio(level, divg, tracer):
    """
    Compute ratio of integrals in the style of gross moist stability and
    related quantities.
    """
    import numpy as np
    # Integrate divergence over lower tropospheric layer
    div_v2 = divg_int_max(level, divg)
    # Integrate tracer*divergence over whole column and divide.
    tracerdiv = np.sum(int_dp_g(divg*tracer, level), axis=0)
    return tracerdiv / div_v2

def moisture_strat(level, divg, sphum):
    """Gross moisture stratification, in horizontal divergence form."""
    from aospy.constants import L_v
    return L_v*gms_like_ratio(level, divg, sphum)

def gross_dry_stab(level, divg, temp, hght):
    """Gross dry stability, in horizontal divergence form."""
    return -gms_like_ratio(level, divg, dse(temp, hght))

def gross_moist_stab(level, divg, temp, hght, sphum):
    """Gross moist stability, in horizontal divergence form."""
    return -gms_like_ratio(level, divg, mse(temp, hght, sphum))

def divergence(omega, level, ps):
    """Horizontal divergence computed using omega & continuity equation."""
    import numpy as np
    level = np.tile(level[:,np.newaxis,np.newaxis], ps.shape)
    lev_bot = np.amax(np.where(level <= ps, level, ps))
    k_bot = np.argmax(level <= ps, axis=0)
    domega_dp = np.empty(omega.shape)
    print level.shape, k_bot.shape, domega_dp.shape, omega[k_bot].shape
    # One-sided differencing from ps (where omega=0) to the next level above.
    domega_dp[0] = omega[k_bot] / (level[k_bot] - ps)
    # One-sided differencing from topmost p to p=0, where omega=0.
    domega_dp[-1] = omega[-1] / level[-1]
    domega_dp[1:-1] = (omega[2:-1] - omega[1:-2]) / (level[2:-1] - level[1:-2])
    return np.ma.array(domega_dp, mask=omega.mask)
    
# def divergence(ucomp, vcomp):
#     """Horizontal divergence computed using the windspharm package."""
#     import numpy as np
#     from windspharm.standard import VectorWind
#     w = VectorWind(ucomp, vcomp)
#     return w.divergence()

def msf(lats, levs, vcomp):
    """Meridional mass streamfunction."""
    import numpy as np
    # Compute half level boundaries and widths.
    p_top = 5.; p_bot = 1005.
    p_half = 0.5*(levs[1:] + levs[:-1])
    p_half = np.insert(np.append(p_half, p_top), 0, p_bot)
    dp = 100.*(p_half[:-1] - p_half[1:]) # convert from hPa to Pa
    # Integrate from TOA down to surface.
    strmfunc = 2.*np.pi*r_e/grav*(np.cos(np.deg2rad(lats))[np.newaxis,:] *
                                  np.cumsum((vcomp.mean(axis=-1) *
                                  dp[:,np.newaxis])[::-1], axis=0))[::-1]
    # Average the values calculated at half levels; flip sign by convention.
    strmfunc[:-1] = -0.5*(strmfunc[1:] + strmfunc[:-1])
    # Uppermost level goes to 0 hPa (so divide by 2); surface value is zero.
    strmfunc[-1]*=0.5; strmfunc[0] = 0.
    return strmfunc

def msf_max(lats, levs, vcomp):
    """Maximum meridional mass streamfunction magnitude at each latitude."""
    import numpy as np
    from aospy.calcs import msf
    strmfunc = msf(lats, levs, vcomp)
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

def qu(sphum, ucomp):
    """"Zonal moisture flux."""
    return sphum*ucomp

def qv(sphum, vcomp):
    """Meridional moisture flux."""
    return sphum*vcomp

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
