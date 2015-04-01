def hadley_bounds(lats, levs, vcomp):
    """Poleward extent of Hadley Cell."""
    import numpy as np
    from aospy.calcs import msf
    # Get meridional mass streamfunction at 500 hPa.
    p_ind = np.where(levs == 500.)
    sf = np.squeeze(msf(lats, levs, vcomp)[p_ind])
    zero_cross = np.where(np.diff(np.sign(sf)))[0]
    # Hadley Cell center is the zero crossing nearest the equator.
    ind_cent = np.argmin(np.abs(zero_cross - lats.size/2))
    # Hadley cell poleward boundaries are the zero crossings on either side.
    ind = [zero_cross[ind_cent + i] for i in range(-1,2)]
    # Linearly interpolate to true zero crossings.
    return np.array([lats[ind[i]+1] - (lats[ind[i]+1] - lats[ind[i]]) / 
                (sf[ind[i]+1] - sf[ind[i]])*sf[ind[i]+1] 
                for i in range(3)])

def had_bounds(strmfunc, lat, return_max=False):
    """Hadley cell poleward extent and center location."""
    import numpy as np
    # Get data max and min values and indices, such that min is north of max.
    z_max_ind, y_max_ind = np.where(strmfunc == strmfunc.max())
    z_max_ind = z_max_ind[0]; y_max_ind=y_max_ind[0]
    z_min_ind, y_min_ind = np.where(strmfunc[:,y_max_ind:] ==
                                    strmfunc[:,y_max_ind:].min())
    z_min_ind = z_min_ind[0]; y_min_ind = y_min_ind[0] + y_max_ind
    min_north = strmfunc[z_min_ind, y_min_ind]
    # Return locations and values of Hadley cells' strengths.
    if return_max:
        return [z_min_ind, y_min_ind, strmfunc[z_min_ind, y_min_ind],
                z_max_ind, y_max_ind, strmfunc[z_max_ind, y_max_ind]]
    # Find latitude where streamfunction at its level of maximum decreases
    # to 10% of that maximum value. 
    had_max = np.where(np.diff(np.sign(strmfunc[z_max_ind] - 
                                       0.1*strmfunc.max())))[0]
    had_min = np.where(np.diff(np.sign(strmfunc[z_min_ind] -
                                       0.1*strmfunc.min())))[0]
    had_lims = (np.intersect1d(range(y_max_ind), had_max)[-1],
               np.intersect1d(range(y_min_ind, lat.size), had_min)[0])
    # Center is streamfunction zero crossing between the two cells at level
    # halfway between the levels of the two maxima.
    p_ind = 0.5*(z_min_ind + z_max_ind)
    zero_cross = np.where(np.diff(np.sign(np.squeeze(
        strmfunc[p_ind,y_max_ind:y_min_ind]))))[0]
    had_center = zero_cross[0] + y_max_ind
    return np.array([lat[had_lims[0]], lat[had_center], lat[had_lims[1]]])

def had_bounds500(strmfunc, lat):
    """Hadley cells extent based on 500 hPa streamfunction zero crossings."""
    import numpy as np
    # Find latitudes where streamfunction at 500 hPa changes sign.
    strmfunc = strmfunc[5]
    zero_cross = np.where(np.diff(np.sign(strmfunc)))[0]
    # Hadley Cell center is the zero crossing nearest the equator.
    ind_cent = np.argmin(np.abs(zero_cross - lat.size/2))
    # Hadley cell poleward boundaries are the zero crossings on either side.
    ind = [zero_cross[ind_cent + i] for i in range(-1,2)]
    # Linearly interpolate to true zero crossings.
    return np.array([lat[ind[i]+1] - (lat[ind[i]+1] - lat[ind[i]]) / 
                    (strmfunc[ind[i]+1] - strmfunc[ind[i]])*strmfunc[ind[i]+1] 
                    for i in range(3)])

def thermal_equator(flux, lat):
    """Location of zero-crossing of energy flux."""
    import numpy as np
    # Find latitude indices where flux changes sign.
    zero_cross = np.where(np.diff(np.sign(flux)))[0]
    # Thermal equator is the zero crossing nearest the equator.
    ind = zero_cross[np.argmin(np.abs(zero_cross - lat.size/2))]
    # Linearly interpolate to true zero crossing latitude.
    return (lat[ind+1] - (lat[ind+1] - lat[ind]) / 
            (flux[ind+1] - flux[ind])*flux[ind+1])

def itcz_pos(precip, lat, return_indices=False):
    """Calculate ITCZ location."""
    import numpy as np
    # Find the latitude index with maximum precip and calculate dP/d(latitude)
    # for the adjacent grid latitudes.
    ind_max = precip.argmax()
    if return_indices:
        return ind_max
    dP = np.gradient(precip[ind_max-1:ind_max+2])
    phi = np.deg2rad(lat)[ind_max-1:ind_max+2]
    dphi = np.gradient(phi)
    dP_dphi = dP/dphi
    # Find on which side of the maximum does dP/d(latitude) change sign.
    i = np.where(np.diff(np.sign(dP)))[0][0]
    # Linearly interpolate to determine the latitude of max precip.
    zero_interp = phi[i] - (dP_dphi[i]*(phi[i+1] - phi[i]) / 
                            (dP_dphi[i+1] - dP_dphi[i]))
    return np.rad2deg(zero_interp)

def itcz_loc(lats, precip):
    """Calculate ITCZ location."""
    import numpy as np
    # Find the latitude index with maximum precip and calculate dP/d(latitude)
    # for the adjacent grid latitudes.
    ind_max = precip.argmax()
    dP = np.gradient(precip[ind_max-1:ind_max+2])
    phi = np.deg2rad(lats)[ind_max-1:ind_max+2]
    dphi = np.gradient(phi)
    dP_dphi = dP/dphi
    # Find on which side of the maximum does dP/d(latitude) change sign.
    i = np.where(np.diff(np.sign(dP)))[0][0]
    # Linearly interpolate to determine the latitude of max precip.
    return np.deg2rad(phi[i] - (dP_dphi[i]*(phi[i+1] - phi[i]) / 
                               (dP_dphi[i+1] - dP_dphi[i])))

def prec_centroid(precip, lat, lat_max=20.):
    """
    Calculate ITCZ location as the centroid of the area weighted zonal-mean P.
    """
    import numpy as np
    # Interpolate zonal mean precip to a 0.1 degree latitude grid for 20S-20N
    trop = np.where(np.abs(lat) < lat_max)
    lat_interp = np.arange(-lat_max, lat_max + 0.01, 0.1)
    prec = np.interp(lat_interp, lat[trop], precip[trop])
    # Integrate area-weighted precip and find median.
    prec_int = np.cumsum(prec*np.abs(np.cos(np.deg2rad(lat_interp))))
    return lat_interp[np.argmin(np.abs(prec_int - 0.5*prec_int[-1]))]

def precip_centroid(lats, precip, lat_max=20.):
    """
    Calculate ITCZ location as the centroid of the area weighted zonal-mean P.
    """
    import numpy as np
    # Interpolate zonal mean P to 0.1 degree latitude grid over desired extent. 
    trop = np.where(np.abs(lats) < lat_max)
    precip = precip.mean(axis=-1)[trop]
    lat_interp = np.arange(-lat_max, lat_max + 0.01, 0.1)
    prec = np.interp(lat_interp, lats[trop], precip)
    # Integrate area-weighted precip and find median.
    prec_int = np.cumsum(prec*np.abs(np.cos(np.deg2rad(lat_interp))), axis=0)
    return lat_interp[np.argmin(np.abs(prec_int - 0.5*prec_int[-1]), axis=0)]

def ang_mom(lats, ucomp):
    """Angular momentum per unit mass."""
    import numpy as np
    from aospy.constants import r_e, Omega
    cos_lat = np.cos(np.deg2rad(lats[np.newaxis,:,np.newaxis]))
    return (Omega*r_e*cos_lat + ucomp)*r_e*cos_lat
