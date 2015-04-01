def weighted_variance(vals, weights):
    """Calculate weighted spatial variance of a time-series array."""
    import numpy as np
    from numpy import ma
    # Assume computing at single timestep.
    try:
        return (ma.sum(weights*(vals - ma.average(vals, weights=weights))**2)
                / weights.sum())
    # If not, compute for whole timeseries.
    except:
        return (ma.sum(weights*(vals - ma.average(vals, weights=weights, 
                axis=1)[:,np.newaxis])**2, axis=1) / weights.sum())

def region_avg(vals, region, model, data_type='ts', is_znl=False):
    """Calculate time average and variance over desired regions."""
    import numpy as np
    from numpy import ma
    # Interpolate region to model grid.
    reg_mask = region.make_mask(model)
    # Mask input values where region mask is zero.
    if data_type == 'ts':
        try:
            vals = ma.masked_where(np.tile(reg_mask == 0., 
                                   (vals.shape[0], 1, 1)), vals)
        except:
            vals = ma.masked_where(np.tile(reg_mask == 0., (vals.shape[0],
                                   vals.shape[1], 1, 1)), vals)
    elif data_type == 'av':
        try:
            vals = ma.masked_where(reg_mask == 0., vals)
        except:
            vals = ma.masked_where(np.tile(reg_mask == 0., 
                                   (vals.shape[0], 1, 1)), vals)
    # Where nonzero, weight by surface area times region mask value.
    weights = ma.masked_where(reg_mask == 0, model.sfc_area*reg_mask)
    # Time average at each point is used for computing spatial variance.
    if data_type == 'ts':
        loc_av = vals.mean(axis=0)
    elif data_type == 'av':
        loc_av = vals
    # Syntax differs if data is defined meridionally or not (e.g. zonal mean).
    if is_znl:
        vals = vals.mean(axis=-1)
        weights = weights.mean(axis=-1)
        if data_type == 'av':
            return (weighted_variance(vals, weights),
                    weighted_variance(loc_av.mean(axis=-1), weights))
        else:
            return weighted_variance(vals, weights)
    else:
        # Syntax differs if time-series data is defined vertically or not.
        if data_type == 'ts':
            if vals.ndim == 3:
                vals = vals.reshape(vals.shape[0], -1)
                weights = np.tile(weights.ravel(), (vals.shape[0], 1))
            elif vals.ndim == 4:
                vals = vals.reshape(vals.shape[0], vals.shape[1], -1)
                weights = np.tile(weights.ravel(), 
                                  (vals.shape[0], vals.shape[1], 1))
        # Return the region average time-series, mean, standard deviation,
        # and spatial variance time-series and mean.
        if data_type == 'ts':
            reg_ts = ma.average(vals, weights=weights, axis=-1)
            return (reg_ts, reg_ts.mean(axis=0), reg_ts.std(axis=0),
                    weighted_variance(vals, weights),
                    weighted_variance(loc_av.ravel(), weights[0]))
        elif data_type == 'av':
            try:
                av = ma.average(vals, weights=weights)
            except:
                av = ma.average(vals.reshape(vals.shape[0], -1), 
                                weights=weights.ravel(), axis=1)
            try:
                wt_var = weighted_variance(vals, weights)
            except:
                wt_var = weighted_variance(vals.reshape(vals.shape[0], -1),
                                           weights=weights.ravel())
            return av, wt_var

def region_vals(vals, level, proj, model, data_type, is_znl):
    """
    Create alphabetized dicts for each region calculation linking the region
    name to its value for the given calculation.
    """
    import numpy as np
    names = [reg.name for reg in proj.regions]
    # Make alphabetized dict of regions and their values.
    alph_dict = lambda n, v: dict(sorted(zip(n, v)))
    data = [region_avg(vals, reg, model, is_znl=is_znl, data_type=data_type)
            for reg in proj.regions]
    dicts = [alph_dict(names, dat) for dat in np.transpose(data)]
    # Save differences between certain regions.
    for dictionary in dicts:
        if 'nh' in names and 'sh' in names:
            dictionary.update({'nh-sh': dictionary['nh'] - dictionary['sh']})
        if 'epac' in names and 'wpwp' in names:
            dictionary.update({'epac-wpwp':
                               dictionary['epac'] - dictionary['wpwp']})
        if 'burls_epac' in names and 'burls_wpac' in names:
            dictionary.update({'burls_wpac-epac': dictionary['burls_wpac'] - 
                               dictionary['burls_epac']})
        if 'burls_trop_pac' in names and 'burls_ext_pac' in names:
            dictionary.update({'burls_merid_pac': 
                               dictionary['burls_ext_pac'] - 
                               dictionary['burls_trop_pac']})
    return dicts
