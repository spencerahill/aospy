def prep_plot_data(plot_type, proj, model_name, run_name, region, ens_mem,
                   var_name, level, intvl, yr_range, **kwargs):
    """Load and prep data for plotting."""
    import numpy as np
    from scipy.stats import t
    from aospy.io import (load_plot_data, _proj_inst, _model_inst, 
                          _run_inst, _var_inst)
    # Create aospy objects.
    proj = _proj_inst(proj)
    model = _model_inst(model_name, proj)
    run = _run_inst(run_name, model, proj)
    var = _var_inst(var_name)
    # Load desired data.
    if plot_type == 'time_vert':
        intvl = 'ann_cycle'
        dtype = 'reg.av'
        dtype_std = 'reg.std'
    elif plot_type == 'lon_lat':
        dtype = 'av'
        dtype_std = 'std'
        region = False
    elif plot_type == 'lon_vert':
        dtype = 'av'
        dtype_std = 'std'
    elif plot_type == 'time_1d':
        if kwargs.get('xlim', False) == 'ann_cycle':
            intvl = 'ann_cycle'
            dtype = 'reg.av'
            dtype_std = 'reg.std'
        else:
            intvl = intvl
            dtype = 'reg.ts'
            dtype_std = 'reg.std'
    elif plot_type == '1d_vert':
        dtype = 'reg.av'
        dtype_std = 'reg.std'
    # Grab only lats/lons to be plotted.
    map_corners = kwargs.get('map_corners', False)
    if map_corners:
        lats = (map_corners['llcrnrlat'], map_corners['urcrnrlat'])
        lons = (map_corners['llcrnrlon'], map_corners['urcrnrlon'])
    else:
        lats, lons = False, False
    data = load_plot_data(proj, model, run, ens_mem, var, level, 
                          intvl, dtype, yr_range, do_znl_mean=False,
                          region=region, lats=lats, lons=lons)
    # Apply desired transforms.
    if kwargs.get('do_zasym', False):
        data = data - data.mean(axis=1)[:,np.newaxis]
    if kwargs.get('do_subtract_mean', False):
        data = data - ma.average(data, weights=model.sfc_area)
    print "%s max and min: %0.2f %0.2f " % (run, np.max(data), np.min(data))
    # Get year range.
    if yr_range == 'default':
        try:
            start_yr, end_yr = run.default_yr_range
        except AttributeError:
            start_yr, end_yr = model.default_yr_range
    else:
        start_yr, end_yr = yr_range
    num_yr = end_yr - start_yr + 1
    # Load control run data, if desired.
    cont_run = kwargs.get('cont_run', 'cont')
    do_cont_cntrs = kwargs.get('do_cont_cntrs', False)
    if do_cont_cntrs:
        data_cont = load_plot_data(
            proj, model, cont_run, ens_mem, var, level, intvl, dtype, 
            yr_range, do_znl_mean=False, region=region
        )
    else:
        data_cont = False    
    # Load statistical masking/hatching/shading data, if desired.
    do_stat_mask = kwargs.get('do_stat_mask', False)
    do_plot_shading = kwargs.get('do_plot_shading', False)
    if do_stat_mask or do_plot_shading:
        data_std = load_plot_data(
            proj, model, run, ens_mem, var, level, intvl, dtype_std, 
            yr_range, do_znl_mean=False, region=region
        )
        # Either use stdev directly or compute t-test.
        if do_plot_shading or do_stat_mask in ['std', 'stdev']:
            data_stat = data_std
        elif do_stat_mask == 'ttest':
            data_stat = t.pdf(data/data_std, num_yr - 1)
        # Mask data if desired.
        stat_mask_hatch = kwargs.get('stat_mask_hatch', False) 
        if stat_mask_hatch == 'mask':
            ax_fill_color = '0.85'
            data = np.ma.masked_where(data_stat > .05, data)
        else:
            ax_fill_color = None
    else:
        ax_fill_color = None
        data_stat = False
    # # If pert minus control, make `run` the control.
    # try:
    #     run = run[-1]
    # except TypeError:
    #     pass
    # Normalize by specified value, if desired.
    do_norm_cont = kwargs.get('do_norm_cont', False)
    do_normalize = kwargs.get('do_normalize', False)
    if do_normalize or do_norm_cont:
        norm_region = kwargs.get('norm_region', region)
        if norm_region:
            dtype = 'reg.av'
        if do_norm_cont:
            norm_run = cont_run
            norm_var = var
        else:
            norm_run = kwargs.get('norm_run', cont_run)
            if not norm_run:
                norm_run = cont_run
            norm_var = kwargs.get('norm_var', var)
            if not norm_var:
                norm_var = var
        data_norm = load_plot_data(
            proj, model, norm_run, ens_mem, norm_var, level,
            intvl, dtype, yr_range, do_znl_mean=False,
            region=norm_region
        )
        try:
            data /= data_norm
        except ValueError:
            data /= data_norm[:,np.newaxis]
    return (proj, model, run, var, data, data_cont, data_stat,
            ax_fill_color, start_yr, end_yr)

def _set_axes_props(ax, var, **kwargs):
    from matplotlib import pyplot as plt
    """Set the properties of the matplotlib Axes instance."""
    ylim = kwargs.get('ylim', False)
    yticks = kwargs.get('yticks', False)
    yticklabels = kwargs.get('yticklabels', False)
    ylabel = kwargs.get('ylabel', 'units')
    if ylim:
        ax.set_ylim(ylim)
    if yticks:
        ax.set_yticks(yticks)
    if yticklabels:
        ax.set_yticklabels(yticklabels, fontsize='small')
    if ylabel:
        if ylabel == 'units': ylabel = var.plot_units
        ax_ylabel = ax.set_ylabel(ylabel, fontsize='small', labelpad=-2)
    xlim = kwargs.get('xlim', False)
    xticks = kwargs.get('xticks', False)
    xticklabels = kwargs.get('xticklabels', False)
    xlabel = kwargs.get('xlabel', False)
    if xlim:
        if xlim == 'ann_cycle':
            xlim = (1,12)
            xticks = range(1,13)
            xticklabels = tuple('JFMAMJJASOND')
        ax.set_xlim(xlim)
    if xticks:
        ax_xticks = ax.set_xticks(xticks)
    if xticklabels:
        ax_xticklabels = ax.set_xticklabels(xticklabels, fontsize='x-small')
    if xlabel:
        if xlabel == 'units': xlabel = var.plot_units
        ax_xlabel = ax.set_xlabel(xlabel, fontsize='small', labelpad=1)
    plt.tick_params(labelsize='x-small')
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    # Retain right edge ticks for vertical-1D plots.
    if 'do_vert_centroid' not in kwargs.keys():
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')

def auto_extend_cbar(data, cntr_lvls):
    """
    Determine if colorbar needs to extend past the min and/or max
    designated contour level.
    """
    import numpy as np
    if np.min(data) < cntr_lvls[0]:
        if np.max(data) > cntr_lvls[-1]:
            extend = 'both'
        else:
            extend = 'min'
    else:
        if np.max(data) > cntr_lvls[-1]:
            extend = 'max'
        else:
            extend = 'neither'
    return extend

def make_colorbar(contours, label, for_all_panels=True, 
                  ticks=False, ticklabels=False, extend='both',
                  ax_lim=[0.31,0.08,0.4,0.02]):
    """Create colorbar for multi panel plots."""
    from matplotlib import pyplot as plt
    import numpy as np
    # Goes at bottom center if for all panels.
    if for_all_panels:
        ax_cbar = plt.gcf().add_axes(ax_lim)
	cbar = plt.gcf().colorbar(contours, cax=ax_cbar,
				  orientation='horizontal', drawedges=False,
				  spacing='proportional', extend=extend)
    # Goes within current axis if only for that axis
    else:
        cbar = plt.gcf().colorbar(contours, ax=plt.gca(),
                                  drawedges=False,
                                  orientation='vertical',
                                  spacing='proportional',
                                  extend=extend, fraction=0.1,
                                  aspect=16, pad=0.03)
    # Set tick and label properties.
    if np.any(ticks):
        cbar.set_ticks(ticks)
    if ticklabels not in (None, False):
        cbar.set_ticklabels(ticklabels)
    cbar.ax.tick_params(labelsize='x-small')
    cbar.set_label(label, fontsize='x-small', labelpad=0.5)

def plot_rectangle(bmap, lonmin, lonmax, latmin, latmax):
    """Plot a lon-lat 'rectangle' in Basemap.  Taken from:
http://stackoverflow.com/questions/23941738/python-draw-rectangle-in-basemap
    """
    xs = [lonmin, lonmax, lonmax, lonmin, lonmin]
    ys = [latmin, latmin, latmax, latmax, latmin]
    return bmap.plot(xs, ys, latlon=True)

def plot_lon_vert(ax, proj_name, model_name, run_name, ens_mem, var_name,
                  intvl, yr_range, lon, lev, **kwargs):
    """Plot monthly- and zonal-mean time averaged data as a contour plot."""
    import numpy as np
    from matplotlib import pyplot as plt
    # Load and prep the data.
    proj, model, run, var, data, data_cont, data_stat, \
        ax_fill_color, start_yr, end_yr = prep_plot_data(
            'lon_vert', proj_name, model_name, run_name, False,
            ens_mem, var_name, False, intvl, yr_range, **kwargs
        )
    # Specify contour parameters.
    max_cntr = kwargs.get('max_cntr', 10.)
    min_cntr = kwargs.get('min_cntr', -10.)
    cntr_lvls = np.linspace(min_cntr, max_cntr, kwargs.get('num_cntr', 12) + 1)
    col_map = kwargs.get('col_map', 'RdBu_r')
    cbar_extend = kwargs.get('cbar_extend', 'both')
    lon_bounds = kwargs.get('lon_bounds', [0.,360.])
    lev_type = kwargs.get('lev_type', 'linear')
    lev_bounds = kwargs.get('lev_bounds', [1000., 10.])
    # Specify lon, and level paramters.
    if lon == 'lon':
        lons = model.lon
    # Pivot lon indices in order to span array edge if necessary.
    if lon_bounds[0] < 0 and lons[0] > 0:
        pivot = len(lons)/2
        lons[pivot:] -= 360
        lons = np.roll(lons, pivot, axis=None)
    # Extend by 5 degrees lon beyond bounds to ensure whole region contoured.
    lon_indices = np.where((lons >= lon_bounds[0] - 5) &
                           (lons <= lon_bounds[1] + 5))
    lons = lons[lon_indices]
    if lev == 'level':
        levs = model.level
    if lev_type == 'exp':
        levs = np.log10(levs)
        lev_bounds = np.log10(lev_bounds)
    # Average data over desired latitude range.
    lat_avg_range = kwargs.get('lat_avg_range', (-90, 90))
    lats = np.where((model.lat >= lat_avg_range[0]) &
                    (model.lat <= lat_avg_range[1]))
    def av_lats(data, lats, sfc_area):
        return (np.average(np.squeeze(data[:,lats]), axis=1,
                           weights=np.squeeze(sfc_area[lats,0])))
    data = np.roll(data, pivot, axis=-1)
    data = av_lats(data[:,:,lon_indices], lats, model.sfc_area)
    do_cont_cntrs = kwargs.get('do_cont_cntrs', False)
    if do_cont_cntrs:        
        data_cont = np.roll(data_cont, pivot, axis=-1)            
        data_cont = av_lats(data_cont[:,:,lon_indices], lats, model.sfc_area)
    do_hatch = kwargs.get('do_hatch', False)
    if do_hatch:
        data_stat = np.roll(data_stat, pivot, axis=-1)
        data_stat = av_lats(data_stat[:,:,lon_indices], lats, model.sfc_area)
    # Make contours.
    X, Z = np.meshgrid(lons, levs)
    cs1 = ax.contourf(X, Z, data, cntr_lvls, extend='both')
    cs1.set_cmap(col_map)
    # Make control contours, if desired.
    do_cont_cntrs = kwargs.get('do_cont_cntrs', False)
    if do_cont_cntrs:
        cont_cntr_levs = kwargs.get('cont_cntr_levs', False)
        if cont_cntr_levs:
            cs_cont = ax.contour(X, Z, data_cont, cont_cntr_levs,
                                 colors='0.5', linewidths=0.7)
        else:
            cs_cont = ax.contour(X, Z, data_cont, colors='0.5',
                                 linewidths=0.7)
        # Make control contour labels, if desired.
        cont_cntr_labels = kwargs.get('cont_cntr_labels', False)
        if cont_cntr_labels:
            if cont_cntr_labels is not True:
                fmt=cont_cntr_labels
            else:
                fmt='%1.1f'
            cont_labels = plt.clabel(cs_cont, fontsize=7, fmt=fmt,
                                     inline=True, inline_spacing=0)
    # # Add hatching if desired.
    # do_stat_mask = kwargs.get('do_stat_mask', False)
    # if do_stat_mask not in [False, 'mask', ['mask'], ('mask',)]:
    #     stat_mask_color = kwargs.get('stat_mask_color', 'none')
    #     # 2014-11-14: edgecolor part not working.  Intent is to specify
    #     # color of the hatches.
    #     cs_std = ax.contourf(X, Z, data_stat, stat_mask_cntrs, 
    #                          hatches=stat_mask_hatch, colors='none',
    #                          edgecolor=stat_mask_color)
    # Set axis properties.
    yticks = kwargs.get('yticks', False)
    # If log-p vertical spacing, convert tick locations if needed.
    if lev_type in ['exp', 'log'] and np.min(yticks) > 3:
        yticks = np.log10(yticks)
    yticklabels = kwargs.get('yticklabels', False)
    if type(yticks) not in (None, bool):
        ax_yticks = ax.set_yticks(yticks)
    if yticklabels:
        ax_yticklabels = ax.set_yticklabels(yticklabels, fontsize='x-small')
    ax_ylabel = ax.set_ylabel(kwargs.get('ylabel', ''), 
                              fontsize='small', labelpad=2)
    ax_ylim = ax.set_ylim(lev_bounds)

    ax.set_xlim(lon_bounds)
    xticks = kwargs.get('xticks', False)
    xticklabels = kwargs.get('xticklabels', False)
    if xticks:
        ax.set_xticks(xticks)
    if xticklabels:
        ax.set_xticklabels(xticklabels, fontsize='x-small')
    ax.set_xlabel(kwargs.get('xlabel', ''), fontsize='small', labelpad=2)
    plt.tick_params(labelsize='small')
    # Colorbar.
    do_cbar = kwargs.get('do_colorbar', False)
    if do_cbar:
        cbar_label = kwargs.get('cbar_label', 'units')
        if cbar_label == 'units':
            cbar_label = var.plot_units
        all_panel = True if do_cbar == 'all' else False
        make_colorbar(cs1, cbar_label, for_all_panels=all_panel, 
                      ticks=kwargs.get('cbar_ticks', False),
                      ticklabels=kwargs.get('cbar_ticklabels', False),
                      extend=cbar_extend,
                      ax_lim=kwargs.get('ax_cbar'))
    return cs1

def plot_lat_vert(ax, pre, model, run, var, yr_range, time_label, 
                  lat, lev, **kwargs):
    """Plot latitude v. pressure data as a contour plot."""
    # Load and prep the data.
    proj, model, run, var, data, data_cont, data_stat, \
        ax_fill_color, start_yr, end_yr = prep_plot_data(
            'lat_vert', proj, model_name, run_name, region,
            ens_mem, var_name, level, intvl, yr_range, **kwargs
        )
    # # Get the specified data.
    # data = load_znl_vert(var, pre, model, run, yr_range, time_label)
    # # Convert to plotting units.
    # conv, units = display_params(var)
    # data*=conv
    # # Get control run data if desired.
    # do_cont = kwargs.get('do_cont', False)
    # if do_cont:
    #     data_cont = load_znl_vert(var, pre, model, 'cont', yr_range, time_label)
    #     data_cont*=conv
    # else:
    #     data_cont = 0.
    # Get angular momentum data if desired.
    do_ang_mom = kwargs.get('do_ang_mom', False)
    if do_ang_mom:
        # Use CONT values for run differences, otherwise use that run's values.
        if run not in ['aero', 'aero_tm', 'aero_mtm', 
                       'gas', 'gas_tm', 'gas_mtm']:
            file_ang_mom = (pre + model + '/cont/ang_mom.av.znl.' + 
                            model + '.cont.' + time_label  + '.' + 
                            yr_range + '.p')
        else:
            file_ang_mom = (path + 'ang_mom.av.znl.' + model + '.' + run + 
                            '.'  + time_label  + '.' + yr_range + '.p')
        data_ang_mom = load(open(file_ang_mom, 'r'))/(Omega*r_e**2)

    # Specify contour parameters.
    max_cntr = kwargs.get('max_cntr', 10.)
    min_cntr = kwargs.get('min_cntr', -10.)
    cntr_lvls = np.linspace(min_cntr, max_cntr, kwargs.get('num_cntr', 12) + 1)
    col_map = kwargs.get('col_map', 'RdBu_r')
    lat_type = kwargs.get('lat_type', 'reg')
    lat_bounds = kwargs.get('lat_bounds', [-90.,90.])
    lev_type = kwargs.get('lev_type', 'linear')
    lev_bounds = kwargs.get('lev_bounds', [1000., 10.])
    if lat_type == 'reg':
        lats = lat
    elif lat_type == 'sin':
        lats = np.sin(np.deg2rad(lat))
        lat_bounds = np.sin(np.deg2rad(lat_bounds))
    if lev_type == 'linear':
        levs = lev
    elif lev_type == 'exp':
        levs = np.log10(lev)
        lev_bounds = np.log10(lev_bounds)
    # Make contours.
    Y, Z = np.meshgrid(lats, levs)
    cs1 = ax.contourf(Y, Z, data, cntr_lvls, extend='both')
    cs1.set_cmap(col_map)
    cs2 = ax.contour(Y, Z, data, cntr_lvls)
    cs2.set_cmap(col_map)
    # Compute and plot Hadley cell maximum values if desired.
    do_had_max = kwargs.get('do_had_max', False)
    if do_had_max:
        # Use CONT values for run differences, otherwise use that run's values.
        if run not in ['aero', 'aero_tm', 'aero_mtm', 
                       'gas', 'gas_tm', 'gas_mtm']:
            file_had_bnd = (pre + model + '/cont/msf.av.' + 
                            model + '.cont.' + time_label  + '.' + 
                            yr_range + '.p')
        else:
            file_had_bnd = (path + 'msf.av.' + model + '.' + run + '.'  +
                            time_label  + '.' + yr_range + '.p')
        strmfunc = load(open(file_had_bnd, 'r'))
        # Compute.
        had_max_data = had_bounds(strmfunc, lat, return_max=True)
        z_min = levs[had_max_data[0]]; y_min = lats[had_max_data[1]]
        had_min = had_max_data[2]*1e-10
        z_max = levs[had_max_data[3]]; y_max = lats[had_max_data[4]]
        had_max = had_max_data[5]*1e-10
        # Plot.
        loc_max = ax.plot(y_max, z_max, '.k')
        loc_min = ax.plot(y_min, z_min, '.k')
        dat_max = ax.text(y_max, z_max, '  %0d' % had_max, fontsize='small')
        dat_min = ax.text(y_min, z_min, '  %0d' % had_min, fontsize='small')
    # Plot control contours if desired.
    if do_cont:
        ctr_cont = ax.contour(lats, levs, data_cont, colors='k')
        ctr_cont.clabel(fontsize=8, fmt='%0.f', inline=True, inline_spacing=0)
    # Plot angular momentum contours if desired.
    if do_ang_mom:
        ctr_ang_mom = ax.contour(Y, Z, data_ang_mom, np.linspace(0., 1., 11),
                                 colors='grey', linewidths=0.2)
    # Set axis properties.
    ax.set_ylim(lev_bounds)
    yticks = kwargs.get('yticks', False)
    yticklabels = kwargs.get('yticklabels', False)
    if yticks:
        ax.set_yticks(yticks)
    if yticklabels:
        ax.set_yticklabels(yticklabels, fontsize='small')
    ax.set_ylabel(kwargs.get('ylabel', ''), fontsize='small', labelpad=2)
    ax.set_xlim(lat_bounds)
    xticks = kwargs.get('xticks', False)
    xticklabels = kwargs.get('xticklabels', False)
    if xticks:
        ax.set_xticks(xticks)
    if xticklabels:
        ax.set_xticklabels(xticklabels, fontsize='small')
    ax.set_xlabel(kwargs.get('xlabel', ''), fontsize='small', labelpad=2)
    plt.tick_params(labelsize='small')
    # Colorbar.
    do_cbar = kwargs.get('do_colorbar', False)
    if do_cbar:
        all_panel = True if do_cbar == 'all' else False
        make_colorbar(units, for_all_panel=all_panel,
                      ticks=kwargs.get('cbar_ticks', False),
                      extend=kwargs.get('cbar_extend', 'both'),
                      ax_lim=kwargs.get('ax_cbar'))
    return cs1
        
def plot_lat_ssnl(ax, proj, model, run, ens_mem, var, level, 
		      yr_range, **kwargs):
    """Plot monthly- and zonal-mean time averaged data as a contour plot."""
    import numpy as np
    from matplotlib import pyplot as plt
    from aospy.calcs import (had_bounds, had_bounds500, thermal_equator,
                             itcz_pos, prec_centroid)
    # Load and prep the data.
    proj, model, run, var, data, data_cont, data_stat, \
        ax_fill_color, start_yr, end_yr = prep_plot_data(
            'lat_ssnl', proj, model_name, run_name, region,
            ens_mem, var_name, level, intvl, yr_range, **kwargs
        )
    # # Get aospy objects if only the name was passed.
    # var, run, model, proj = _aospy_inst(var=var, run=run, 
    #                                     model=model, proj=proj)
    # # Load each month's data. If multiple runs take the difference.
    # if type(run) not in [list, tuple]:
    #     run = [run]
    # data = [load_ann_cycle(proj, model, rn, ens_mem, var, level, 'av',  
    #                        yr_range, do_znl_mean=True) for rn in run]
    # if len(data) == 1:
    #     data = data[0]
    # else:
    #     data = np.sum(data[:-1], axis=0) - data[-1]*(len(data[:-1]))
    # run = run[-1]
    # # Convert
    # data = var.convert_to_plot_units(data)
    # Get ITCZ location if desired.
    do_itcz = kwargs.get('do_itcz', False)
    if do_itcz:
        # Use perturbation run values for run differences.
        prec_znl = load_ann_cycle(proj, model, run, ens_mem, 'precip',
				  None, 'av', yr_range, do_znl_mean=True)
        data_itcz = np.array([itcz_pos(prec_znl[t], model.lat) 
			      for t in range(12)])
    # Get Hadley cell boundaries if desired.
    had_bnds_type = kwargs.get('had_bnds_type', False)
    plot_had_bnds = kwargs.get('plot_had_bnds', False)
    mask_hadley = kwargs.get('mask_hadley', False)
    if any(plot_had_bnds) or mask_hadley:
        # if run == 'cont':
        strmfunc = load_ann_cycle(proj, model, run, ens_mem, 'msf',
				  None, 'av', yr_range)
        if had_bnds_type == '500hPa':
            data_had_bnd = np.array([had_bounds500(strmfunc[t], model.lat) 
                                     for t in range(12)]).T
        else:
            data_had_bnd = np.array([had_bounds(strmfunc[t], model.lat) 
                                     for t in range(12)]).T
    # Get thermal equator if desired.
    plot_thermal_eq = kwargs.get('plot_thermal_eq', False)
    if plot_thermal_eq:
        flux_type = plot_thermal_eq
        en_flux = load_ann_cycle(flux_type, pre, model, 
				    'cont', yr_range)
        data_thermal_eq = np.array([thermal_equator(en_flux[t], model.lat) 
                                    for t in range(12)])  
    # Mask values outside Hadley cells and near their inner boundary if desired.
    if mask_hadley:
        if mask_hadley == 'outer':
            had_mask = np.array([np.where((model.lat <= data_had_bnd[0,i])
                                | (model.lat >=  data_had_bnd[2,i]) , 1, 0) 
                                for i in range(12)])
        else: 
            had_mask = np.array([np.where((model.lat <= data_had_bnd[0,i])
                                | (model.lat >=  data_had_bnd[2,i])
                                | (np.abs(model.lat - data_had_bnd[1,i]) < 6.),
					  1, 0) for i in range(12)])
        data = ma.array(data, mask=had_mask)
        ax.set_axis_bgcolor('0.8')
    # Specify contour parameters.
    max_cntr = kwargs.get('max_cntr', 10.)
    min_cntr = kwargs.get('min_cntr', -10.)
    cntr_lvls = np.linspace(min_cntr, max_cntr, kwargs.get('num_cntr', 12) + 1)
    col_map = kwargs.get('col_map', 'RdBu_r')
    lat_type = kwargs.get('lat_type', 'reg')
    lat_bounds = kwargs.get('lat_bounds', [-90.,90.])
    if lat_type == 'reg':
        lats = model.lat
    elif lat_type == 'sin':
        lats = np.sin(np.deg2rad(model.lat))
        lat_bounds = np.sin(np.deg2rad(lat_bounds))
	if do_itcz:
            data_itcz = np.sin(np.deg2rad(data_itcz))
	if any(plot_had_bnds):
            data_had_bnd = np.sin(np.deg2rad(data_had_bnd))
    # Plot ITCZ, Hadley bounds, and/or thermal equator if desired.
    if do_itcz:
        itcz_plot = ax.plot(range(1,13), data_itcz, ':b', linewidth=2)
    for i, bnd in enumerate(plot_had_bnds):
        if bnd:
            ax.plot(range(1,13), data_had_bnd[i], 'grey')
    if plot_thermal_eq:
        thermal_eq_plot = ax.plot(range(1,13), data_thermal_eq, '--k')
    # Create contours.
    cs1 = ax.contourf(range(1,13), lats, data.T, cntr_lvls, extend='both')
    cs1.set_cmap(col_map)
    cs2 = ax.contour(range(1,13), lats, data.T, cntr_lvls, extend='both')
    cs2.set_cmap(col_map)
    # Hatching for statistical significance if desired.
    hatch = kwargs.get('do_hatch', False)
    if hatch:
        # data_ttest = load_ann_cycle(var, pre, model, run, yr_range, 
				    # data_type='ttest')
        cs_ttest = ax.contour(range(1,13), lats, data_ttest.T, 
                              [hatch])
        hatch_colors = 'k'
        hatch_patterns='x'
        #hatch_contour(ax, cs_ttest, hatch_colors, hatch_patterns,
        #              remove_contour=False)
    # Set axis properties.
    ax.set_ylim(lat_bounds)
    yticks = kwargs.get('yticks', False)
    yticklabels = kwargs.get('yticklabels', False)
    if yticks:
        ax.set_yticks(yticks)
    if yticklabels:
        ax.set_yticklabels(yticklabels, fontsize='small')
    ax.set_ylabel(kwargs.get('ylabel', ''), fontsize='small', labelpad=2)
    xticks = kwargs.get('xticks', False)
    xticklabels = kwargs.get('xticklabels', False)
    if xticks:
        ax.set_xticks(xticks)
    if xticklabels:
        ax.set_xticklabels(xticklabels, fontsize='small')
    ax.set_xlabel(kwargs.get('xlabel', ''), fontsize='small', labelpad=2)
    plt.tick_params(labelsize='small')
    # Colorbar.
    do_cbar = kwargs.get('do_colorbar', False)
    if do_cbar:
        all_panel = True if do_cbar == 'all' else False
        make_colorbar(cs1, var.plot_units, for_all_panels=all_panel, 
                      ticks=kwargs.get('cbar_ticks', False),
                      extend=kwargs.get('cbar_extend', 'both'),
                      ax_lim=kwargs.get('ax_cbar'))
    return cs1

def plot_time_1d(ax, proj, model_name, run_name, region, ens_mem,
                 var_name, level, intvl, yr_range, **kwargs):
    """Plot 1d data w/ time as the x-axis."""
    import numpy as np
    from matplotlib import pyplot as plt
    # Load and prep the data.
    proj, model, run, var, data, data_cont, data_stat, \
        ax_fill_color, start_yr, end_yr = prep_plot_data(
            'time_1d', proj, model_name, run_name, region,
            ens_mem, var_name, level, intvl, yr_range, **kwargs
        )
    # Plot data.
    c = kwargs.get('color', 'b')
    ls = kwargs.get('linestyle', '-')
    lw = kwargs.get('linewidth', 2)
    xlim = kwargs.get('xlim', False)
    if xlim == 'ann_cycle':
        xdata = range(1,13)
    else:
        xdata = range(start_yr, end_yr + 1)
    plot_handle = ax.plot(xdata, data, color=c, linestyle=ls, linewidth=lw)
    # Set axis properties.
    _set_axes_props(ax, var, **kwargs)
    # # Horizontal zero line.
    # if ylim and ylim[0] < 0 < ylim[1]:
    #     ax_hline = ax.axhline(0., color='grey')
    return plot_handle


def plot_lon_lat(ax, proj_name, model_name, run_name, ens_mem, var_name, 
                 level, intvl, data_type, yr_range, **kwargs):
    """Plot lat-lon data as a map."""
    import numpy as np
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import addcyclic, shiftgrid, Basemap
    from aospy.io import (load_file, load_plot_data, _var_inst, 
                          _proj_inst, _run_inst, _model_inst)
    from aospy.plotting import make_colorbar
    # Load and prep the data.
    proj, model, run, var, data, data_cont, data_stat, \
        ax_fill_color, start_yr, end_yr = prep_plot_data(
            'lon_lat', proj_name, model_name, run_name, False,
            ens_mem, var_name, level, intvl, yr_range, **kwargs
        )
    # Load data for wind (or wind-like) quiver plot, if desired.
    do_quiver = kwargs.get('do_quiver', False)
    if do_quiver:
        lev_quiv = kwargs.get('level_quiv', 850)
        xvar_quiv = kwargs.get('xvar_quiv', 'ucomp')
        xrun_quiv = kwargs.get('xrun_quiv', False)
        yvar_quiv = kwargs.get('yvar_quiv', 'vcomp')
        yrun_quiv = kwargs.get('yrun_quiv', False)
        if not xrun_quiv:
            xrun_quiv = run
        if not yrun_quiv:
            yrun_quiv = run
        uquiv = load_plot_data(
            proj, model, xrun_quiv, ens_mem, 'ucomp', lev_quiv, intvl, 'av',
            yr_range)
        vquiv = load_plot_data(
            proj, model, yrun_quiv, ens_mem, 'vcomp', lev_quiv, intvl, 'av',
            yr_range)
        scalar_quiv = kwargs.get('scalar_quiv', None)
        if scalar_quiv:
            scalar_run_quiv = kwargs.get('scalar_run_quiv', False)
            if not scalar_run_quiv:
                scalar_run_quiv = run
            squiv = load_plot_data(
                proj, model, scalar_run_quiv, ens_mem, scalar_quiv, lev_quiv,
                intvl, 'av', yr_range)
            uquiv *= squiv
            vquiv *= squiv

    # Specify contour parameters.
    max_cntr = kwargs.get('max_cntr', 10.)
    min_cntr = kwargs.get('min_cntr', -10.)
    cntr_lvls = np.linspace(min_cntr, max_cntr, kwargs.get('num_cntr', 10) + 1)
    col_map = kwargs.get('col_map', 'RdBu_r')
    cbar_extend = kwargs.get('cbar_extend', 'both')
    # Map and contours.
    map_proj = kwargs.get('map_proj', 'moll')
    map_res = kwargs.get('map_res', 'c')
    left_lon = kwargs.get('left_lon', 0.)
    map_corners = kwargs.get('map_corners', {})
    # Shift data as desired.
    if map_proj == 'moll':
        data, lon = shiftgrid(model.lon[0]+180., data, model.lon)
    elif map_proj == 'cyl':
        data, lon = shiftgrid(181., data, model.lon, start=False)
    else:
        lon = model.lon
    bmap = Basemap(projection=map_proj, #lon_0=left_lon,
                   resolution=map_res, ax=ax, **map_corners)
    x, y = bmap(*np.meshgrid(lon, model.lat))
    # Make the contour plot.
    cs1 = bmap.contourf(x, y, data, cntr_lvls, extend=cbar_extend)
    cs1.set_cmap(col_map)
    # Add contours from control run if desired.
    do_cont_cntrs = kwargs.get('do_cont_cntrs', False)
    if do_cont_cntrs:
        if map_proj == 'moll':
            data_cont, lon = shiftgrid(model.lon[0]+180., data_cont, model.lon)
        elif map_proj == 'cyl':
            data_cont, lon = shiftgrid(181., data_cont, model.lon, start=False)
        else:
            lon = model.lon
        cont_cntr_levs = kwargs.get('cont_cntr_levs', False)
        if cont_cntr_levs:
            cs_cont = bmap.contour(x, y, data_cont, cont_cntr_levs,
                                   colors='0.5', linewidths=0.7)
        else:
            cs_cont = bmap.contour(x, y, data_cont, colors='0.5',
                                   linewidths=0.7)
        # Make control contour labels, if desired.
        cont_cntr_labels = kwargs.get('cont_cntr_labels', False)
        if cont_cntr_labels:
            if cont_cntr_labels is not True:
                fmt=cont_cntr_labels
            else:
                fmt='%1.1f'
            cont_labels = plt.clabel(cs_cont, fontsize=7, fmt=fmt,
                                     inline=True, inline_spacing=0)
    # Add hatching if desired.
    do_stat_mask = kwargs.get('do_stat_mask', False)
    stat_mask_hatch = kwargs.get('stat_mask_hatch', ['///'])
    if do_stat_mask and (stat_mask_hatch != 'mask'):
        stat_mask_cntrs = kwargs.get('stat_mask_cntrs', [0.05, 1])
        data_stat, lon = shiftgrid(181., data_stat, model.lon, start=False)
        cs_std = bmap.contourf(x, y, data_stat, stat_mask_cntrs, 
                               hatches=stat_mask_hatch, colors='none')
    # Add wind/flux arrows if desired.
    if do_quiver:
        if map_proj == 'moll':
            uquiv, lon = shiftgrid(model.lon[0]+180., uquiv, model.lon)
            vquiv, lon = shiftgrid(model.lon[0]+180., vquiv, model.lon)
        elif map_proj == 'cyl':
            uquiv, lon = shiftgrid(181., uquiv, model.lon, start=False)
            vquiv, lon = shiftgrid(181., vquiv, model.lon, start=False)
        else:
            lon = model.lon
        qskip = kwargs.get('quiv_skip', 1)
        if type(qskip) is not int:
            qskip = 1
        quiv_kw = kwargs.get('quiv_kwargs', {})
        cs_quiv = bmap.quiver(x[::qskip], y[::qskip], uquiv[::qskip],
                              vquiv[::qskip], **quiv_kw)
        do_quiv_key = kwargs.get('do_quiv_key', False)
        if do_quiv_key:
            scale_qk = kwargs.get('scale_quiv_key', 1)
            lbl_qk = kwargs.get('label_quiv_key', '1')
            quivkey = plt.quiverkey(cs_quiv, -0.1, 0.5, scale_qk, lbl_qk,
                                    fontproperties={'size': 'x-small'})
    # Add any desired lat-lon rectangles, e.g. region boundaries.
    draw_latlon_rect = kwargs.get('draw_latlon_rect', False)
    if draw_latlon_rect:
        plt_rect = plot_rectangle(bmap, *draw_latlon_rect)
    # Fill or outline continents, map border, etc.
    if kwargs.get('fill_continents', False):
        bmap.fillcontinents()
    else:
        bmap.drawcoastlines(linewidth=0.3)
    bmap.drawmapboundary(linewidth=0.1, fill_color=ax_fill_color)
    # Colorbar.
    do_cbar = kwargs.get('do_colorbar', False)
    if do_cbar:
        cbar_label_norm = 'dimensionless'
        cbar_label = kwargs.get('cbar_label', 'units')
        if cbar_label:
            if cbar_label == 'units':
                do_normalize = kwargs.get('do_normalize', False)
                do_norm_cont = kwargs.get('do_norm_cont', False)
                if do_normalize or do_norm_cont:
                    cbar_label = cbar_label_norm
                else:
                    cbar_label = var.plot_units
        else:
            cbar_label = ''
        all_panel = True if do_cbar == 'all' else False
        make_colorbar(cs1, cbar_label, for_all_panels=all_panel, 
                      ticks=kwargs.get('cbar_ticks', False),
                      ticklabels=kwargs.get('cbar_ticklabels', False),
                      extend=cbar_extend,
                      ax_lim=kwargs.get('ax_cbar'))
    return cs1

def plot_lat_1d(ax, proj, model, run, ens_mem, var, x_data, level, 
                intvl, yr_range, **kwargs):
    """Plot 1d data."""
    import numpy as np
    from matplotlib import pyplot as plt
    from aospy.io import load_file
    from aospy.io import _var_inst, _proj_inst, _run_inst, _model_inst
    # Get aospy objects if only the name was passed.
    proj = _proj_inst(proj)
    model = _model_inst(model, proj)
    run = _run_inst(run, model, proj)
    var = _var_inst(var)    
    # Get year range.
    if yr_range == 'default':
        yr_range = model.default_yr_range
    num_yr = yr_range[1] - yr_range[0] + 1
    # Load  data. If multiple runs take the difference.
    if type(run) not in [list, tuple]:
        data = load_file(proj, model, run, ens_mem, var, 
			 level, intvl, 'av', yr_range)
    else:
        data = [load_file(proj, model, rn, ens_mem, var, 
			  level, intvl, 'av', yr_range)
		for rn in run]
	data = data[0] - data[1]
	run = run[0]
    if var.def_lon:
        data = np.mean(data, axis=-1)
    # Convert to plotting units.
    try:
        data *= var.plot_units_conv
    except:
        pass
    # Get lat data.
    if x_data == 'lat':
        x_data = model.lat
        if kwargs.get('lat_type', 'reg') == 'sin':
            x_data = np.sin(np.deg2rad(x_data))
    # Get parameters for linestyle, color, etc.
    c = kwargs.get('color', 'b')
    ls = kwargs.get('linestyle', '-')
    lw = kwargs.get('linewidth', 2)
    # Plot confidence interval/errorbars if desired.
    if kwargs.get('do_errorbar', False):
        conf_thresh = kwargs.get('conf_thresh', 0.95)
        interv = t.interval(conf_thresh, num_yr - 1)[1]
        try:
            name_std = var_label + '.znl.std.' + suffix
            std = load(open(path_in + name_std, 'r'))
        except:
            name_std = var_label + '.std.' + suffix
            std = load(open(path_in + name_std, 'r'))
        try:
            conf_int = std*interv*var.plot_units_conv
        except:
            conf_int = std*interv
        # Semi transparent shaded confidence interval.
        ax.fill_between(x_data, data-conf_int, data+conf_int, 
                        alpha=0.5, color=c[i])
    # Plot data.
    plot_handle = ax.plot(x_data, data, color=c, linestyle=ls, linewidth=lw)
    # Set axis properties.
    ylim = kwargs.get('ylim', False)
    yticks = kwargs.get('yticks', False)
    yticklabels = kwargs.get('yticklabels', False)
    ylabel = kwargs.get('ylabel', 'units')
    if ylim:
        ax.set_ylim(ylim)
    if yticks:
        ax.set_yticks(yticks)
    if yticklabels:
        ax.set_yticklabels(yticklabels, fontsize='small')
    if ylabel:
        if ylabel == 'units': ylabel = var.plot_units
        ax_ylabel = ax.set_ylabel(ylabel, fontsize='small', labelpad=-2)
    xticks = kwargs.get('xticks', False)
    xticklabels = kwargs.get('xticklabels', False)
    xlim = kwargs.get('xlim', False)
    if xticks:
        ax_xticks = ax.set_xticks(xticks)
    if xticklabels:
        ax_xticklabels = ax.set_xticklabels(xticklabels, fontsize='small')
    if xlim:
        ax_xlim = ax.set_xlim(xlim)
    ax_xlabel = ax.set_xlabel(kwargs.get('xlabel', ''), 
                              fontsize='small', labelpad=2)
    plt.tick_params(labelsize='small')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    # Create legend if desired.
    if kwargs.get('do_legend', False):
        legend_labels = kwargs.get('legend_labels', False)
        leg = ax.legend(legend_labels, frameon=False, loc=0, fontsize='small')
    # Horizontal and vertical zero lines.
    if not xlim or xlim[0] < 0 < xlim[1]:
        ax_vline = ax.axvline(0., color='grey')
    if not ylim or ylim[0] < 0 < ylim[1]:
        ax_hline = ax.axhline(0., color='grey')
    return plot_handle

def plot_1d_vert(ax, proj, model_name, run_name, region_name, ens_mem, 
                 var_name, z_data, level, intvl, yr_range, **kwargs):
    """Plot 1d data."""
    import numpy as np
    from matplotlib import pyplot as plt
    from aospy import calc_levs_thick
    from aospy.calcs import vert_centroid
    # Load and prep the data.
    proj, model, run, var, data, data_cont, data_stat, \
        ax_fill_color, start_yr, end_yr = prep_plot_data(
            '1d_vert', proj, model_name, run_name, region_name,
            ens_mem, var_name, level, intvl, yr_range, **kwargs
        )
    # Get vertical coordinate data, converting from Pa to hPa if needed.
    if z_data == 'level':
        z_data = model.level
        if np.max(z_data) >= 1e4:
            z_data *= 1e-2
    # Get parameters for linestyle, color, etc.
    c = kwargs.get('color', 'b')
    ls = kwargs.get('linestyle', '-')
    lw = kwargs.get('linewidth', 2)
    # Plot error/stdev shading, if desired.
    do_plot_shading = kwargs.get('do_plot_shading', False)
    if do_plot_shading:
        shade_handle = ax.fill_betweenx(z_data, data-data_stat, data+data_stat,
                                        color=c, alpha=0.2)
    # Plot data.
    plot_handle = ax.plot(data, z_data, color=c, linestyle=ls, linewidth=lw)
    # Plot vertical centroid of data on left axis, if desired.
    do_vert_centroid = kwargs.get('do_vert_centroid', False)
    if type(do_vert_centroid) in (tuple, list):
        p_top, p_bot = np.sort(do_vert_centroid)
    else:
        p_top, p_bot = 100., 925.
    if do_vert_centroid:
        centroid_xpos = kwargs.get('centroid_xpos', 0)
        centroid = vert_centroid(data, z_data, p_top=p_top, p_bot=p_bot)
        centroid_handle = ax.plot(centroid_xpos, centroid, '_', color=c,
                               markersize=15)
    # Plot average of data on bottom axis, if desired.
    do_plot_avg = kwargs.get('do_plot_avg', False)
    if do_plot_avg:
        # Average over specified vertical range.
        z_small, z_large = sorted(kwargs.get('plot_avg_levels', (150, 850)))
        avg_levels = np.where((z_data >= z_small) & (z_data <= z_large))
        dz = calc_levs_thick(z_data)
        x = np.ma.average(data[avg_levels], weights=dz[avg_levels])
        ax.plot(x, ylim[0], 'o', color=c)
    zero_line = ax.axvline(0., color='grey')
    _set_axes_props(ax, var, **kwargs)
    return plot_handle

def plot_time_vert(ax, proj, model_name, run_name, region, ens_mem,
                   var_name, level, yr_range, **kwargs):
    """Plot vertically defined data with time as x axis."""
    import numpy as np
    from matplotlib import pyplot as plt
    from aospy.plotting import make_colorbar
    # Load and prep the data.
    time_lim = kwargs.get('time_lim', 'ann_cycle')
    proj, model, run, var, data, data_cont, data_stat, \
        ax_fill_color, start_yr, end_yr = prep_plot_data(
            'time_vert', proj, model_name, run_name, region,
            ens_mem, var_name, level, time_lim, yr_range, **kwargs
        )
    # Specify contour parameters.
    max_cntr = kwargs.get('max_cntr', 10.)
    min_cntr = kwargs.get('min_cntr', -10.)
    cntr_lvls = np.linspace(min_cntr, max_cntr, kwargs.get('num_cntr', 12) + 1)
    col_map = kwargs.get('col_map', 'RdBu_r')
    lev_type = kwargs.get('lev_type', 'linear')
    lev_bounds = kwargs.get('lev_bounds', [1000., 10.])
    if level == None:
        level = model.level
    if level[0] > 1e4:
        level*=1e-2             # Convert Pa to hPa.
    if lev_type == 'linear':
        levels = level
    elif lev_type in ['exp', 'log']:
        levels = np.log10(level)
        lev_bounds = np.log10(lev_bounds)
    cbar_extend = kwargs.get('cbar_extend', 'both')
    if cbar_extend == 'auto':
        cbar_extend = auto_extend_cbar(data, cntr_lvls)
    # Make filled contours of main data.
    Y, Z = np.meshgrid(range(1,13), levels)
    cs1 = ax.contourf(Y, Z, data.T, cntr_lvls, extend=cbar_extend)
    cs1.set_cmap(col_map)
    # Make control contours, if desired.
    do_cont_cntrs = kwargs.get('do_cont_cntrs', False)
    if do_cont_cntrs:
        cont_cntr_levs = kwargs.get('cont_cntr_levs', False)
        if cont_cntr_levs:
            cs_cont = ax.contour(Y, Z, data_cont.T, cont_cntr_levs,
                                 colors='0.5', linewidths=0.7)
        else:
            cs_cont = ax.contour(Y, Z, data_cont.T, colors='0.5',
                                 linewidths=0.7)
        # Make control contour labels, if desired.
        cont_cntr_labels = kwargs.get('cont_cntr_labels', False)
        if cont_cntr_labels:
            if cont_cntr_labels is not True:
                fmt=cont_cntr_labels
            else:
                fmt='%1.1f'
            cont_labels = plt.clabel(cs_cont, fontsize=7, fmt=fmt,
                                     inline=True, inline_spacing=0)
    # Add hatching if desired.
    do_stat_mask = kwargs.get('do_stat_mask', False)
    stat_mask_hatch = kwargs.get('stat_mask_hatch', ['///'])
    if do_stat_mask and (stat_mask_hatch != 'mask'):
        stat_mask_cntrs = kwargs.get('stat_mask_cntrs', [0.05, 1])
        stat_mask_color = kwargs.get('stat_mask_color', 'none')
        # 2014-11-14: edgecolor part not working.  Intent is to specify
        # color of the hatches.
        cs_std = ax.contourf(Y, Z, data_stat.T, stat_mask_cntrs, 
                             hatches=stat_mask_hatch, colors='none',
                             edgecolor=stat_mask_color)
    yticks = kwargs.get('yticks', False)
    # If log-p vertical spacing, convert tick locations if needed.
    if lev_type in ['exp', 'log'] and np.min(yticks) > 3:
        yticks = np.log10(yticks)
    yticklabels = kwargs.get('yticklabels', False)
    if type(yticks) not in (None, bool):
        ax_yticks = ax.set_yticks(yticks)
    if yticklabels:
        ax_yticklabels = ax.set_yticklabels(yticklabels, fontsize='small')
    ax_ylabel = ax.set_ylabel(kwargs.get('ylabel', ''), 
                              fontsize='small', labelpad=2)
    ax_ylim = ax.set_ylim(lev_bounds)
    ax.set_xlim([1,12])
    xticks = kwargs.get('xticks', False)
    xticklabels = kwargs.get('xticklabels', False)
    if xticks:
        ax.set_xticks(xticks)
    if xticklabels:
        ax.set_xticklabels(xticklabels, fontsize='small')
    ax.set_xlabel(kwargs.get('xlabel', ''), fontsize='small', labelpad=2)
    plt.tick_params(labelsize='small')
    # Colorbar.
    do_cbar = kwargs.get('do_colorbar', False)
    if do_cbar:
        cbar_label = kwargs.get('cbar_label', 'units')
        if cbar_label == 'units':
            cbar_label = var.plot_units
        all_panel = True if do_cbar == 'all' else False
        make_colorbar(cs1, cbar_label, for_all_panels=all_panel, 
                      ticks=kwargs.get('cbar_ticks', False),
                      ticklabels=kwargs.get('cbar_ticklabels', False),
                      extend=cbar_extend,
                      ax_lim=kwargs.get('ax_cbar'))
    return cs1
