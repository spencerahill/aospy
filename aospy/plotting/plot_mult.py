import numpy as np
from aospy.plotting import Figure
# from aospy.plotting import (
#     plot_lon_lat, plot_lon_vert, plot_lat_1d, plot_lat_vert, plot_lat_ssnl, 
#     plot_time_1d, plot_time_vert, plot_1d_vert
# )
fig = Figure(
    n_row=1,
    n_col=1,
    n_ax='all',
    n_plot=1,
    n_data=1,

    row_size=1.25,
    col_size=4.5,

    # proj=['obs'] + ['aero_3agcm']*3,
    proj='aero_3agcm',
    # model=[[['cru', 'merra']], 'am2', 'am3', 'hiram'],
    # model=[[['merra', 'cru']], 'am2', 'am3', 'hiram'],
    model='am2',
    intvl='jas',
    # run=[
        # [['v3.22', 'merra']], 
        # [['merra', 'v3.22']], 
        # [['amip', 'amip'], [('reyoi+2K', 'reyoi_cont'), ('reyoi+2K', 'reyoi_cont')]], 
        # [['amip', 'amip'], [('hurrell+2K', 'hurrell_cont'), ('hurrell+2K', 'hurrell_cont')]], 
        # [['amip', 'amip'], [('gas_tm', 'cont'), ('gas_tm', 'cont')]]
    # ],
    run='reyoi_cont',
    # ens_mem=[None]*2 + [['avg', None]]*2,
    ens_mem=None,
    # var=[[['sphum', 't_surf']]],
    # var=[[['precip', 'sphum']]],
    var='toa_rad',
    dtype='av',
    # dtype=['reg.ts'] + [[['reg.ts', 'reg.ts'], ['reg.av','reg.av']]]*3,
    level=None,
    # level=[[[300, None]]],
    # level=[[[None, 300]]],
    region=None,
    # region='sahel',
    yr_range='default',
    # yr_range=[[[(1979, 2011), 'default']]] + ['default']*3,
    # yr_range=[[['default', (1979, 2011)]]] + ['default']*3,
    plot_type='map',

    ## Titles and labels
    fig_title=r'AM2 control JAS TOA radiation',
    # fig_title=r'Sahel JAS $T_s$ v. 300 hPa $q$',    
    # ax_title=['MERRA & CRU', 'AM2.1', 'AM3', 'HiRAM'],
    ax_title=False,
    ax_left_labels=False,
    ax_left_label_rot='vertical', # 'horizontal', 'vertical'
    ax_right_labels=False,

    ## Axis limits
    xlim=(-180,180),
    # xlim=(-2.5, 2.5),
    # xlim=False,
    xticks=False,
    xticklabels=False,
    # xlabel=r'Kelvin',
    # xlabel=r'W m$^{-2}$',
    # xlabel=r'g kg$^{-1}$',
    # xlabel=r'mm day$^{-1}$',
    ylim=(-30,30),
    # ylim=(-2.5, 2.5),
    # ylim=(-5, 5),
    # ylim=False,
    yticks=False,
    yticklabels=False,
    # ylabel=r'Kelvin',
    # ylabel=r'g kg$^{-1}$',
    # ylabel='W m$^{-2}$',
    # ylabel=r'mm day$^{-1}$',
    # ylabel='',

    min_cntr=-150,
    max_cntr=150,
    num_cntr=15,
    contourf_extend=True,
    do_colorbar='column',              # 'all' 'column' False True
    # cbar_ticks=False,
    # cbar_ticks=range(-7,8,2),
    col_map='RdBu_r',
    # 'RdBu_r' 'coolwarm' 'Greys' 'BrBG'
    # col_map=[['RdBu_r'], ['BrBG']],
    # cbar_extend: 'auto' 'neither' 'min' 'max' 'both'
    # cbar_extend='both',
    # cbar_label='units',
    # cbar_label=r'Left panels: Kelvin. Right panels: percent'
    # cbar_ticklabels=False,
    # do_cont_cntrs=True,
    # cont_run=([['reyoi_cont']] + [['hurrell_cont']] + ,
    #             [['cont']] + [['ming0']]),
    # cntrs_color='0.85',
    # cont_cntr_levs=False,
    # cont_cntr_levs=range(-60, 61, 10), # omega; cre_net/sw/lw; shflx,
    # cont_cntr_levs=range(-30, 31, 10), # ucomp,
    # cont_cntr_levs=range(-30, 31, 3), # vort, pv,
    # cont_cntr_levs=range(-20, 21, 4), # p-e,
    # cont_cntr_levs=range(-10, 51, 10), # gms_each_level,
    # cont_cntr_levs=range(-4, 5, 1), # vcomp,
    # cont_cntr_levs=list(np.arange(-2, 2.1,0.5)),
    # cont_cntr_levs=list(np.arange(-2, 2.51, 0.5)), # tdt_conv/ls,
    # cont_cntr_levs=list(np.arange(0.5, 6.6, 1)), # mc,
    # cont_cntr_levs=range(1, 20, 3), # sphum,
    # cont_cntr_levs=range(1, 12, 2), # evap,
    # cont_cntr_levs=range(2,26,3), # precip,
    # cont_cntr_levs=range(5, 36, 5), # cld_amt,
    # cont_cntr_levs=range(10, 101, 20), # rh,
    # cont_cntr_levs=range(210,311,20), # temp,
    # cont_cntr_levs=range(200,331,5), # t_surf,
    # cont_cntr_levs=range(200, 361, 10), # equiv_pot_temp; virt_pot_temp; MSE,
    # cont_cntr_levs=range(990,1041,3), # slp,
    # cont_cntr_levs=range(200, 301, 20), # temp,
    # cont_cntr_levs=[[range(200,331,5)],[range(10,101,10)]],
    # cont_cntr_labels: 'True', '%d' for integer, '%0.1f' for decimal,
    # cont_cntr_labels='%d',

    ## Latitude plot parameters.
    # x_data=[['lat']],
    # lat_bounds=[-90,90],
    # lat_ticks=range(-90,91,30),
    # lat_ticks=list(np.sin(np.deg2rd(np.arange(-90., 91., 10.)))),
    # lat_labels=[r'90$^\circ$S', '', '', '', '', '', '30$^\circ$S', '', '', 'EQ', ,
               # '', '', '30$^\circ$N', '', '', '', '', '', '90$^\circ$N']
    # lat_labels=[r'90$^\circ$S', '60$^\circ$S', '30$^\circ$S', 'EQ',,
    #               '30$^\circ$N', '60$^\circ$N', '90$^\circ$N'],

    ## Longitude plot parameters.
    # x_data=[['lon']],
    # lon_bounds=(-90,60),
    # lat_avg_range=(10, 20),
    # lon_ticks=range(Fig.lon_bounds[0], Fig.lon_bounds[1]+1, 30),
    # lon_ticks=False,

    # def make_lon_ticklabels(lon_ticks, circ=True):,
    #     lon_labels=[],
    #     c=r'$^\circ$' if circ else '',
    #     for tick in lon_ticks:,
    #         if tick < 0:,
    #             lon_labels.append(r'%d%sW' % (abs(tick), c)),
    #         elif tick > 0:,
    #             lon_labels.append(r'%d%sE' % (abs(tick), c)),
    #         else:,
    #             lon_labels.append(r"0%s" % c),
    #     return lon_labels,

    # lon_labels=make_lon_ticklabels(lon_ticks),
    # lon_labels=[str(x) for x in lon_ticks],
    # lon_labels=False,

    ## Zonal-seasonal and zonal-1d plot parameters
    # do_itcz=False,
    # had_bnds_type='500hPa',
    # plot_had_bnds=[False, False, False],
    # mask_hadley=True,
    # lat_type='reg',

    ## Vert plot parameters.
    # z_data=[['level']],
    # lev_type='linear',             # 'linear' 'exp',
    # lev_bounds=[1000, 100],
    # lev_ticks=range(1000, 99, -100),
    # lev_ticklabel=['1000', '', '800', '', '600', '', '400', '', '200', ''],
    # lev_label='hPa',

    # do_vert_centroid=(100, 1000),
    # centroid_xpos=0,
    # centroid_xpos=[[xlim[0]]*4 + [xlim[0] + 5]*2],
    # do_plot_avg=False,
    # do_plot_shading=True,

    ## Time-1D plot paramters
    # time_lim='ann_cycle',          # 'ann_cycle',
    # time_ticks=False,
    # time_ticklabels=False,

    ## Map plot parameters.
    map_proj='cyl',                # 'moll',
    # map_corners=False,
    # map_corners={'llcrnrlon': -60, 'llcrnrlat': -10,,
                   # 'urcrnrlon': 80, 'urcrnrlat': 40}
    map_res='c',
    # left_lon=0,
    shiftgrid_start=False,
    shiftgrid_cyclic=360.0,
    # draw_latlon_rect=[-18, 40, 10, 20],

    ## Quiver (i.e. arrow) plot overlayed on maps.
    do_quiver=False,
    # level_quiv=level,
    # xvar_quiv='ucomp',
    # yvar_quiv='vcomp',
    # xrun_quiv=cont_run,
    # xrun_quiv=([[['reyoi+2K', 'reyoi_cont']]] +,
            # [[['hurrell+2K', 'hurrell_cont']]] +
            # [[['gas_tm', 'cont']]] +
            # [[['ming0_p2K', 'ming0']]])
    # yrun_quiv=xrun_quiv,
    # scalar_quiv=None,
    # quiv_skip=[[1],[1],[10],[1]],
    # quiv_skip=False,
    # quiv_kwargs={},
    # quiv_kwargs={'color': 'red', 'edgecolors': 'red'},
    # do_quiv_key=True,
    # scale_quiv_key=10,
    # label_quiv_key=str(scale_quiv_key) + r' g kg$^{-1}$ m s$^{-1}$',
    # label_quiv_key=str(scale_quiv_key) + r' m s$^{-1}$',

    ## Scatter plot parameters
    # marker_size=10,
    # marker_color='k',
    # marker_shape='.',
    # marker_size=[10] + [[10, 30]]*3,
    # marker_color=['k'] + [['k', 'r']]*3,
    # marker_shape=['.'] + [['.', 's']]*3,
    # do_best_fit_line=[True] + [[True, False]]*3,
    # print_best_fit_slope=[True] + [[True, False]]*3,
    # print_corr_coeff=[True] + [[True, False]]*3,
    
    ## Statistical masking parameters.
    # do_stat_mask='ttest'             # False 'ttest',
    # do_stat_mask=[[False], ['ttest']],
    # stat_mask_cntrs=[0, 0.05] # Mark where significant.,
    # stat_mask_cntrs=[0.05, 1] # Mark where insignificant.,
    # stat_mask_hatch=['///'] # ['...'] ['XXX'],
    # stat_mask_hatch='mask',
    # stat_mask_color='b', # 2015-02-27: not yet working,

    ## Transformations to apply to data.
    # do_subtract_mean=[True, [True, False], [True, False], [True, False]],
    do_subtract_mean=False,
    # do_zasym=False,
    # do_norm_cont=False,
    # do_normalize=False,
    # norm_region='sahel',
    # norm_var='precip,' # False
    # norm_run=cont_run, # runs
    # norm_run=([[['reyoi+2K', 'reyoi_cont']]] +
            # [[['hurrell+2K', 'hurrell_cont']]] +
            # [[['gas_tm', 'cont']]] +
            # [[['ming0_p2K', 'ming0']]])

    ## 1d plot parameters
    # colors=['brown', 'orange', 'cyan', 'blue', 'gray'],
    # linestyles='-', # [':', '-', '--'] 
    # linewidths=[1.5],
    # do_legend=True,             # 'first', True
    # legend_labels='models',
    # legend_labels=['AM2.1', 'AM3', 'HiRAM', 'c48-HiRAM', 'ERA-I'],
    # legend_labels=['ERA-I', 'MERRA', 'NCEP-CFSR'],
    # Dummy kwarg so that comma on last un-commented kwarg above doesn't
    # cause syntax error.
    dummy=None
)

## Make the plot.
if __name__ == '__main__':
    fig.create_fig()
    fig.make_plots()
    plt.show()

####################################

# znl_kwargs_l = {
#     'do_errorbar': False, 'lat_type': lat_type, 'xlim': lat_bounds,
#      'xticks': lat_ticks, 'xticklabels': ' ', 'xlabel': ' ', 'ylim':
#      ylim, 'yticks': False, 'yticklabels': False, 'ylabel': 'units',
#      'color': colors, 'linestyle': linestyles, 'linewidth': linewidths
# }
# znl_kwargs_bl = {
#     'do_errorbar': False, 'lat_type': lat_type, 'xlim': lat_bounds,
#       'xticks': lat_ticks, 'xticklabels': lat_labels, 'xlabel': ' ',
#       'ylim': ylim, 'yticks': False, 'yticklabels': False, 'ylabel':
#       'units', 'color': colors, 'linestyle': linestyles, 'linewidth':
#       linewidths
# }
# znl_kwargs_b = {
#     'do_errorbar': False, 'lat_type': lat_type, 'xlim': lat_bounds,
#      'xticks': lat_ticks, 'xticklabels': lat_labels, 'xlabel': ' ',
#      'ylim': ylim, 'yticks': False, 'yticklabels': ' ', 'ylabel':
#      None, 'color': colors, 'linestyle': linestyles, 'linewidth':
#      linewidths
# }
# znl_kwargs_int = {
#     'do_errorbar': False, 'lat_type': lat_type, 'xlim': lat_bounds,
#      'xticks': lat_ticks, 'xticklabels': ' ', 'xlabel': ' ', 'ylim':
#      ylim, 'yticks': False, 'yticklabels': ' ', 'ylabel': None,
#      'color': colors, 'linestyle': linestyles, 'linewidth': linewidths
# }
# # 1d vertical plots
# vert_1d_kwargs_l = {
#     'do_errorbar': False, 'lev_type': lev_type, 'xlim': xlim,
#     'xticks': False, 'xticklabels': ' ', 'xlabel': ' ', 'ylim':
#     lev_bounds, 'yticks': lev_ticks, 'yticklabels': lev_ticklabel,
#     'ylabel': lev_label, 'color': colors, 'linestyle': linestyles,
#     'linewidth': linewidths, 'do_vert_centroid': do_vert_centroid,
#     'centroid_xpos': centroid_xpos, 'do_plot_avg': do_plot_avg,
#     'do_plot_shading': do_plot_shading
# }
# vert_1d_kwargs_bl = {
#     'do_errorbar': False, 'lev_type': lev_type, 'xlim': xlim,
#     'xticks': False, 'xticklabels': False, 'xlabel': xlabel, 'ylim':
#     lev_bounds, 'yticks': lev_ticks, 'yticklabels': lev_ticklabel,
#     'ylabel': lev_label, 'color': colors, 'linestyle': linestyles,
#     'linewidth': linewidths, 'do_vert_centroid': do_vert_centroid,
#     'centroid_xpos': centroid_xpos, 'do_plot_avg': do_plot_avg,
#     'do_plot_shading': do_plot_shading
# }
# vert_1d_kwargs_b = {
#     'do_errorbar': False, 'lev_type': lev_type, 'xlim': xlim,
#      'xticks': False, 'xticklabels': False, 'xlabel': xlabel, 'ylim':
#      lev_bounds, 'yticks': lev_ticks, 'yticklabels': ' ', 'ylabel':
#      '', 'color': colors, 'linestyle': linestyles, 'linewidth':
#      linewidths, 'do_vert_centroid': do_vert_centroid,
#      'centroid_xpos': centroid_xpos, 'do_plot_avg': do_plot_avg,
#     'do_plot_shading': do_plot_shading
# }
# vert_1d_kwargs_int = {
#     'do_errorbar': False, 'lev_type': lev_type, 'xlim': xlim,
#      'xticks': False, 'xticklabels': ' ', 'xlabel': ' ', 'ylim':
#      lev_bounds, 'yticks': lev_ticks, 'yticklabels': ' ', 'ylabel':
#      '', 'color': colors, 'linestyle': linestyles, 'linewidth':
#      linewidths, 'do_vert_centroid': do_vert_centroid,
#      'centroid_xpos': centroid_xpos, 'do_plot_avg': do_plot_avg,
#     'do_plot_shading': do_plot_shading
# }
# # Lat-lon (e.g. map) plots.
# map_kwargs_l = {
#     'is_SST': False, 'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr,
#     'col_map': col_map, 'do_colorbar': False, 'cbar_label': cbar_label,
#     'cbar_ticks': cbar_ticks, 'cbar_extend': cbar_extend, 'do_stat_mask':
#     do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs, 'stat_mask_hatch':
#     stat_mask_hatch, 'do_cont_cntrs': do_cont_cntrs, 'cont_cntr_levs':
#     cont_cntr_levs, 'cont_cntr_labels': cont_cntr_labels, 'cont_run': cont_run,
#     'stat_mask_color': stat_mask_color, 'map_corners': map_corners, 'map_proj':
#     map_proj, 'left_lon': left_lon, 'do_subtract_mean': do_subtract_mean,
#     'do_zasym': do_zasym, 'do_norm_cont': do_norm_cont, 'do_normalize':
#     do_normalize, 'norm_region': norm_region, 'norm_run': norm_run,
#     'draw_latlon_rect': draw_latlon_rect, 'do_quiver': do_quiver, 'xvar_quiv':
#     xvar_quiv, 'yvar_quiv': yvar_quiv, 'xrun_quiv': xrun_quiv, 'yrun_quiv':
#     yrun_quiv, 'level_quiv': level_quiv, 'scalar_quiv': scalar_quiv,
#     'quiv_skip': quiv_skip, 'quiv_kwargs': quiv_kwargs, 'do_quiv_key':
#     do_quiv_key, 'scale_quiv_key': scale_quiv_key, 'label_quiv_key':
#     label_quiv_key
# }
# map_kwargs_bl = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'do_colorbar': 'all', 'cbar_label': cbar_label, 'cbar_ticks': cbar_ticks,
#     'cbar_extend': cbar_extend, 'do_stat_mask': do_stat_mask,
#     'stat_mask_cntrs': stat_mask_cntrs, 'stat_mask_hatch': stat_mask_hatch,
#     'do_cont_cntrs': do_cont_cntrs, 'cont_cntr_levs': cont_cntr_levs,
#     'cont_cntr_labels': cont_cntr_labels, 'cont_run': cont_run,
#     'stat_mask_color': stat_mask_color, 'map_corners': map_corners, 'map_proj':
#     map_proj, 'left_lon': left_lon, 'do_subtract_mean': do_subtract_mean,
#     'do_zasym': do_zasym, 'do_norm_cont': do_norm_cont, 'do_normalize':
#     do_normalize, 'norm_region': norm_region, 'norm_run': norm_run,
#     'draw_latlon_rect': draw_latlon_rect, 'do_quiver': do_quiver, 'xvar_quiv':
#     xvar_quiv, 'yvar_quiv': yvar_quiv, 'xrun_quiv': xrun_quiv, 'yrun_quiv':
#     yrun_quiv, 'level_quiv': level_quiv, 'scalar_quiv': scalar_quiv,
#     'quiv_skip': quiv_skip, 'quiv_kwargs': quiv_kwargs, 'do_quiv_key':
#     do_quiv_key, 'scale_quiv_key': scale_quiv_key, 'label_quiv_key':
#     label_quiv_key
# }
# map_kwargs_b = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'do_colorbar': False, 'cbar_label': cbar_label, 'cbar_ticks': cbar_ticks,
#     'cbar_extend': cbar_extend, 'do_stat_mask': do_stat_mask,
#     'stat_mask_cntrs': stat_mask_cntrs, 'stat_mask_hatch': stat_mask_hatch,
#     'stat_mask_color': stat_mask_color, 'do_cont_cntrs': do_cont_cntrs,
#     'cont_cntr_levs': cont_cntr_levs, 'cont_cntr_labels': cont_cntr_labels,
#     'cont_run': cont_run, 'map_corners': map_corners, 'map_proj': map_proj,
#     'left_lon': left_lon, 'do_subtract_mean': do_subtract_mean, 'do_zasym':
#     do_zasym, 'do_norm_cont': do_norm_cont, 'do_normalize': do_normalize,
#     'norm_region': norm_region, 'norm_run': norm_run, 'draw_latlon_rect':
#     draw_latlon_rect, 'do_quiver': do_quiver, 'xvar_quiv': xvar_quiv,
#     'yvar_quiv': yvar_quiv, 'xrun_quiv': xrun_quiv, 'yrun_quiv': yrun_quiv,
#     'level_quiv': level_quiv, 'scalar_quiv': scalar_quiv, 'quiv_skip':
#     quiv_skip, 'quiv_kwargs': quiv_kwargs, 'do_quiv_key': do_quiv_key,
#     'scale_quiv_key': scale_quiv_key, 'label_quiv_key': label_quiv_key
# }
# map_kwargs_int = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'do_colorbar': 'all', 'cbar_label': cbar_label, 'cbar_ticks': cbar_ticks,
#     'cbar_extend': cbar_extend, 'do_stat_mask': do_stat_mask,
#     'stat_mask_cntrs': stat_mask_cntrs, 'stat_mask_hatch': stat_mask_hatch,
#     'stat_mask_color': stat_mask_color, 'do_cont_cntrs': do_cont_cntrs,
#     'cont_cntr_levs': cont_cntr_levs, 'cont_cntr_labels': cont_cntr_labels,
#     'cont_run': cont_run, 'map_corners': map_corners, 'map_proj': map_proj,
#     'left_lon': left_lon, 'do_subtract_mean': do_subtract_mean, 'do_zasym':
#     do_zasym, 'do_norm_cont': do_norm_cont, 'do_normalize': do_normalize,
#     'norm_region': norm_region, 'norm_run': norm_run, 'draw_latlon_rect':
#     draw_latlon_rect, 'do_quiver': do_quiver, 'xvar_quiv': xvar_quiv,
#     'yvar_quiv': yvar_quiv, 'xrun_quiv': xrun_quiv, 'yrun_quiv': yrun_quiv,
#     'level_quiv': level_quiv, 'scalar_quiv': scalar_quiv, 'quiv_skip':
#     quiv_skip, 'quiv_kwargs': quiv_kwargs, 'do_quiv_key': do_quiv_key,
#     'scale_quiv_key': scale_quiv_key, 'label_quiv_key': label_quiv_key
# }
# # Meridional-vertical plots.
# vert_kwargs_l = {
#     'do_ang_mom': True, 'do_had_max': True, 'min_cntr': minc,
#     'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lat_type': lat_type, 'lat_bounds': lat_bounds, 'lev_type':
#     lev_type, 'lev_bounds': lev_bounds, 'xticks': lat_ticks,
#     'xticklabels': None, 'yticks': lev_ticks, 'yticklabels': False,
#     'do_colorbar': False, 'cbar_ticks': cbar_ticks
# }
# vert_kwargs_bl = {
#     'do_ang_mom': True, 'do_had_max': True, 'min_cntr': minc,
#     'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lat_type': lat_type, 'lat_bounds': lat_bounds, 'lev_type':
#     lev_type, 'lev_bounds': lev_bounds, 'xticks': lat_ticks,
#     'xticklabels': lat_labels, 'yticks': lev_ticks, 'yticklabels':
#     False, 'do_colorbar': True, 'cbar_ticks': cbar_ticks
# }
# vert_kwargs_b = {
#     'do_ang_mom': True, 'do_had_max': True, 'min_cntr': minc,
#     'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lat_type': lat_type, 'lat_bounds': lat_bounds, 'lev_type':
#     lev_type, 'lev_bounds': lev_bounds, 'xticks': lat_ticks,
#     'xticklabels': lat_labels, 'yticks': lev_ticks, 'yticklabels':
#     '', 'do_colorbar': False, 'cbar_ticks': cbar_ticks
# }
# vert_kwargs_int = {
#     'do_ang_mom': True, 'do_had_max': True, 'min_cntr': minc,
#     'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lat_type': lat_type, 'lat_bounds': lat_bounds, 'lev_type':
#     lev_type, 'lev_bounds': lev_bounds, 'xticks': lat_ticks,
#     'xticklabels': None, 'yticks': lev_ticks, 'yticklabels': ' ',
#     'do_colorbar': False, 'cbar_ticks': cbar_ticks
# }
# # Zonal-mean annual cycle contour plots.
# ssnl_kwargs_l = {
#     'do_itcz': do_itcz, 'had_bnds_type': had_bnds_type,
#     'plot_had_bnds': plot_had_bnds, 'mask_hadley': mask_hadley,
#     'plot_thermal_eq': False, 'min_cntr': minc, 'max_cntr': maxc,
#     'num_cntr': nctr, 'col_map': col_map, 'lat_type': lat_type,
#     'lat_bounds': lat_bounds, 'xticks': range(1,13), 'xticklabels': '',
#     'yticks': lat_ticks, 'yticklabels': lat_labels, 'do_colorbar':
#     False, 'cbar_ticks': False, 'do_hatch': False
# }
# ssnl_kwargs_bl = {
#     'do_itcz': do_itcz, 'had_bnds_type': had_bnds_type,
#     'plot_had_bnds': plot_had_bnds, 'mask_hadley': mask_hadley,
#     'plot_thermal_eq': False, 'min_cntr': minc, 'max_cntr': maxc,
#     'num_cntr': nctr, 'col_map': col_map, 'lat_type': lat_type,
#     'lat_bounds': lat_bounds, 'xticks': range(1,13), 'xticklabels':
#     tuple('JFMAMJJASOND'), 'yticks': lat_ticks, 'yticklabels':
#     lat_labels, 'do_colorbar':'all', 'cbar_ticks': False, 'do_hatch':
#     False
# }
# ssnl_kwargs_b = {
#     'do_itcz': do_itcz, 'had_bnds_type': had_bnds_type,
#     'plot_had_bnds': plot_had_bnds, 'mask_hadley': mask_hadley,
#     'plot_thermal_eq': False, 'min_cntr': minc, 'max_cntr': maxc,
#     'num_cntr': nctr, 'col_map': col_map, 'lat_type': lat_type,
#     'lat_bounds': lat_bounds, 'xticks': range(1,13), 'xticklabels':
#     tuple('JFMAMJJASOND'), 'yticks': lat_ticks, 'yticklabels': ' ',
#     'do_colorbar':False, 'cbar_ticks': False, 'do_hatch': False
# }
# ssnl_kwargs_int = {
#     'do_itcz': do_itcz, 'had_bnds_type': had_bnds_type,
#     'plot_had_bnds': plot_had_bnds, 'mask_hadley': mask_hadley,
#     'plot_thermal_eq': False, 'min_cntr': minc, 'max_cntr': maxc,
#     'num_cntr': nctr, 'col_map': col_map, 'lat_type': lat_type,
#     'lat_bounds': lat_bounds, 'xticks': range(1,13), 'xticklabels': '',
#     'yticks': lat_ticks, 'yticklabels': ' ', 'do_colorbar':False,
#     'cbar_ticks': False, 'do_hatch': False
# }
# # 1-D annual cycle plots
# time_1d_kwargs_l = {
#     'xlim': time_lim, 'xticks': time_ticks, 'xticklabels': ' ', 'ylim': ylim,
#     'yticks': yticks, 'yticklabels':yticklabels, 'color': colors, 'linestyle':
#     linestyles, 'linewidth': linewidths
# }
# time_1d_kwargs_bl = {
#     'xlim': time_lim, 'xticks': time_ticks, 'xticklabels': time_ticklabels,
#     'ylim': ylim, 'yticks': yticks, 'yticklabels': yticklabels, 'color':
#     colors, 'linestyle': linestyles, 'linewidth': linewidths
# }
# time_1d_kwargs_b = {
#     'xlim': time_lim, 'xticks': time_ticks, 'xticklabels': time_ticklabels,
#     'ylim': ylim, 'ylabel': None, 'yticks': yticks, 'yticklabels': '', 'color':
#     colors, 'linestyle': linestyles, 'linewidth': linewidths
# }
# time_1d_kwargs_int = {
#     'xlim': time_lim, 'xticks': time_ticks, 'xticklabels': ' ', 'ylim': ylim,
#     'ylabel': None, 'yticks': yticks, 'yticklabels': ' ', 'color': colors,
#     'linestyle': linestyles, 'linewidth': linewidths 
# }
# # Vertical profile annual cycle contour plots.
# time_vert_kwargs_l = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map':
#     col_map, 'lev_type': lev_type, 'lev_bounds': lev_bounds, 'xticks':
#     range(1,13), 'xticklabels': ' ', 'yticks': lev_ticks, 'ylabel':
#     lev_label, 'yticklabels': lev_ticklabel, 'do_colorbar': False,
#     'cbar_label': cbar_label, 'cbar_ticks': cbar_ticks,
#     'cbar_ticklabels': cbar_ticklabels, 'cbar_extend': cbar_extend,
#     'do_hatch': False, 'do_cont_cntrs': do_cont_cntrs,
#     'cont_cntr_levs': cont_cntr_levs, 'cont_cntr_labels':
#     cont_cntr_labels, 'cont_run': cont_run,
#     'do_stat_mask': do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs,
#     'stat_mask_hatch': stat_mask_hatch, 'stat_mask_color':
#     stat_mask_color, 'do_norm_cont': do_norm_cont, 'do_normalize': do_normalize,
#     'norm_region': norm_region, 'norm_run': norm_run, 'norm_var': norm_var
# }
# time_vert_kwargs_bl = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map':
#     col_map, 'lev_type': lev_type, 'lev_bounds': lev_bounds, 'xticks':
#     range(1,13), 'xticklabels': tuple('JFMAMJJASOND'), 'yticks':
#     lev_ticks, 'ylabel': lev_label, 'yticklabels': lev_ticklabel,
#     'do_colorbar': 'all', 'cbar_label': cbar_label, 'cbar_ticks':
#     cbar_ticks, 'cbar_ticklabels': cbar_ticklabels, 'cbar_extend':
#     cbar_extend, 'do_hatch': False, 'do_cont_cntrs': do_cont_cntrs,
#     'cont_cntr_levs': cont_cntr_levs, 'cont_cntr_labels':
#     cont_cntr_labels, 'cont_run': cont_run,
#     'do_stat_mask': do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs,
#     'stat_mask_hatch': stat_mask_hatch, 'stat_mask_color':
#     stat_mask_color, 'do_norm_cont': do_norm_cont, 'do_normalize': do_normalize,
#     'norm_region': norm_region, 'norm_run': norm_run, 'norm_var': norm_var
# }
# time_vert_kwargs_b = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map':
#     col_map, 'lev_type': lev_type, 'lev_bounds': lev_bounds, 'xticks':
#     range(1,13), 'xticklabels': tuple('JFMAMJJASOND'), 'yticks':
#     lev_ticks, 'ylabel': ' ', 'yticklabels': ' ', 'do_colorbar':
#     False, 'cbar_label': cbar_label, 'cbar_ticks': cbar_ticks,
#     'cbar_ticklabels': cbar_ticklabels, 'cbar_extend': cbar_extend,
#     'do_hatch': False, 'do_cont_cntrs': do_cont_cntrs,
#     'cont_cntr_levs': cont_cntr_levs, 'cont_cntr_labels':
#     cont_cntr_labels, 'cont_run': cont_run,
#     'do_stat_mask': do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs,
#     'stat_mask_hatch': stat_mask_hatch, 'stat_mask_color':
#     stat_mask_color, 'do_norm_cont': do_norm_cont, 'do_normalize': do_normalize,
#     'norm_region': norm_region, 'norm_run': norm_run, 'norm_var': norm_var
# }
# time_vert_kwargs_int = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map':
#     col_map, 'lev_type': lev_type, 'lev_bounds': lev_bounds, 'xticks':
#     range(1,13), 'xticklabels': ' ', 'yticks': lev_ticks, 'ylabel':
#     '', 'yticklabels': ' ', 'do_colorbar': False, 'cbar_label':
#     cbar_label, 'cbar_ticks': cbar_ticks, 'cbar_ticklabels':
#     cbar_ticklabels, 'cbar_extend': cbar_extend, 'do_hatch': False,
#     'do_cont_cntrs': do_cont_cntrs, 'cont_cntr_levs': cont_cntr_levs,
#     'cont_cntr_labels': cont_cntr_labels, 'cont_run':
#     cont_run, 'do_stat_mask': do_stat_mask, 'stat_mask_cntrs':
#     stat_mask_cntrs, 'stat_mask_hatch': stat_mask_hatch,
#     'stat_mask_color': stat_mask_color, 'do_norm_cont': do_norm_cont,
#     'do_normalize': do_normalize, 'norm_region': norm_region,
#     'norm_run': norm_run, 'norm_var': norm_var
# }
# # Longitude-vertical plots.
# lon_vert_kwargs_l = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lon_bounds': lon_bounds, 'lev_type': lev_type, 'lev_bounds': lev_bounds,
#     'xticks': lon_ticks, 'xticklabels': ' ', 'yticks': lev_ticks,
#     'yticklabels': lev_ticklabel, 'ylabel': lev_label, 'do_colorbar': False,
#     'cbar_ticks': cbar_ticks, 'cbar_ticklabels': cbar_ticklabels,
#     'cbar_extend': cbar_extend, 'lat_avg_range': lat_avg_range, 'do_hatch':
#     False, 'do_cont_cntrs': do_cont_cntrs, 'cont_cntr_levs': cont_cntr_levs,
#     'cont_cntr_labels': cont_cntr_labels, 'cont_run': cont_run,
#     'do_stat_mask': do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs,
#     'stat_mask_hatch': stat_mask_hatch, 'stat_mask_color': stat_mask_color,
#     'do_norm_cont': do_norm_cont, 'do_normalize': do_normalize, 'norm_region':
#     norm_region, 'norm_run': norm_run, 'norm_var': norm_var
# }
# lon_vert_kwargs_bl = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lon_bounds': lon_bounds, 'lev_type': lev_type, 'lev_bounds': lev_bounds,
#     'xticks': lon_ticks, 'xticklabels': lon_labels, 'yticks': lev_ticks,
#     'yticklabels': lev_ticklabel, 'ylabel': lev_label, 'do_colorbar': 'all',
#     'cbar_ticks': cbar_ticks, 'cbar_ticklabels': cbar_ticklabels,
#     'cbar_extend': cbar_extend, 'lat_avg_range': lat_avg_range, 'do_hatch':
#     False, 'do_cont_cntrs': do_cont_cntrs, 'cont_cntr_levs': cont_cntr_levs,
#     'cont_cntr_labels': cont_cntr_labels, 'cont_run': cont_run,
#     'do_stat_mask': do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs,
#     'stat_mask_hatch': stat_mask_hatch, 'stat_mask_color': stat_mask_color,
#     'do_norm_cont': do_norm_cont, 'do_normalize': do_normalize, 'norm_region':
#     norm_region, 'norm_run': norm_run, 'norm_var': norm_var
# }
# lon_vert_kwargs_b = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lon_bounds': lon_bounds, 'lev_type': lev_type, 'lev_bounds': lev_bounds,
#     'xticks': lon_ticks, 'xticklabels': lon_labels, 'yticks': lev_ticks,
#     'yticklabels': '', 'ylabel': False, 'do_colorbar': False, 'cbar_ticks':
#     cbar_ticks, 'cbar_ticklabels': cbar_ticklabels, 'cbar_extend': cbar_extend,
#     'lat_avg_range': lat_avg_range, 'do_hatch': False, 'do_cont_cntrs':
#     do_cont_cntrs, 'cont_cntr_levs': cont_cntr_levs, 'cont_cntr_labels':
#     cont_cntr_labels, 'cont_run': cont_run, 'do_stat_mask':
#     do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs, 'stat_mask_hatch':
#     stat_mask_hatch, 'stat_mask_color': stat_mask_color, 'do_norm_cont':
#     do_norm_cont, 'do_normalize': do_normalize, 'norm_region': norm_region,
#     'norm_run': norm_run, 'norm_var': norm_var
# }
# lon_vert_kwargs_int = {
#     'min_cntr': minc, 'max_cntr': maxc, 'num_cntr': nctr, 'col_map': col_map,
#     'lon_bounds': lon_bounds, 'lev_type': lev_type, 'lev_bounds': lev_bounds,
#     'xticks': lon_ticks, 'xticklabels': ' ', 'yticks': lev_ticks,
#     'yticklabels': ' ', 'ylabel': False, 'do_colorbar': False, 'cbar_ticks':
#     cbar_ticks, 'cbar_ticklabels': cbar_ticklabels, 'cbar_extend': cbar_extend,
#     'lat_avg_range': lat_avg_range, 'do_hatch': False, 'do_cont_cntrs':
#     do_cont_cntrs, 'cont_cntr_levs': cont_cntr_levs, 'cont_cntr_labels':
#     cont_cntr_labels, 'cont_run': cont_run, 'do_stat_mask':
#     do_stat_mask, 'stat_mask_cntrs': stat_mask_cntrs, 'stat_mask_hatch':
#     stat_mask_hatch, 'stat_mask_color': stat_mask_color, 'do_norm_cont':
#     do_norm_cont, 'do_normalize': do_normalize, 'norm_region': norm_region,
#     'norm_run': norm_run, 'norm_var': norm_var
# }

                # Get panel-specific keyword arguments.
                # tmp_kwargs = panel_kwargs.copy()
                # for key, val in panel_kwargs.iteritems():
                    # if type(val) in [tuple, list] and type(val[0]) in [tuple, list]:
                        # tmp_kwargs.update({key: val[n][m]})
                # Plot each data.
                # if self.plot_type[n] == 'lat_1d':
                #     self.plot_handles.append(plot_lat_1d(
                #         self.ax[n], proj[n][m], models[n][m], runs[n][m], ens_mems[n][m],
                #         variables[n][m], x_data[n][m], level[n][m], intvls[n][m],
                #         yr_ranges[n][m], **tmp_kwargs
                #     ))
                # elif self.plot_type[n] == 'lat_vert':
                #     self.plot_handles.append(plot_lat_vert(
                #         self.ax[n], proj[n][m], models[n][m], runs[n][m], ens_mems[n][m],
                #         variables[n][m], level[n][m], intvls[n][m], yr_ranges[n][m],
                #         **tmp_kwargs
                #     ))
                # elif self.plot_type[n] == 'lon_vert':
                #     self.plot_handles.append(plot_lon_vert(
                #         self.ax[n], proj[n][m], models[n][m], runs[n][m], ens_mems[n][m],
                #         variables[n][m], intvls[n][m], yr_ranges[n][m],
                #         x_data[n][m], z_data[n][m], **tmp_kwargs
                #     ))
                # elif self.plot_type[n] in ('lon_lat', 'map'):
                #      self.plot_handles.append(plot_lon_lat(
                #          self.ax[n], proj[n][m], models[n][m], runs[n][m], ens_mems[n][m],
                #          variables[n][m], level[n][m], intvls[n][m], dtype[n][m],
                #          yr_ranges[n][m], **tmp_kwargs
                #      ))
                # elif self.plot_type[n] == 'ssnl_lat':
                #     self.plot_handles.append(plot_ssnl_lat(
                #         self.ax[n], proj[n][m], models[n][m], runs[n][m], ens_mems[n][m],
                #         variables[n][m], level[n][m], yr_ranges[n][m], **tmp_kwargs
                #     ))
                # elif self.plot_type[n] == 'time_1d':
                #     self.plot_handles.append(plot_time_1d(
                #         self.ax[n], proj[n][m], models[n][m], runs[n][m], regions[n][m],
                #         ens_mems[n][m], variables[n][m], level[n][m], intvls[n][m],
                #         yr_ranges[n][m], **tmp_kwargs
                #     ))
                # elif self.plot_type[n] == 'time_vert':
                #     self.plot_handles.append(plot_time_vert(
                #         self.ax[n], proj[n][m], models[n][m], runs[n][m], regions[n][m],
                #         ens_mems[n][m], variables[n][m], level[n][m], yr_ranges[n][m],
                #         **tmp_kwargs
                #     ))
                # elif self.plot_type[n] == '1d_vert':
                #     self.plot_handles.append(plot_1d_vert(
                #         self.ax[n], proj[n][m], models[n][m], runs[n][m], regions[n][m],
                #         ens_mems[n][m], variables[n][m], z_data[n][m], level[n][m],
                #         intvls[n][m], yr_ranges[n][m], **tmp_kwargs
                #     ))
        # Axis legend.
        # if self.do_legend and '1d' in self.plot_type[n]:
            # self.plot_handles = [ph[0] for ph in self.plot_handles]
            # if (self.do_legend == 'first' and n == 0) or (self.do_legend is True):
                # if legend_labels == 'models':
                    # legend_labels = models[n]
                # leg = self.ax[n].legend(self.plot_handles, legend_labels,
                                   # frameon=False, loc=0, fontsize='x-small')
        # Axis vertical labels.
        # if self.ax_left_labels:
        #     if self.ax_left_label_rot == 'horizontal':
        #         horiz_frac = -0.17
        #     else:
        #         horiz_frac = -0.05
        #     self.ax[n].text(
        #         horiz_frac, 0.5, self.ax_left_labels[n], verticalalignment='center',
        #         horizontalalignment='left',
        #         rotation=self.ax_left_label_rot, fontsize='small', 
        #         transform=self.ax[n].transAxes
        #     )
        # if self.ax_right_labels:
        #     horiz_frac = 1.03
        #     self.ax[n].text(
        #         horiz_frac, 0.5, self.ax_right_labels[n], verticalalignment='center',
        #         rotation='vertical', fontsize='small', transform=self.ax[n].transAxes
        #     )
        # Update colorbar axes limits.
        # if panel_kwargs.get('do_colorbar') == 'all':
            # panel_kwargs.update({'ax_cbar': self.ax_cbar})

    ## Contour plot settings
    # minc=0; maxc=1; nctr=10; cbar_ticks=False,
    # minc=0; maxc=1.2; nctr=12; cbar_ticks=False,
    # minc=-0.005; maxc=.145; nctr=15; cbar_ticks=False,
    # minc=0; maxc=2; nctr=16; cbar_ticks=False,
    # minc=0; maxc=3; nctr=12; cbar_ticks=False,
    # minc=0; maxc=5; nctr=10; cbar_ticks=False,
    # minc=0; maxc=15; nctr=10; cbar_ticks=False,
    # minc=0; maxc=18; nctr=18; cbar_ticks=False,
    # minc=0; maxc=20; nctr=10; cbar_ticks=False,
    # minc=0; maxc=25; nctr=10; cbar_ticks=False,
    # minc=0; maxc=30; nctr=10; cbar_ticks=False,
    # minc=0; maxc=80; nctr=10; cbar_ticks=False,
    # minc=0; maxc=100; nctr=10; cbar_ticks=False,
    # minc=210; maxc=310; nctr=10; cbar_ticks=False,
    # minc=-0.02; maxc=0.02; nctr=20; cbar_ticks=False,
    # minc=-0.2; maxc=0.2; nctr=16; cbar_ticks=False,
    # minc=-0.3; maxc=0.3; nctr=12; cbar_ticks=False,
    # minc=-0.5; maxc=0.5; nctr=10; cbar_ticks=False,
    # minc=-0.75; maxc=0.75; nctr=15; cbar_ticks=np.arange(-.75, .76, .25),
    # minc=-1; maxc=1; nctr=20; cbar_ticks=False,
    # minc=-1.4; maxc=1.4; nctr=14; cbar_ticks=False,
    # minc=-1.5; maxc=1.5; nctr=10; cbar_ticks=False,
    # minc=-2; maxc=2; nctr=16; cbar_ticks=False,
    # minc=-2.5; maxc=2.5; nctr=20; cbar_ticks=range(-2,3),
    # minc=-3; maxc=3; nctr=12; cbar_ticks=False,
    # minc=-3.5; maxc=3.5; nctr=14; cbar_ticks=False,
    # minc=-4; maxc=4; nctr=16; cbar_ticks=range(-4,5,2),
    # minc=-4.5; maxc=4.5; nctr=18; cbar_ticks=False,
    # minc=-5; maxc=5; nctr=10; cbar_ticks=range(-4,5,2),
    # minc=-6; maxc=6; nctr=12; cbar_ticks=range(-6,7,2),
    # minc=-6.5; maxc=6.5; nctr=13; cbar_ticks=range(-6,7,2),
    # minc=-7; maxc=7; nctr=14; cbar_ticks=range(-6,7,2),
    # minc=-7.5; maxc=7.5; nctr=15; cbar_ticks=range(-7,8,2),
    # minc=-8; maxc=8; nctr=16; cbar_ticks=range(-8,9,2),
    # minc=-10; maxc=10; nctr=10; cbar_ticks=False,
    # minc=-12; maxc=12; nctr=12; cbar_ticks=False,
    # minc=-16; maxc=16; nctr=16; cbar_ticks=range(-16, 17, 4),
    # minc=-20; maxc=20; nctr=20; cbar_ticks=False,
    # minc=-25; maxc=25; nctr=20; cbar_ticks=range(-25, 26, 10),
    # minc=-30; maxc=30; nctr=15; cbar_ticks=range(-30, 31, 12),
    # minc=-35; maxc=35; nctr=14; cbar_ticks=False,
    # minc=-40; maxc=40; nctr=16; cbar_ticks=False,
    # minc=-50; maxc=50; nctr=20; cbar_ticks=False,
