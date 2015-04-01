import cPickle, os
from scipy.stats import scoreatpercentile, t
import numpy as np
from numpy import ma
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset, MFDataset
from aospy.av_stat import grid_sfc_area
from aospy.data_io import _time_label, _ens_label, _yr_label, run_label
from aospy.plotting import subcolorbar, trunc_col_map

def make_plots(proj, model, run, ens_mem, var, level, intvl,
               intvl_type, data_type, yr_range):
    """Make and save plots of desired variables/calculations."""
    # Create string label for file names of start and end years.
    time_label, intvl = _time_label(intvl)
    start_yr, end_yr = yr_range
    yr_label = _yr_label(start_yr, end_yr)
    num_yr = end_yr - start_yr + 1
    ens_label = _ens_label(ens_mem)
    var_label = _var_label(var, level)
    # # Set statistics parameters.
    # conf_thresh = 0.95
    # interv = t.interval(conf_thresh, num_yr - 1)[1]
    # Create output label based on time interval.
    path_out = (proj.direc_out + '/figs/' + model.name + '/' +
                run.name + '/' + ens_label + '/')
    suffix = (model.name + '.' + run.name + '.' + ens_label + '.' + time_label +
              '.' + yr_label + '.p')
    # Tags used to ID different data types in their file names.
    tags = ('av', 'av.znl', 'av.z_asym', 'ts.znl')
    # Upload each variable.
    for tag in tags:
        name = var_label + '.' + tag + '.' + suffix
        try:
            data_in = open(path_in + name, 'r')
            data = cPickle.load(data_in)
        except:
            print 'no data'
        if False:
            path_cont = pre_in + '/data/' + model + '/cont/'
            name_cont = (var_label + '.' + tag+ '.' + model +
                         '.cont.' + time_label + '.' +
                         yr_label + '.p')
            file_cont = path_cont + name_cont
            data_in_cont = open(file_cont, 'r')
            data_cont = cPickle.load(data_in_cont)
            data_in_cont.close()
        # Plot.
        fig = plt.figure(figsize=(5,3))
        ax = fig.add_subplot(111)
        plot_title = (model + ' ' + run + ' ' + var_label +
                      ' ' + tag + ' ' + time_label + ' ' +
                      yr_label)
        plt.title(plot_title, fontsize=11)
        # For sin(latitude) plots.
        sin_lat = np.sin(np.deg2rad(lats))
        sin_range = np.arange(-90., 91., 10.)
        sin_ticks = np.sin(np.deg2rad(sin_range))
        sin_lbl = ('90S', '', '', '60S', '', '', '30S', '', '',
                   'EQ', '', '','30N', '', '', '60N', '', '', '90N')
        # Plot type depends on data.
        name_std = var_label + name_tag + suffix
        file_in = path_in + name_std
        file_in = open(file_in, 'r')
        stdev = cPickle.load(file_in)
        conf = stdev*interv
        file_in.close()
        ax.fill_between(lats, data+conf, data-conf, alpha=0.5)
        ax.plot(lats, data, 'k')
        ax.fill_between(sin_lat, data+conf, data-conf, alpha=0.5)
        ax.plot(sin_lat, data, 'k')
        plt.xlim(-1, 1)
        plt.xticks(sin_ticks, sin_lbl, fontsize=10)
        ax.axhline(0., color='grey')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_xlim((-90, 90))
        ax.set_xticks(range(-90, 91, 10))
        plt.xticks(fontsize=10)
        ax.set_xlim((-1., 1.))
        ax.set_xticks(sin_ticks)
        ax.set_xticklabels(sin_lbl, fontsize=10)
        plt.yticks(fontsize=10)
        if tag in ['av', 'av.z_asym']:
            # Contour intervals.
            n_cntr = 21
            thresh = 1.
            if is_diff:
                cmap, cntr_lvls = trunc_col_map(data, thresh, n_cntr, col_map)
            else:
                cmap = col_map
                cntr_lvls = np.linspace(data.min(), data.max(), n_cntr)
            # Map and contours.
            center_lon = -160
            data, lons_plot = shiftgrid(center_lon+181, data,
                                        lons, start=False)
            map = Basemap(projection='moll', lon_0=center_lon,
                          resolution='c', ax=ax)
            lon, lat = np.meshgrid(lons_plot, lats)
            X, Y = map(lon, lat)
            cs1 = map.contourf(X, Y, data, cntr_lvls, extend='both')
            cs1.set_cmap(cmap)
            cs2 = map.contour(X, Y, data, cntr_lvls, extend='both')
            cs2.set_cmap(cmap)
            if is_diff:
                data_cont, lons_plot = shiftgrid(center_lon+181, data_cont, 
                                                 lons, start=False)
                cs = map.contour(X, Y, data_cont, 10,
                                 colors='gray', linewidths=0.2)
            plt.clabel(cs, fontsize=6, fmt='%1.1f',
                       inline=True, inline_spacing=0)
            map.drawcoastlines(linewidth=0.1)
            map.drawmapboundary(linewidth=0.5)
            cbar = fig.colorbar(cs1, orientation='horizontal', 
                                shrink=0.8, pad=0.1)
            fig.text(0.19, 0.19, '%.1f' % data.min(),
                     horizontalalignment='right')
            fig.text(0.83, 0.19, '%.1f' % data.max(),
                     horizontalalignment='left')

        elif tag == 'ts.znl':
            time = range(data.shape[0])
            X, Y = np.meshgrid(time, sin_lat)
        # Contour intervals: center colorbar on zero, span all data.
            n_cntr = 21
            dat_max = ma.max(data)
            dat_min = ma.min(data)
            thresh = 1.
            cmap, cntr_lvls = trunc_col_map(data, thresh,n_cntr, col_map)
            cax = plt.contourf(X, Y, data.T, cntr_lvls, extend='both')
            cax.set_cmap(cmap)
            cax2 = plt.contour(X, Y, data.T, cntr_lvls, extend='both')
            cax2.set_cmap(cmap)
            cbar = fig.colorbar(cax, orientation='horizontal',
                                shrink=0.8, pad=0.1, format='%d')
            fig.text(0.19, 0.19, '%d' % dat_min,
                     horizontalalignment='right')
            fig.text(0.83, 0.19, '%d' % dat_max,
                     horizontalalignment='left')
#           if is_diff:
#               cs = plt.contour(X, Y, data_cont.T, colors='gray',
#                                linewidths=0.1)
#               plt.clabel(cs, fontsize=6, fmt='%0.1f',
#                          inline=True, inline_spacing=0)
            plt.ylim(-1, 1)
            plt.yticks(sin_ticks, sin_lbl, fontsize=9)
            plt.yticks(fontsize=9)    
        else:
            Y, Z = np.meshgrid(sin_lat, levs)
            if data.ndim == 3:
                data = data.mean(axis=-1)
                # Contour intervals: center colorbar on zero; span all data.
                n_cntr = 21
                dat_max = ma.max(data)
                dat_min = ma.min(data)
                thresh = 1.
                cmap, cntr_lvls = trunc_col_map(data, thresh,
                                                n_cntr, col_map)
                cax = plt.contourf(Y, Z, data, cntr_lvls,
                                   extend='both')
                cax.set_cmap(cmap)
                cax2 = plt.contour(Y, Z, data, cntr_lvls,
                                   extend='both')
                cax2.set_cmap(cmap)
                cbar = fig.colorbar(cax, orientation='horizontal',
                                    shrink=0.8, pad=0.1, format='%d')
                fig.text(0.19, 0.19, '%d' % dat_min,
                         horizontalalignment='right')
                fig.text(0.83, 0.19, '%d' % dat_max,
                         horizontalalignment='left')
                if is_diff:
                    if data_cont.ndim == 3:
                    data_cont = data_cont.mean(axis=-1)
                    cs = plt.contour(Y, Z, data_cont, 9,
                                     colors='gray', linewidths=0.1)
                    plt.clabel(cs, fontsize=6, fmt='%0.1f',
                               inline=True, inline_spacing=0)
                    plt.ylim(1000, 10)
                    plt.ylabel('Pressure (hPa)', fontsize=10)
                    plt.xlim(-1, 1)
                    plt.xticks(sin_ticks, sin_lbl, fontsize=9)
                    plt.yticks(fontsize=9)
        # Save plot.
        path_out = (pre_out + '/figs/' + model + '/' +
                    run + '/')
        name = (var_label + '.' + tag + '.' + model + '.' +
                run + '.' + time_label + '.' + yr_label)
        file_out = path_out + name
        fig.savefig(file_out, format='pdf')
        plt.close(fig)
        print name
        else:
            print ('Warning: ' + file_in + ' does not exist.' +
                   ' Skipping this plot.')
