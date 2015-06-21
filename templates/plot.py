#! /usr/bin/env python
import numpy as np
import aospy.plotting

def plot():
    fig = aospy.plotting.Fig(
        n_row=1,
        n_col=1,
        n_ax='all',
        n_plot=1,
        n_data=2,

        row_size=4,
        col_size=5,
        subplot_lims={'left': 0.05, 'right': 0.95, 'wspace': 0.1,
                      'bottom': 0.12, 'top': 0.88, 'hspace': 0.15},

        min_cntr=-3.5,
        max_cntr=3.5,
        num_cntr=13,
        contourf_extend='both', # 'auto' 'neither' 'min' 'max' 'both'
        col_map='BrBG',
        do_colorbar='all',      # 'all' 'column' 'row' False True
        cbar_ax_lim = (0.1, 0.08, 0.8, 0.03),
        cbar_ticks=False,
        cbar_ticklabels=False,
        cbar_label='units',

        proj='aero_3agcm',
        model='am2',
        run='reyoi_cont',
        ens_mem=None,
        var=['precip', 'evap'],
        intvl_in='monthly',
        intvl_out='jas',
        dtype_in_time='ts',
        dtype_in_vert=False,
        dtype_out_time='reg.av',
        dtype_out_vert=False,
        level=None,
        region='sahel',
        yr_range='default',

        plot_type='scatter',
        x_dim=False,
        y_dim=False,

        ## Titles and labels
        fig_title=False,
        ax_title=False,
        ax_left_label=False,
        ax_right_label=False,

        # Axis limits, ticks, and labels
        x_lim=False,
        x_ticks=False,
        x_ticklabels=False,
        x_label=False,

        y_lim=False,
        y_ticks=False,
        y_ticklabels=False,
        y_label=False,

        lat_lim=(-45, 45),
        lat_ticks=False,
        lat_ticklabels=False,
        lat_label=False,

        lon_lim=(-180, 180),
        lon_ticks=False,
        lon_ticklabels=False,
        lon_label=False,

        p_lim=(1000, 100),
        p_ticks=False,
        p_ticklabels=False,
        p_label=False,
        
        sigma_lim=(1, 0.1),
        sigma_ticks=False,
        sigma_ticklabels=False,
        sigma_label=False,

        time_lim='ann_cycle',
        time_ticks=False,
        time_ticklabels=False,
        time_label=False,

        ## Map plot parameters.
        map_proj='cyl',
        map_res='c',
        left_lon=0,
        shiftgrid_start=False,
        shiftgrid_cyclic=360.0,
        latlon_rect=False,
        ## Quiver (i.e. arrow) plot overlayed on maps.
        do_quiver=False,
        ## Line plot parameters.
        line_color='k',
        linestyle='-',
        ## Scatter plot parameters
        marker_size=10,
        marker_color='k',
        marker_shape='.',
        ## Transformations to apply to data.
        do_subtract_mean=False,
    )

    fig.create_fig()
    fig.make_plots()
    plt.show()
    return fig
    
if __name__ == '__main__':
    fig = plot_mult()
