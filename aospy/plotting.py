"""Classes for creating multi-panel figures using data generated via aospy."""
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap
import xarray as xr

from .__config__ import default_colormap, PFULL_STR
from .calc import Calc, CalcInterface
from .io import to_dup_list
from .operator import Operator
from .utils import to_hpa

fig_specs = (
    'fig_title', 'n_row', 'n_col', 'row_size', 'col_size', 'n_ax',
    'subplot_lims', 'cbar_ax_lim', 'cbar_ticks',
    'cbar_ticklabels', 'cbar_label', 'cbar_label_kwargs',
    'cbar_left_label', 'cbar_left_label_coords', 'cbar_left_label_kwargs',
    'cbar_right_label', 'cbar_right_label_coords', 'cbar_right_label_kwargs',
    'verbose'
)
ax_specs = (
    'n_plot', 'ax_title', 'ax_label', 'ax_label_coords',
    'ax_left_label', 'ax_left_label_coords', 'ax_left_label_kwargs',
    'ax_right_label', 'ax_right_label_coords', 'ax_right_label_kwargs',
    'map_proj', 'map_corners',
    'map_res', 'shiftgrid_start', 'shiftgrid_cyclic', 'do_legend',
    'legend_labels', 'legend_loc', 'x_dim', 'x_lim', 'x_ticks', 'x_ticklabels',
    'x_label', 'y_dim', 'y_lim', 'y_ticks', 'y_ticklabels', 'y_label',
    'lat_lim', 'lat_ticks', 'lat_ticklabels', 'lat_label', 'lon_lim',
    'lon_ticks', 'lon_ticklabels', 'lon_label', 'p_lim', 'p_ticks',
    'p_ticklabels', 'p_label', 'sigma_lim', 'sigma_ticks', 'sigma_ticklabels',
    'sigma_label', 'time_lim', 'time_ticks', 'time_ticklabels', 'time_label',
    'do_mark_x0', 'do_mark_y0'
)
plot_specs = (
    'plot_type', 'do_best_fit_line', 'print_best_fit_slope',
    'print_corr_coeff', 'cntr_lvls', 'colormap', 'min_cntr', 'max_cntr',
    'num_cntr', 'contours_extend', 'latlon_rect', 'do_mask_oceans',
    'contour_labels', 'contour_kwargs', 'contourf_kwargs', 'plot_kwargs',
    'quiver_kwargs', 'quiver_n_lon', 'quiver_n_lat', 'do_quiverkey',
    'quiverkey_args', 'quiverkey_kwargs', 'scatter_kwargs', 'do_colorbar'
)
data_specs = (
    'proj', 'model', 'run', 'ens_mem', 'var', 'level', 'region', 'date_range',
    'intvl_in', 'intvl_out', 'dtype_in_time', 'dtype_in_vert',
    'dtype_out_time', 'dtype_out_vert', 'do_subtract_mean', 'do_mult_factor',
    'mask_unphysical'
)
specs = fig_specs + ax_specs + plot_specs + data_specs


class Fig(object):
    """Class for producing figures with one or more panels."""
    def __init__(self, fig_params, n_ax=1, n_plot=1, n_data=1, n_row=1,
                 n_col=1, date_range=None, intvl_in=None, intvl_out=None,
                 dtype_in_time=None, dtype_in_vert=None, dtype_out_time=None,
                 dtype_out_vert=None, level=None, **kwargs):
        self.__dict__ = vars(fig_params)

        self.n_ax = n_ax
        self.n_plot = n_plot
        self.n_data = n_data
        self.n_row = n_row
        self.n_col = n_col
        self._set_n_ax_plot_data()

        self.date_range = date_range
        self.intvl_in = intvl_in
        self.intvl_out = intvl_out
        self.dtype_in_time = dtype_in_time
        self.dtype_in_vert = dtype_in_vert
        self.dtype_out_time = dtype_out_time
        self.dtype_out_vert = dtype_out_vert
        self.level = level

        self.ax = []

        # Accept all other keyword arguments passed in as attrs.
        for key, val in kwargs.items():
            setattr(self, key, val)

        self._expand_attrs_over_tree()
        self._make_ax_objs()

    def _set_n_ax_plot_data(self):
        """Set the number of axes, plots, and data."""
        # Set the number of Axes.
        if self.n_ax == 'all':
            self.n_ax = self.n_row*self.n_col
        else:
            assert self.n_ax <= self.n_row*self.n_col
        # Distribute n_plot across the Axes.
        self.n_plot = to_dup_list(self.n_plot, self.n_ax)
        # Distribute n_data across the Plots.
        self.n_data = to_dup_list(self.n_data, self.n_ax)
        for i, ndt in enumerate(self.n_data):
            self.n_data[i] = to_dup_list(ndt, self.n_plot[i])

    def _traverse_child_tree(self, value, level='data'):
        """
        Traverse the "tree" of child Ax, Plot, and Data objects by creating a
        nested list of three levels, each level cascading one child object
        class down (e.g. from Ax to Plot to Data).
        """
        assert level in ('fig', 'ax', 'plot', 'data')
        if level in ('ax', 'plot', 'data'):
            value = to_dup_list(value, self.n_ax)
            if level in ('plot', 'data'):
                for i, vi in enumerate(value):
                    value[i] = to_dup_list(vi, self.n_plot[i])
                    if level == 'data':
                        for j, vij in enumerate(value[i]):
                            value[i][j] = to_dup_list(
                                vij, self.n_data[i][j], single_to_list=False
                            )
        return value

    def _expand_specs(self, specs, spec_name):
        for attr in specs:
            value = getattr(self, attr, False)
            setattr(self, attr, self._traverse_child_tree(value, spec_name))

    def _expand_attrs_over_tree(self):
        """
        Replicate attrs such that each Ax-level attr is a list with length
        equal to the number of Axes.  Similarly, Plot-level attrs become
        two-level lists, the first level being corresponding to the Axes and
        the second level to the Plots within each Ax.  And likewise for
        Data-level attrs.
        """
        for specs, name in ((fig_specs, 'fig'),
                            (ax_specs, 'ax'),
                            (plot_specs, 'plot'),
                            (data_specs, 'data')):
            self._expand_specs(specs, name)

    def _locate_ax(self, panel_num):
        """Determine if the panel is interior, left, bottom-left, or bottom."""
        if panel_num == 0 or panel_num == self.n_col*(self.n_row - 1):
            return 'bottomleft'
        elif panel_num > self.n_col*(self.n_row - 1):
            return 'bottom'
        elif not panel_num % self.n_col:
            return 'left'
        else:
            return 'interior'

    def _make_ax_objs(self):
        """Create the Ax obj for each panel."""
        self.ax = [Ax(self, n, self._locate_ax(n))
                   for n in range(self.n_ax)]

    @staticmethod
    def __add_text(ax, coords, string, kwargs):
        ax.text(coords[0], coords[1], string, **kwargs)

    def _make_colorbar(self, ax):
        """Create colorbar for multi panel plots."""
        # Don't make if already made.
        if hasattr(self, 'cbar'):
            return
        # Goes at bottom center if for all panels.
        self.cbar_ax = self.fig.add_axes(self.cbar_ax_lim)
        self.cbar = self.fig.colorbar(
            ax.Plot[0].handle, cax=self.cbar_ax, orientation='horizontal',
            drawedges=False, spacing='proportional',
            extend=self.contours_extend
        )
        # Set tick properties.
        if np.any(self.cbar_ticks):
            self.cbar.set_ticks(self.cbar_ticks)
        if self.cbar_ticklabels not in (None, False):
            self.cbar.set_ticklabels(self.cbar_ticklabels)
        self.cbar.ax.tick_params(labelsize='x-small')
        # Add center, left, and right labels as desired.
        if self.cbar_label:
            var = self.var[0][0][0]
            if self.cbar_label == 'units':
                if self.dtype_out_vert[0][0][0] == 'vert_int':
                    label = var.units.vert_int_plot_units
                else:
                    label = var.units.plot_units
            else:
                label = self.cbar_label
            self.cbar.set_label(label, **self.cbar_label_kwargs)

        def make_cbar_label_args(obj, string):
            return [obj.cbar_ax] + [
                    getattr(obj, 'cbar_' + string + '_label' + suffix, False)
                    for suffix in ('_coords', '', '_kwargs')
            ]

        if self.cbar_left_label:
            self.__add_text(*make_cbar_label_args(self, 'left'))
        if self.cbar_right_label:
            self.__add_text(*make_cbar_label_args(self, 'right'))

    def create_fig(self):
        """Create the figure and set up the subplots."""
        self.fig = plt.figure(figsize=(self.n_col*self.col_size,
                                       self.n_row*self.row_size))
        self.fig.subplots_adjust(**self.subplot_lims)
        if self.fig_title:
            self.fig.suptitle(self.fig_title, fontsize=12)

        for n in range(self.n_ax):
            self.ax[n].ax = self.fig.add_subplot(self.n_row, self.n_col, n+1)

    def make_plots(self):
        """Render the plots in every Ax."""
        for n in range(self.n_ax):
            self.ax[n].make_plots()

    def savefig(self, *args, **kwargs):
        """Save the Fig using matplotlib's built-in 'savefig' method."""
        self.fig.savefig(*args, **kwargs)

    def draw(self):
        """Call the matplotlib method canvas.draw() to re-render the figure."""
        self.fig.canvas.draw()


class AxInterface(object):
    pass


class Ax(object):
    # Which labels to include based on position in the figure.
    labels = {
        'left': {'x_ticklabels': ' ', 'y_ticklabels': True,
                 'x_label': False, 'y_label': True, 'do_colorbar': False},
        'interior': {'x_ticklabels': ' ', 'y_ticklabels': ' ',
                     'x_label': False, 'y_label': False, 'do_colorbar': False},
        'bottomleft': {'x_ticklabels': True, 'y_ticklabels': True,
                       'x_label': True, 'y_label': True, 'do_colorbar': True},
        'bottom': {'x_ticklabels': True, 'y_ticklabels': ' ',
                   'x_label': True, 'y_label': False, 'do_colorbar': True}
    }

    def __init__(self, Fig, ax_num, ax_loc):
        self.fig = Fig
        self.ax_num = ax_num
        self.ax_loc = ax_loc
        self.n_plot = Fig.n_plot[ax_num]
        self.n_data = Fig.n_data[ax_num]
        self.Plot = []

        self._copy_attrs_from_fig()
        # self._set_ax_loc_specs()
        self._set_xy_attrs_to_coords()

    def _traverse_child_tree(self, value, level='data'):
        """Traverse the "tree" of child Plot, and Data objects."""
        assert level in ('ax', 'plot', 'data')
        if level in ('plot', 'data'):
            value = to_dup_list(value, self.n_plot)
            if level in ('data'):
                for i, vi in enumerate(value):
                    value[i] = to_dup_list(vi, self.n_data[i])
        return value

    def _copy_attrs_from_fig(self):
        """Copy the attrs of the parent Fig that correspond to this Ax."""
        for attr in ax_specs:
            value = getattr(self.fig, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'ax'))

        for attr in plot_specs:
            value = getattr(self.fig, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'plot'))

        for attr in data_specs:
            value = getattr(self.fig, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'data'))

    def _set_ax_loc_specs(self):
        """Set attrs that depend on Ax location within the Fig."""
        # Take the Fig's attr value if it's neeeded; otherwise set False.
        for key, val in self.labels[self.ax_loc].items():
            if val:
                if val == ' ':
                    new_val = ' '
                else:
                    new_val = getattr(self.fig, key, False)

            else:
                new_val = False
            setattr(self, key, new_val)

    def _set_xy_attrs_to_coords(self):
        """
        Set the x and y axis dimensions and related attributes to the values
        specified by the 'x_dim' and 'y_dim' attributes.  E.g. if self.x_dim =
        'lat', then set 'x_lim', etc. equal to 'lat_lim', etc.
        """
        for l, dim in zip(('x', 'y'), ('x_dim', 'y_dim')):
            prefix = getattr(self, dim)
            # prefix being False implies to use the actual x_lim, x_ticks, etc.
            if prefix is False:
                prefix = l
            for attr in ('lim', 'ticks', 'ticklabels', 'label'):
                setattr(self, '_'.join([l, attr]),
                        getattr(self, '_'.join([prefix, attr])))

    def _set_axes_props(self):
        """Set the properties of the matplotlib Axes instance."""
        if self.x_lim:
            if self.x_lim == 'ann_cycle':
                self.x_lim = (1, 12)
                self.x_ticks = range(1, 13)
                self.x_ticklabels = tuple('JFMAMJJASOND')
            self.ax.set_xlim(self.x_lim)
        if self.do_mark_y0:
            self.ax.axhline(color='0.5')
        if self.x_ticks:
            self.ax.set_xticks(self.x_ticks)
        if self.x_ticklabels:
            self.ax.set_xticklabels(self.x_ticklabels, fontsize='x-small')
        if self.x_label:
            self.ax.set_xlabel(self.x_label,
                               fontsize='x-small', labelpad=1)
        if self.y_lim:
            self.ax.set_ylim(self.y_lim)
        if self.do_mark_x0:
            self.ax.axvline(color='0.5')
        if self.y_ticks:
            self.ax.set_yticks(self.y_ticks)
        if self.y_ticklabels:
            self.ax.set_yticklabels(self.y_ticklabels, fontsize='x-small')
        if self.y_label:
            self.ax.set_ylabel(self.y_label, fontsize='x-small', labelpad=-2)

        self.ax.tick_params(labelsize='x-small')
        if not (self.x_dim == 'lon' and self.y_dim == 'lat'):
            self.ax.spines['right'].set_visible(False)
            self.ax.spines['top'].set_visible(False)
            self.ax.xaxis.set_ticks_position('bottom')
            self.ax.yaxis.set_ticks_position('left')

    def _set_axes_labels(self):
        """Create the axis title and other text labels."""
        # Axis title.
        if self.ax_title:
            self.ax.set_title(self.ax_title, fontsize='small')
        # Axis panel labels, i.e. (a), (b), (c), etc.
        if self.ax_label:
            if self.ax_label == 'auto':
                text = '({})'.format(tuple('abcdefghijklmnop')[self.ax_num])
            else:
                text = self.ax_label
            self.panel_label = self.ax.text(
                self.ax_label_coords[0], self.ax_label_coords[1],
                text, fontsize='small', transform=self.ax.transAxes
            )
        # Labels to left and/or right of Axis.
        if self.ax_left_label:
            # if self.ax_left_label_rot == 'horizontal':
                # horiz_frac = -0.17
            # else:
            self.ax.text(
                self.ax_left_label_coords[0], self.ax_left_label_coords[1],
                self.ax_left_label, transform=self.ax.transAxes,
                **self.ax_left_label_kwargs
            )
        if self.ax_right_label:
            self.ax.text(
                self.ax_right_label_coords[0], self.ax_right_label_coords[1],
                self.ax_right_label, transform=self.ax.transAxes,
                **self.ax_right_label_kwargs
            )

    def _make_plot_objs(self):
        """Create the Plot object for each plotted element."""
        self.Plot = []
        for n, plot_type in zip(range(self.n_plot), self.plot_type):
            plot_interface = PlotInterface(ax=self, plot_num=n,
                                           plot_type=plot_type)
            if plot_type == 'scatter':
                self.Plot.append(Scatter(plot_interface))
            elif plot_type == 'contour':
                self.Plot.append(Contour(plot_interface))
            elif plot_type == 'contourf':
                self.Plot.append(Contourf(plot_interface))
            elif plot_type == 'quiver':
                self.Plot.append(Quiver(plot_interface))
            elif plot_type == 'line':
                self.Plot.append(Line(plot_interface))
            else:
                raise TypeError("Plot type '%s' not recognized."
                                % plot_type)

    def make_plots(self):
        """Call the matplotlib plotting command for each Plot."""
        self._make_plot_objs()
        self._set_axes_props()
        self._set_axes_labels()

        for n in range(self.n_plot):
            self.Plot[n].plot()

        if self.do_legend:
            self.ax.legend(self.legend_labels, loc=self.legend_loc,
                           frameon=False, fontsize='small')


class PlotInterface(object):
    """Interface to Plot class."""
    def __init__(self, ax=False, plot_num=1, plot_type=False):
        self.ax = ax
        self.plot_num = plot_num
        self.plot_type = plot_type
        self._copy_attrs_from_ax(ax)

    def _copy_attrs_from_ax(self, ax):
        """Copy the attrs of the parent Ax corresponding to this Plot."""
        self.fig = ax.fig
        try:
            self.n_data = ax.n_data[self.plot_num]
        except AttributeError:
            self.n_data = 1

        for attr in plot_specs:
            value = getattr(ax, attr)[self.plot_num]
            tree_value = self._traverse_child_tree(value, 'plot')
            setattr(self, attr, tree_value)

        for attr in data_specs:
            value = getattr(ax, attr)[self.plot_num]
            tree_value = self._traverse_child_tree(value, 'data')
            setattr(self, attr, tree_value)

    def _traverse_child_tree(self, value, level='data'):
        """Traverse the "tree" of child Data objects."""
        assert level in ('plot', 'data')
        if level == 'data':
            value = to_dup_list(value, self.n_data)
        return value


class Plot(object):
    def __init__(self, plot_interface):
        """Class for plotting single data element."""
        self.__dict__ = vars(plot_interface)

        if isinstance(self.var[0], tuple):
            self.var = self.var[0]
            self.n_data = len(self.var)

        self.calc = self._make_calc_obj()
        # Pass the calcs to the _load_data method.
        self.data = [self._load_data(calc, n)
                     for n, calc in enumerate(self.calc)]
        # Strip extra dimensions as necessary.
        data_shape = np.shape(self.data)
        if data_shape[0] == 1:
            self.data = self.data[0]
        elif data_shape[0] == 2 and data_shape[1] == 1:
            self.data = np.squeeze(self.data)
        # _set_coord_arrays() below needs Calc, not Operator, objects.
        for i, calc in enumerate(self.calc):
            if isinstance(calc, Operator):
                self.calc[i] = calc.objects[0]

        # self._apply_data_transforms()

        self._set_coord_arrays()

        if self.ax.x_dim == 'lon' and self.ax.y_dim == 'lat':
            self.basemap = self._make_basemap()
            self.backend = self.basemap
        else:
            self.backend = self.ax.ax

    def plot(self):
        return NotImplementedError("plot() is an abstract method: it is only "
                                   "implemented for the classes that "
                                   "inherit from Plot.")

    def _make_calc_obj(self):
        calc_obj = []
        for i in range(self.n_data):
            if isinstance(self.run[i], Operator):
                calcs = [Calc(CalcInterface(
                    proj=self.proj[i], model=self.model[i],
                    run=run,
                    ens_mem=self.ens_mem[i], var=self.var[i],
                    date_range=self.date_range[i], region=self.region[i],
                    intvl_in=self.intvl_in[i], intvl_out=self.intvl_out[i],
                    dtype_in_time=self.dtype_in_time[i],
                    dtype_in_vert=self.dtype_in_vert[i],
                    dtype_out_time=self.dtype_out_time[i],
                    dtype_out_vert=self.dtype_out_vert[i], level=self.level[i],
                    verbose=self.fig.verbose
                )) for run in self.run[i].objects]
                calc_obj.append(Operator(self.run[i].operator, calcs))
            else:
                calc_obj.append(Calc(CalcInterface(
                    proj=self.proj[i], model=self.model[i], run=self.run[i],
                    ens_mem=self.ens_mem[i], var=self.var[i],
                    date_range=self.date_range[i], region=self.region[i],
                    intvl_in=self.intvl_in[i], intvl_out=self.intvl_out[i],
                    dtype_in_time=self.dtype_in_time[i],
                    dtype_in_vert=self.dtype_in_vert[i],
                    dtype_out_time=self.dtype_out_time[i],
                    dtype_out_vert=self.dtype_out_vert[i], level=self.level[i],
                    verbose=self.fig.verbose
                )))
        return calc_obj

    def _set_coord_arrays(self):
        """Set the arrays holding the x- and y-coordinates."""
        array_names = {'lat': 'lat', 'lon': 'lon', 'p': 'level',
                       'sigma': 'pfull', 'time': 'time', 'x': 'x', 'y': 'y'}
        if self.n_data == 1:
            mod_inds = (0, 0)
        else:
            mod_inds = (0, 1)
        for dim, data, lim, i in zip(('x_dim', 'y_dim'), ('x_data', 'y_data'),
                                     ('x_lim', 'y_lim'), mod_inds):
            array_key = getattr(self.ax, dim)

            if array_key in ('x', 'y'):
                if len(self.data) == 1:
                    array = self.data[0]
                else:
                    array = self.data
            else:
                try:
                    array = getattr(self.calc[i].model[0],
                                    array_names[array_key])
                except AttributeError:
                    array = getattr(self.calc[i][0].model[0],
                                    array_names[array_key])

            if array_key == 'p' and array is not None:
                array = to_hpa(array)

            if array_key == 'time' and lim == 'ann_cycle':
                array = np.arange(1, 13)

            setattr(self, data, array)

    def _make_basemap(self):
        if self.ax.x_lim:
            llcrnrlon, urcrnrlon = self.ax.x_lim
        else:
            llcrnrlon, urcrnrlon = -180, 180
        if self.ax.y_lim:
            llcrnrlat, urcrnrlat = self.ax.y_lim
        else:
            llcrnrlat, urcrnrlat = -90, 90

        return mpl_toolkits.basemap.Basemap(
            projection=self.ax.map_proj, resolution=self.ax.map_res,
            ax=self.ax.ax, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat
        )

    @staticmethod
    def regrid_to_avg_coords(dim, *arrs):
        """Average the coordinate arrays of two DataArrays or Dataset."""
        template = arrs[0][dim]
        avg = xr.DataArray(np.zeros(template.shape), dims=template.dims,
                           coords=template.coords)
        for arr in arrs:
            avg += arr[dim]
        avg /= len(arrs)
        for arr in arrs:
            arr[dim] = avg
        return arrs

    @classmethod
    def _perform_oper(cls, arr1, arr2, operator, region=False):
        if region:
            try:
                arr1, arr2 = cls.regrid_to_avg_coords(PFULL_STR, arr1, arr2)
            except KeyError:
                arr1, arr2 = cls.regrid_to_avg_coords(PFULL_STR + '_ref',
                                                      arr1, arr2)
        return eval('arr1' + operator + 'arr2')

    def _load_data(self, calc, n):
        if isinstance(calc, Operator):
            region = self.region[n]
            data = tuple(
                [cl.load(self.dtype_out_time[n],
                         dtype_out_vert=self.dtype_out_vert[n],
                         region=self.region[n], time=False, vert=self.level[n],
                         lat=False, lon=False, plot_units=True,
                         mask_unphysical=self.mask_unphysical)
                 for cl in calc.objects]
            )
            ans = self._perform_oper(data[0], data[1], calc.operator,
                                     region=region,)
            # Combine masks of the two inputs.
            try:
                joint_mask = (np.ma.mask_or(data[0].mask, data[1].mask),)
            except AttributeError:
                return ans
            else:
                return (np.ma.array(ans, mask=joint_mask),)

        if isinstance(calc, (list, tuple)):
            return tuple(
                [cl.load(dto, dtype_out_vert=dtv, region=reg, time=False,
                         vert=lev, lat=False, lon=False, plot_units=True,
                         mask_unphysical=self.mask_unphysical)
                 for cl, dto, dtv, reg, lev in zip(
                         calc, self.dtype_out_time, self.dtype_out_vert,
                         self.region, self.level
                 )]
            )
        return calc.load(self.dtype_out_time[n],
                         dtype_out_vert=self.dtype_out_vert[n],
                         region=self.region[n], time=False,
                         vert=self.level[n],
                         lat=False,
                         lon=False,
                         plot_units=True,
                         mask_unphysical=self.mask_unphysical)

    def _subtract_mean(self, data):
        return np.subtract(data, np.mean(data))

    def _apply_data_transforms(self):
        """Apply any specified transformations to the data once loaded."""
        transforms = {'do_subtract_mean': self._subtract_mean}
        for attr, method in transforms.items():
            for data, do_method in zip(['x_data', 'y_data'],
                                       getattr(self, attr)):
                if do_method:
                    setattr(self, data, method(getattr(self, data)))
        return data

    def prep_data_for_basemap(self):
        # if self.ax.shiftgrid_start:
        #     lon0 = 181.25
        # else:
        #     lon0 = self.ax.shiftgrid_start
        lon0 = 181.25
        self.lons = self.x_data
        self.lats = self.y_data

        self.plot_data = []

        if isinstance(self.data, xr.Dataset):
            loop_data = [self.data[self.var[0].name].values]
        else:
            loop_data = self.data
        for data in loop_data:
            pd, self.plot_lons = mpl_toolkits.basemap.shiftgrid(
                lon0, data, self.x_data,
                start=self.ax.shiftgrid_start, cyclic=self.ax.shiftgrid_cyclic
            )
            self.plot_data.append(pd)
        if len(self.plot_data) == 1:
            self.plot_data = self.plot_data[0]

        self.x_data, self.y_data = self.basemap(*np.meshgrid(self.plot_lons,
                                                             self.y_data))
        if self.do_mask_oceans:
            self.plot_data = mpl_toolkits.basemap.maskoceans(
                self.x_data, self.y_data, self.plot_data,
                inlands=False, resolution='c'
            )

    def plot_rectangle(self, x_min, x_max, y_min, y_max, **kwargs):
        xs = [x_min, x_max, x_max, x_min, x_min]
        ys = [y_min, y_min, y_max, y_max, y_min]
        return self.backend.plot(xs, ys, latlon=True, linewidth=0.8, **kwargs)

    def apply_basemap(self, basemap):
        """Apply basemap extras: coastlines, etc."""
        basemap.drawcoastlines(linewidth=0.1, color='k')
        basemap.drawmapboundary(linewidth=0.1, fill_color='0.85')
        if self.latlon_rect:
            self.plot_rectangle(*self.latlon_rect)

    def apply_colormap(self, colormap):
        if colormap == 'default':
            try:
                colormap = self.var[0].colormap
            except AttributeError:
                colormap = default_colormap
        self.handle.set_cmap(colormap)

    def corr_coeff(self, x_data, y_data, print_corr=True, print_pval=True):
        """Compute the Pearson correlation coefficient and plot it."""
        pearsonr, p_val = scipy.stats.pearsonr(x_data, y_data)
        if print_corr:
            self.ax.ax.text(0.3, 0.95, r'$r=$ %.2f' % pearsonr,
                            transform=self.ax.ax.transAxes, fontsize='x-small')
        if print_pval:
            self.ax.ax.text(0.3, 0.9, r'$p=$ %.3f' % p_val,
                            transform=self.ax.ax.transAxes, fontsize='x-small')
        return pearsonr, p_val

    def best_fit_line(self, print_slope=True, print_y0=True):
        """Plot the best fit line to the data."""
        # Enforce a dtype of float to ensure data plays nice with np.polyfit
        x_data = np.array(self.x_data).astype('float64')
        y_data = np.array(self.y_data).astype('float64')
        best_fit = np.polyfit(x_data, y_data, 1)
        x_lin_fit = [-1e3, 1e3]

        def lin_fit(m, x, b):
            return [m*xx + b for xx in x]
        self.backend.plot(x_lin_fit, lin_fit(best_fit[0], x_lin_fit,
                                             best_fit[1]), 'k')
        if print_slope:
            self.ax.ax.text(0.3, 0.1, r'slope = %0.2f' % best_fit[0],
                            transform=self.ax.ax.transAxes, fontsize='x-small')
        if print_y0:
            self.ax.ax.text(0.3, 0.05, r'y0 = %0.2f' % best_fit[1],
                            transform=self.ax.ax.transAxes, fontsize='x-small')
        return best_fit


class Contour(Plot):
    def __init__(self, plot_interface):
        """Contour plot."""
        Plot.__init__(self, plot_interface)
        self.plot_func = self.backend.contour
        self.plot_func_kwargs = self.contour_kwargs
        self.cntr_lvls = np.linspace(self.min_cntr, self.max_cntr,
                                     self.num_cntr + 1)

    def _prep_data(self):
        if self.basemap:
            self.prep_data_for_basemap()
        else:
            self.plot_data = self.data[0]
        if self.do_mult_factor[0]:
            self.plot_data = np.multiply(float(self.do_mult_factor[0]),
                                         self.plot_data)

    def plot(self):
        self._prep_data()
        self.handle = self.plot_func(self.x_data, self.y_data, self.plot_data,
                                     self.cntr_lvls, **self.plot_func_kwargs)
        if self.contour_labels:
            plt.gca().clabel(self.handle, fontsize=7, fmt='%1d')
        if self.colormap:
            self.apply_colormap(self.colormap)
        if self.do_colorbar:
            self.fig._make_colorbar(self.ax)
        if self.basemap:
            self.apply_basemap(self.basemap)


class Contourf(Contour):
    def __init__(self, plot_interface):
        """Filled contour ('contourf') plot."""
        Contour.__init__(self, plot_interface)
        self.plot_func = self.backend.contourf
        self.plot_func_kwargs = self.contourf_kwargs

class Line(Plot):
    def __init__(self, plot_interface):
        """Line plot."""
        Plot.__init__(self, plot_interface)

    def plot(self):
        self.handle = self.backend.plot(self.x_data, self.y_data,
                                        **self.plot_kwargs)


class Scatter(Plot):
    def __init__(self, plot_interface):
        Plot.__init__(self, plot_interface)
        self.x_data = np.squeeze(self.data[0])
        self.y_data = np.squeeze(self.data[1])
        self._apply_data_transforms()

    def plot(self):
        self.handle = self.backend.scatter(self.x_data, self.y_data,
                                           **self.scatter_kwargs)
        if self.do_best_fit_line:
            self.best_fit_line(print_slope=self.print_best_fit_slope)

        if self.print_corr_coeff:
            self.corr_coeff(self.x_data, self.y_data)


class Quiver(Plot):
    def __init__(self, plot_interface):
        """Quiver (i.e. vector) plot."""
        Plot.__init__(self, plot_interface)

    def prep_quiver(self):
        if not self.quiver_n_lon:
            self.quiver_n_lon = self.plot_lons.size
        if not self.quiver_n_lat:
            self.quiver_n_lat = self.lats.size
        return self.backend.transform_vector(
            self.plot_data[0], self.plot_data[1], self.plot_lons,
            self.lats.values, self.quiver_n_lon, self.quiver_n_lat,
            returnxy=True
        )

    def plot(self):
        if self.basemap:
            self.prep_data_for_basemap()
        else:
            self.plot_data = self.data[0]

        u, v, x, y = self.prep_quiver()
        self.handle = self.backend.quiver(x, y, u, v, **self.quiver_kwargs)

        if self.do_quiverkey:
            self.quiverkey = plt.quiverkey(self.handle, *self.quiverkey_args,
                                           **self.quiverkey_kwargs)
        if self.basemap:
            self.apply_basemap(self.basemap)
