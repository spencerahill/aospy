import scipy.stats
import numpy as np
from matplotlib import pyplot as plt
import mpl_toolkits.basemap
from mpl_toolkits.basemap import Basemap, shiftgrid

from .user import var_inst
from .calc import Calc

def _conv_to_dup_list(x, n, single_to_list=True):
    """
    Convert singleton or iterable into length-n list.  If the input is
    a list, with length-1, its lone value gets duplicated n times.  If
    the input is a list with length-n, leave it the same.  If the
    input is any other data type, replicate it as a length-n list.
    """
    if type(x) is list:
        if len(x) == n:
            return x
        elif len(x) == 1:
            return x*n
        else:
            raise ValueError("Input %s must have length 1 or %d : len(%s) = %d"
                             % (x, n, x, len(x)))
    else:
        if n == 1:
            if single_to_list:
                return [x]
            else:
                return x
        else:
            return [x]*n

class Fig(object):
    # Various categories of inputted specifications
    fig_specs = (
        'fig_title', 'n_row', 'n_col', 'row_size', 'col_size', 'n_ax',
        'subplot_lims', 'do_colorbar', 'cbar_ax_lim', 'cbar_ticks',
        'cbar_ticklabels', 'cbar_label', 'verbose'
    )
    ax_specs = (
        'n_plot','ax_title', 'do_ax_label', 'ax_left_label', 'ax_right_label',
        'map_proj', 'map_corners', 'map_res', 'shiftgrid_start',
        'shiftgrid_cyclic',
        'x_dim', 'x_lim', 'x_ticks', 'x_ticklabels', 'x_label',
        'y_dim', 'y_lim', 'y_ticks', 'y_ticklabels', 'y_label',
        'lat_lim', 'lat_ticks', 'lat_ticklabels', 'lat_label',
        'lon_lim', 'lon_ticks', 'lon_ticklabels', 'lon_label',
        'p_lim', 'p_ticks', 'p_ticklabels', 'p_label',
        'sigma_lim', 'sigma_ticks', 'sigma_ticklabels', 'sigma_label',
        'time_lim', 'time_ticks', 'time_ticklabels', 'time_label',
    )
    plot_specs = (
        'plot_type', 'marker_size', 'marker_shape', 'marker_color',
        'line_color', 'linestyle', 'do_best_fit_line', 'print_best_fit_slope',
        'print_corr_coeff', 'cntr_lvls', 'col_map', 'min_cntr', 'max_cntr',
        'num_cntr', 'contourf_extend', 'latlon_rect'
    )
    data_specs = (
        'proj', 'model', 'run', 'ens_mem', 'var', 'level', 'region', 'yr_range',
        'intvl_in', 'intvl_out', 'dtype_in_time', 'dtype_in_vert',
        'dtype_out_time', 'dtype_out_vert', 'do_subtract_mean'
    )
    specs = fig_specs + ax_specs + plot_specs + data_specs

    def __init__(self, n_ax=1, n_plot=1, n_data=1, n_row=1, n_col=1,
                 proj=None, model=None, run=None, ens_mem=None, var=None,
                 yr_range=None, region=None, intvl_in=None,
                 intvl_out=None, dtype_in_time=None, dtype_in_vert=None,
                 dtype_out_time=None, dtype_out_vert=None, level=None,
                 **kwargs):
        """Class for producing figures with one or more panels."""
        self.n_ax = n_ax
        self.n_plot = n_plot
        self.n_data = n_data
        self.n_row = n_row
        self.n_col = n_col
        self._set_n_ax_plot_data()

        self.proj = proj
        self.model = model
        self.run = run
        self.ens_mem = ens_mem
        self.var = var

        self.yr_range = yr_range
        self.region = region
        self.intvl_in = intvl_in
        self.intvl_out = intvl_out
        self.dtype_in_time = dtype_in_time
        self.dtype_in_vert = dtype_in_vert
        self.dtype_out_time = dtype_out_time
        self.dtype_out_vert = dtype_out_vert
        self.level = level

        self.Ax = []

        # Accept all other keyword arguments passed in as attrs.
        for key, val in kwargs.iteritems():
            setattr(self, key, val)
        # Hack solution to resolve mismatch between region specification and
        # dtype, i.e. if dtype doesn't include 'reg'.
        if 'reg' not in self.dtype_out_time and self.region:
            new_dtype = 'reg.' + self.dtype_out_time
            print ("Warning: region '%s' is non-null but dtype '%s' does not "
                   "specify region.  Setting dtype to '%s'."
                   % (self.region, self.dtype_out_time, new_dtype))
            self.dtype_out_time = new_dtype

        self.do_ax_label = True if self.n_ax > 1 else False

        self._expand_attrs_over_tree()
        self._make_ax_objs()
        # self._set_ax_geom()

    def _set_n_ax_plot_data(self):
        """Set the number of axes, plots, and data."""
        # Set the number of Axes.
        if self.n_ax == 'all':
            self.n_ax = self.n_row*self.n_col
        else:
            assert self.n_ax <= self.n_row*self.n_col
        # Distribute n_plot across the Axes.
        self.n_plot = _conv_to_dup_list(self.n_plot, self.n_ax)
        # Distribute n_data across the Plots.
        self.n_data = _conv_to_dup_list(self.n_data, self.n_ax)
        for i, ndt in enumerate(self.n_data):
            self.n_data[i] = _conv_to_dup_list(ndt, self.n_plot[i])

    def _traverse_child_tree(self, value, level='data'):
        """
        Traverse the "tree" of child Ax, Plot, and Data objects by creating a
        nested list of three levels, each level cascading one child object
        class down (e.g. from Ax to Plot to Data).
        """
        assert level in ('fig', 'ax', 'plot', 'data')
        if level in ('ax', 'plot', 'data'):
            value = _conv_to_dup_list(value, self.n_ax)
            if level in ('plot', 'data'):
                for i, vi in enumerate(value):
                    value[i] = _conv_to_dup_list(vi, self.n_plot[i])
                    if level == 'data':
                        for j, vij in enumerate(value[i]):
                            value[i][j] = _conv_to_dup_list(
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
        for specs, name in ((self.fig_specs, 'fig'),
                            (self.ax_specs, 'ax'),
                            (self.plot_specs, 'plot'),
                            (self.data_specs, 'data')):
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
        self.Ax = [Ax(self, n, self._locate_ax(n))
                   for n in range(self.n_ax)]

    def _make_colorbar(self):
        """Create colorbar for multi panel plots."""
        # Goes at bottom center if for all panels.
        self.cbar_ax = self.fig.add_axes(self.cbar_ax_lim)
        self.cbar = self.fig.colorbar(self.Ax[0].Plot[0].handle,
            cax=self.cbar_ax, orientation='horizontal', drawedges=False,
            spacing='proportional', extend=self.contourf_extend
        )
        # Set tick and label properties.
        if np.any(self.cbar_ticks):
            self.cbar.set_ticks(self.cbar_ticks)
        if self.cbar_ticklabels not in (None, False):
            self.cbar.set_ticklabels(self.cbar_ticklabels)
        self.cbar.ax.tick_params(labelsize='x-small')
        if self.cbar_label:
            var = var_inst(self.var[0][0][0])
            if self.cbar_label == 'units':
                if self.dtype_out_vert[0][0][0] == 'vert_int':
                    label = var.vert_int_plot_units
                else:
                    label = var.plot_units
            else:
                label = self.cbar_label
            self.cbar.set_label(label, fontsize='x-small',
                                labelpad=0.5)

    def create_fig(self):
        """Create the figure and set up the subplots."""
        self.fig = plt.figure(figsize=(self.n_col*self.col_size,
                                       self.n_row*self.row_size))
        self.fig.subplots_adjust(**self.subplot_lims)
        if self.fig_title:
            self.fig.suptitle(self.fig_title, fontsize=12)

        for n in range(self.n_ax):
            self.Ax[n].ax = self.fig.add_subplot(self.n_row, self.n_col, n+1)
            self.Ax[n]._set_axes_props()
            self.Ax[n]._set_axes_labels()

    def make_plots(self):
        """Render the plots in every Ax."""
        for n in range(self.n_ax):
            self.Ax[n].make_plots()
        if self.do_colorbar:
            self._make_colorbar()

    def savefig(self, *args, **kwargs):
        """Save the Fig using matplotlib's built-in 'savefig' method."""
        self.fig.savefig(*args, **kwargs)

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
        self.Fig = Fig
        self.ax_specs = self.Fig.ax_specs
        self.plot_specs = self.Fig.plot_specs
        self.data_specs = self.Fig.data_specs
        self.ax_num = ax_num
        self.ax_loc = ax_loc
        self.n_plot = Fig.n_plot[ax_num]
        self.n_data = Fig.n_data[ax_num]
        self.Plot = []

        self._copy_attrs_from_fig()
        self._set_ax_loc_specs()
        self._set_xy_attrs_to_coords()
        self._make_plot_objs()

    def _traverse_child_tree(self, value, level='data'):
        """Traverse the "tree" of child Plot, and Data objects."""
        assert level in ('ax', 'plot', 'data')
        if level in ('plot', 'data'):
            value = _conv_to_dup_list(value, self.n_plot)
            if level in ('data'):
                for i, vi in enumerate(value):
                    value[i] = _conv_to_dup_list(vi, self.n_data[i])
        return value

    def _copy_attrs_from_fig(self):
        """Copy the attrs of the parent Fig that correspond to this Ax."""
        for attr in self.ax_specs:
            value = getattr(self.Fig, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'ax'))

        for attr in self.plot_specs:
            value = getattr(self.Fig, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'plot'))

        for attr in self.data_specs:
            value = getattr(self.Fig, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'data'))

    def _set_ax_loc_specs(self):
        """Set attrs that depend on Ax location within the Fig."""
        # Take the Fig's attr value if it's neeeded; otherwise set False.
        for key, val in self.labels[self.ax_loc].iteritems():
            if val:
                if val == ' ':
                    new_val = ' '
                else:
                    new_val = getattr(self.Fig, key, False)

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
        if self.y_lim:
            self.ax.set_ylim(self.y_lim)
        if self.y_ticks:
            self.ax.set_yticks(self.y_ticks)
        if self.y_ticklabels:
            self.ax.set_yticklabels(self.y_ticklabels, fontsize='small')
        if self.y_label:
            self.ax.set_ylabel(self.y_label, fontsize='small', labelpad=-2)
        if self.x_lim:
            if self.x_lim == 'ann_cycle':
                self.x_lim = (1,12)
                self.x_ticks = range(1,13)
                self.x_ticklabels = tuple('JFMAMJJASOND')
            self.ax.set_xlim(self.x_lim)
        if self.x_ticks:
            self.ax.set_xticks(self.x_ticks)
        if self.x_ticklabels:
            self.ax.set_xticklabels(self.x_ticklabels, fontsize='x-small')
        if self.x_label:
            self.ax.set_xlabel(self.x_label, fontsize='small', labelpad=1)

        plt.tick_params(labelsize='x-small')
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
        if self.do_ax_label:
            self.panel_label = self.ax.text(
                0.04, 0.9, '(%s)' % tuple('abcdefghijklmnopqrs')[self.ax_num],
                fontsize='small', transform=self.ax.transAxes
                )
        # Labels to left and/or right of Axis.
        if self.ax_left_label:
            # if self.ax_left_label_rot == 'horizontal':
                # horiz_frac = -0.17
            # else:
            horiz_frac = -0.05
            vert_frac = 0.5
            self.ax.text(
                horiz_frac, vert_frac, self.ax_left_label,
                verticalalignment='center', horizontalalignment='left',
                rotation='vertical', fontsize='small',
                transform=self.ax.transAxes
            )
        if self.ax_right_label:
            horiz_frac = 1.02; vert_frac = 0.5
            self.ax.text(
                horiz_frac, vert_frac, self.ax_right_label,
                verticalalignment='center', rotation='vertical',
                fontsize='small', transform=self.ax.transAxes
            )

    def _make_plot_objs(self):
        """Create the Plot object for each plotted element."""
        plots = {'scatter': Scatter, 'map': Map, 'contour': Contour,
                 'contourf': Map, 'line': Line}
        self.Plot = []
        for n in range(self.n_plot):
            try:
                p = plots[self.plot_type[n]]
            except KeyError:
                raise TypeError("Plot type '%s' not recognized."
                                % self.plot_type[n])
            else:
                self.Plot.append(p(self, n))

    def make_plots(self):
        """Call the matplotlib plotting command for each Plot."""
        for n in range(self.n_plot):
            self.Plot[n].plot()

class Plot(object):
    def __init__(self, Ax, plot_num):
        """Class for plotting single data element."""
        self.Ax = Ax
        self.Fig = Ax.Fig
        self.plot_specs = self.Ax.plot_specs
        self.data_specs = self.Ax.data_specs
        self.plot_num = plot_num
        self.n_data = Ax.n_data[plot_num]
        self._copy_attrs_from_ax()
        if type(self.var[0]) is tuple:
            self.var = self.var[0]
            self.n_data = len(self.var)
        self.Calc = self._make_calc_obj()
        self._set_coord_arrays()
        self._set_cntr_lvls()
        self.data = self._load_data()
        self._apply_data_transforms()

    def _traverse_child_tree(self, value, level='data'):
        """Traverse the "tree" of child Data objects."""
        assert level in ('plot', 'data')
        if level == 'data':
            value = _conv_to_dup_list(value, self.n_data)
        return value

    def _copy_attrs_from_ax(self):
        """Copy the attrs of the parent Ax that correspond to this Plot."""
        for attr in self.plot_specs:
            value = getattr(self.Ax, attr)[self.plot_num]
            setattr(self, attr, self._traverse_child_tree(value, 'plot'))

        for attr in self.data_specs:
            value = getattr(self.Ax, attr)[self.plot_num]
            setattr(self, attr, self._traverse_child_tree(value, 'data'))

    def _make_calc_obj(self):
        calc_obj = []
        # for i in range(max(self.n_data, len(self.var[0]))):
        for i in range(self.n_data):
            if type(self.run[i]) is dict:
                calc_pair = [Calc(
                    proj=self.proj[i], model=self.model[i],
                    run=self.run[i].keys()[0][j],
                    ens_mem=self.ens_mem[i], var=self.var[i],
                    yr_range=self.yr_range[i], region=self.region[i],
                    intvl_in=self.intvl_in[i], intvl_out=self.intvl_out[i],
                    dtype_in_time=self.dtype_in_time[i],
                    dtype_in_vert=self.dtype_in_vert[i],
                    dtype_out_time=self.dtype_out_time[i],
                    dtype_out_vert=self.dtype_out_vert[i], level=self.level[i],
                    verbose=self.Fig.verbose, skip_time_inds=True
                ) for j in range(2)]
                calc_obj.append(calc_pair)
            else:
                calc_obj.append(Calc(
                    proj=self.proj[i], model=self.model[i], run=self.run[i],
                    ens_mem=self.ens_mem[i], var=self.var[i],
                    yr_range=self.yr_range[i], region=self.region[i],
                    intvl_in=self.intvl_in[i], intvl_out=self.intvl_out[i],
                    dtype_in_time=self.dtype_in_time[i],
                    dtype_in_vert=self.dtype_in_vert[i],
                    dtype_out_time=self.dtype_out_time[i],
                    dtype_out_vert=self.dtype_out_vert[i], level=self.level[i],
                    verbose=self.Fig.verbose, skip_time_inds=True
                ))
        return calc_obj

    def _set_plot_func(self):
        """
        Set the matplotlib plotting function that will render this plot.

        Can be the name of any pyplot plotting function, e.g. 'contour',
        'contourf', 'plot', 'scatter', 'imshow', etc.
        """
        self.plot_func = getattr(self.Ax.ax, self.plot_type)

    def _set_coord_arrays(self):
        """Set the arrays holding the x- and y-coordinates."""
        array_names = {'lat': 'lat', 'lon': 'lon', 'p': 'level',
                       'sigma': 'pfull', 'time': 'time'}
        if self.n_data == 1:
            mod_inds = (0, 0)
        else:
            mod_inds = (0, 1)
        for dim, data, lim, i in zip(('x_dim', 'y_dim'), ('x_data', 'y_data'),
                                     ('x_lim', 'y_lim'), mod_inds):
            array_key = getattr(self.Ax, dim)
            try:
                array = getattr(self.Calc[i].model[0], array_names[array_key])
            except AttributeError:
                array = getattr(self.Calc[i][0].model[0], array_names[array_key])
            # Restrict the array to the specified limits.
            lim = getattr(self.Ax, lim)
            # Annual cycle has 12 evenly spaced points, one per month.
            if array_key == 'time' and lim == 'ann_cycle':
                array = np.arange(1,13)
            else:
                # For sigma coordinates, divide the pfull array by a reference
                # surface pressure of 1000 hPa to get a fractional value.
                if array_key == 'sigma':
                    array /= 1e3
                # if lim:
                    # array = array[np.where((array <= max(lim)) &
                                           # (array >= min(lim)))]
            setattr(self, data, array)

    def _set_cntr_lvls(self):
        self.cntr_lvls = np.linspace(
            self.min_cntr, self.max_cntr, self.num_cntr + 1
        )

    def _load_data(self):
        try:
            return tuple(
                [cl.load(dto, dtype_out_vert=dtv, region=reg, time=False,
                         vert=lev, lat=False, lon=False, plot_units=True)
                 for cl, dto, dtv, reg, lev in zip(self.Calc,
                                                   self.dtype_out_time,
                                                   self.dtype_out_vert,
                                                   self.region, self.level)]
            )
        except AttributeError:
            data = tuple(
                [cl.load(self.dtype_out_time[0],
                         dtype_out_vert=self.dtype_out_vert[0],
                         region=self.region[0], time=False, vert=self.level[0],
                         lat=False, lon=False, plot_units=True)
                 for cl in self.Calc[0]]
            )
            self.Calc = self.Calc[0]
            return data[0] - data[1]
    #     if self.n_data == 1:
    #         self._load_2d_data()
    #     elif self.n_data == 2:
    #         self._load_xy_data()

    # def _load_2d_data(self):
    #     """Load 2D-data and set the corresponding attrs."""
    #     self.data = np.squeeze(load_plot_data(
    #         self.proj, self.model, self.run, self.ens_mem[0],
    #         self.var, self.level[0], self.intvl[0], self.dtype[0],
    #         self.yr_range[0], region=self.region[0]
    #     ))

    # def _load_xy_data(self):
    #     """Load x- and y-data and set the corresponding attrs."""
    #     for i, data in enumerate(('x_data', 'y_data')):
    #         setattr(self, data, load_plot_data(
    #             self.proj[i], self.model[i], self.run[i], self.ens_mem[i],
    #             self.var[i], self.level[i], self.intvl[i], self.dtype[i],
    #             self.yr_range[i], region=self.region[i]
    #         ))

    def _subtract_mean(self, data):
        return np.subtract(data, np.mean(data))

    # def _mult_by_factor(self, data):
        # return np.multiply(data, self.mult_factor)

    def _apply_data_transforms(self):
        """Apply any specified transformations to the data once loaded."""
        transforms = {'do_subtract_mean': self._subtract_mean}
                      #'mult_by_factor': self._mult_by_factor}
        for attr, method in transforms.iteritems():
            for data, do_method in zip(['x_data', 'y_data'], getattr(self, attr)):
                if do_method:
                    setattr(self, data, method(getattr(self, data)))
        return data

class Line(Plot):
    def __init__(self, Ax, plot_num):
        """Line plot."""
        Plot.__init__(self, Ax, plot_num)

    def plot(self):
        handle = self.Ax.ax.plot(
            self.x_data, self.y_data, color=self.line_color,
            linestyle=self.linestyle, marker=self.marker_shape,
            markerfacecolor=self.marker_color, markersize=self.marker_size
        )
        self.handle = handle

class Scatter(Plot):
    def __init__(self, Ax, plot_num):
        Plot.__init__(self, Ax, plot_num)
        assert self.n_data == 2

    def plot(self):
        self.Ax.ax.scatter(
            self.x_data, self.y_data, s=self.marker_size, c=self.marker_color,
            marker=self.marker_shape
        )
        if self.do_best_fit_line:
            self.best_fit_line(print_slope=self.print_best_fit_slope)

        if self.print_corr_coeff:
            self.corr_coeff()

    def corr_coeff(self):
        """Compute the Pearson correlation coefficient and plot it."""
        pearsonr, p_val = scipy.stats.pearsonr(self.x_data, self.y_data)
        self.Ax.ax.text(0.3, 0.89, r'$r=$ %.2f' % pearsonr,
                        transform=self.Ax.ax.transAxes, fontsize='x-small')

    def best_fit_line(self, print_slope=True):
        """Plot the best fit line to the scatterplot data."""
        best_fit = np.polyfit(self.x_data, self.y_data, 1)
        x_lin_fit = [-1e3,1e3]
        def lin_fit(m, x, b):
            return [m*xx + b for xx in x]
        self.Ax.ax.plot(x_lin_fit, lin_fit(best_fit[0], x_lin_fit,
                                           best_fit[1]), 'k')
        if print_slope:
            self.Ax.ax.text(0.3, 0.07, r'slope = %0.2f' % best_fit[0],
                            transform=self.Ax.ax.transAxes, fontsize='x-small')

class Map(Plot):
    def __init__(self, Ax, plot_num):
        Plot.__init__(self, Ax, plot_num)
        # Assume single data source for mapped data.
        self.model = self.Calc[0].model[0]
        if type(self.data) in (list, tuple):
            self.data = self.data[0]
        self.lon = self.model.lon
        self.lat = self.model.lat

    def _make_basemap(self):
        if self.Ax.x_lim:
            llcrnrlon, urcrnrlon = self.Ax.x_lim
        else:
            llcrnrlon, urcrnrlon = -180, 180
        if self.Ax.y_lim:
            llcrnrlat, urcrnrlat = self.Ax.y_lim
        else:
            llcrnrlat, urcrnrlat = -90, 90

        self.basemap = Basemap(
            projection=self.Ax.map_proj, resolution=self.Ax.map_res,
            ax=self.Ax.ax, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat
        )

    def _shiftgrid(self, lon0, datain, lonsin, start=True, cyclic=360.0):
        return shiftgrid(lon0, datain, lonsin, start=start, cyclic=cyclic)

    def plot_latlon_rect(self, lonmin, lonmax, latmin, latmax):
        """Plot a lon-lat rectangle in Basemap.  Taken from:
        http://stackoverflow.com/questions/23941738/python-draw-rectangle-in-basemap
        """
        xs = [lonmin, lonmax, lonmax, lonmin, lonmin]
        ys = [latmin, latmin, latmax, latmax, latmin]
        return self.basemap.plot(xs, ys, latlon=True)

    def plot(self):
        self._make_basemap()
        if self.Ax.shiftgrid_start:
            lon0 = 180
        else:
            lon0 = 180
        plot_data, lons = self._shiftgrid(
            lon0, self.data, self.lon,
            start=self.Ax.shiftgrid_start, cyclic=self.Ax.shiftgrid_cyclic
        )
        x, y = self.basemap(*np.meshgrid(lons, self.lat))
        plot_data = mpl_toolkits.basemap.maskoceans(
            x, y, plot_data, inlands=False, resolution='c'
        )
        self.handle = self.basemap.contourf(
            x, y, plot_data, self.cntr_lvls, extend=self.contourf_extend
        )
        self.handle.set_cmap(self.col_map)

        self.basemap.drawcoastlines(linewidth=0.3)
        self.basemap.drawmapboundary(linewidth=0.3, fill_color='0.85')
        if self.latlon_rect:
            self.plot_latlon_rect(*self.latlon_rect)

class Contour(Plot):
    def __init__(self, Ax, plot_num):
        """Contour or contourf plot."""
        Plot.__init__(self, Ax, plot_num)

    def plot(self):
        # self._set_plot_func()
        self.handle = self.contourf(
            self.x_data, self.y_data, self.data.T, self.cntr_lvls,
            extend=self.contourf_extend
        )
        self.handle.set_cmap(self.col_map)


    # def _set_ax_geom(self):
    #     """Set row and column sizes and subplot spacing."""
    #     if self.n_col == 1:
    #         if self.n_row == 1:
    #             self.subplot_lims = {
    #                 'left': 0.16, 'right': 0.95, 'wspace': 0.05,
    #                 'bottom': 0.2, 'top': 0.85, 'hspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.2, 0.15, 0.7, 0.04]
    #         elif self.n_row == 2:
    #             self.subplot_lims = {
    #                 'left': 0.1, 'right': 0.93, 'hspace': 0.2,
    #                 'bottom': 0.18, 'top': 0.9, 'wspace': 0.05
    #             }
    #             self.cbar_ax_lim = [0.2, 0.1, 0.7, 0.025]
    #         elif self.n_row == 3:
    #             self.subplot_lims = {
    #                 'left': 0.16, 'right': 0.93, 'hspace': 0.15,
    #                 'bottom': 0.16, 'top': 0.93, 'wspace': 0.05
    #             }
    #             self.cbar_ax_lim = [0.18, 0.09, 0.7, 0.015]
    #         elif self.n_row ==4:
    #             self.subplot_lims = {
    #                 'left': 0.1, 'right': 0.93, 'hspace': 0.15,
    #                 'bottom': 0.1, 'top': 0.93, 'wspace': 0.05
    #             }
    #             self.cbar_ax_lim = [0.2, 0.07, 0.7, 0.02]
    #         else:
    #             self.subplot_lims = {
    #                 'left': 0.08, 'right': 0.95, 'hspace': 0.2,
    #                 'bottom': 0.1, 'top': 0.88, 'wspace': 0.05
    #             }
    #             self.cbar_ax_lim = [0.2, 0.05, 0.6, 0.025]
    #     elif self.n_col == 2:
    #         if self.n_row == 1:
    #             self.subplot_lims = {
    #                 'left': 0.1, 'right': 0.95, 'hspace': 0.15,
    #                 'bottom': 0.25, 'top': 0.8, 'wspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.325, 0.11, 0.4, 0.04]
    #         elif self.n_row == 2:
    #             self.subplot_lims = {
    #                 'left': 0.1, 'right': 0.95, 'hspace': 0.2,
    #                 'bottom': 0.15, 'top': 0.88, 'wspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.325, 0.08, 0.4, 0.02]
    #         elif self.n_row == 3:
    #             self.subplot_lims = {
    #                 'left': 0.1, 'right': 0.95, 'hspace': 0.15,
    #                 'bottom': 0.11, 'top': 0.92, 'wspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.325, 0.1, 0.4, 0.02]
    #         elif self.n_row == 4:
    #             self.subplot_lims = {
    #                 'left': 0.1, 'right': 0.95, 'hspace': 0.1,
    #                 'bottom': 0.1, 'top': 0.94, 'wspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.27, 0.1, 0.5, 0.013]
    #         else:
    #             self.subplot_lims = {
    #                 'left': 0.1, 'right': 0.95, 'hspace': 0.2,
    #                 'bottom': 0.1, 'top': 0.88, 'wspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.37, 0.07, 0.3, 0.03]
    #     elif self.n_col == 3:
    #         if self.n_row == 1:
    #             self.subplot_lims = {
    #                 'left': 0.06, 'right': 0.95, 'hspace': 0.2,
    #                 'bottom': 0.25, 'top': 0.85, 'wspace': 0.05
    #             }
    #             self.cbar_ax_lim = [0.35, 0.14, 0.3, 0.04]
    #         elif self.n_row == 2:
    #             self.subplot_lims = {
    #                 'left': 0.08, 'right': 0.95, 'hspace': 0.2,
    #                 'bottom': 0.15, 'top': 0.88, 'wspace': 0.05
    #             }
    #             self.cbar_ax_lim = [0.37, 0.07, 0.3, 0.03]
    #         elif self.n_row == 3:
    #             self.subplot_lims = {
    #                 'left': 0.06, 'right': 0.95, 'hspace': 0.15,
    #                 'bottom': 0.09, 'top': 0.92, 'wspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.37, 0.07, 0.3, 0.03]
    #         elif self.n_row == 4:
    #             self.subplot_lims = {
    #                 'left': 0.06, 'right': 0.95, 'hspace': 0.15,
    #                 'bottom': 0.09, 'top': 0.94, 'wspace': 0.1
    #             }
    #             self.cbar_ax_lim = [0.355, 0.06, 0.3, 0.013]
    #         else:
    #             self.subplot_lims = {
    #                 'left': 0.08, 'right': 0.95, 'hspace': 0.2,
    #                 'bottom': 0.1, 'top': 0.88, 'wspace': 0.05
    #             }
    #             self.cbar_ax_lim = [0.37, 0.07, 0.3, 0.03]
    #     elif self.n_col == 4:
    #         if self.n_row == 1:
    #             self.subplot_lims = {
    #                 'left': 0.07, 'right': 0.98, 'hspace': 0.08,
    #                 'bottom': 0.15, 'top': 0.8, 'wspace': 0.1
    #             }
    #         else:
    #             self.subplot_lims = {
    #                 'left': 0.06, 'right': 0.95, 'hspace': 0.2,
    #                 'bottom': 0.25, 'top': 0.85, 'wspace': 0.05
    #             }
    #     else:
    #         self.subplot_lims = {
    #             'left': 0.08, 'right': 0.95, 'hspace': 0.2,
    #             'bottom': 0.1, 'top': 0.88, 'wspace': 0.05
    #         }
    #         self.cbar_ax_lim = [0.2, 0.05, 0.6, 0.02]
