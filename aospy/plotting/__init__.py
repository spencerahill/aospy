"""Plotting submodule of aospy module."""

def _conv_to_dup_list(x, n):
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
            return [x[0]]*n
        else:
            raise ValueError("Input %s must have length 1 or %d : len(%s) = %d"
                             % (x, n, x, len(x)))
    else:
        return [x]*n

class Figure(object):
    # Various categories of inputted specifications
    fig_specs = (
        'fig_title', 'n_row', 'n_col', 'row_size', 'col_size', 'n_ax'
    )
    ax_specs = (
        'n_plot','ax_title', 'do_ax_label', 'xlim', 'xticks', 'ylim', 'yticks',
        'map_proj', 'map_corners', 'map_res', 'shiftgrid_start',
        'shiftgrid_cyclic'
    )
    plot_specs = (
        'plot_type', 'marker_size', 'marker_shape', 'marker_color',
        'do_best_fit_line', 'print_best_fit_slope', 'print_corr_coeff',
        'cntr_lvls', 'col_map', 'min_cntr', 'max_cntr', 'num_cntr',
        'contourf_extend'
    )
    data_specs = (
        'proj', 'model', 'run', 'ens_mem', 'var', 'region', 'level', 'intvl',
        'dtype', 'yr_range', 'do_subtract_mean'
    )
    specs = fig_specs + ax_specs + plot_specs + data_specs

    def __init__(self, n_ax=1, n_plot=1, n_data=1, n_row=1, n_col=1,
                 proj=None, model=None, run=None, var=None, **kwargs):
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
        self.var = var
        self.Ax = []
        # Accept all other keyword arguments passed in as attrs.
        for key, val in kwargs.iteritems():
            setattr(self, key, val)

        self.do_ax_label = True if self.n_ax > 1 else False

        self._expand_attrs_over_tree()
        self._make_ax_objs()
        self._set_ax_geom()

    def _set_n_ax_plot_data(self):
        """Number of axes must be <= rows*columns."""
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
                                vij, self.n_data[i][j]
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

    def _set_ax_geom(self):
        """Set row and column sizes and subplot spacing."""
        if self.n_col == 1:
            if self.n_row == 1:
                self.panels_adjust = {
                    'left': 0.16, 'right': 0.95, 'wspace': 0.05,
                    'bottom': 0.2, 'top': 0.85, 'hspace': 0.1
                }
                # self.ax_cbar = [0.2, 0.15, 0.7, 0.04]
            elif self.n_row == 2:
                self.panels_adjust = {
                    'left': 0.16, 'right': 0.93, 'hspace': 0.2,
                    'bottom': 0.18, 'top': 0.9, 'wspace': 0.05
                }
                # self.ax_cbar = [0.2, 0.1, 0.7, 0.025]
                self.panels_adjust = {
                    'left': 0.05, 'right': 0.95, 'hspace': 0.15,
                    'bottom': 0.18, 'top': 0.9, 'wspace': 0.05
                }
                # self.ax_cbar = [0.17, 0.13, 0.66, 0.03]
            elif self.n_row == 3:
                self.panels_adjust = {
                    'left': 0.16, 'right': 0.93, 'hspace': 0.15,
                    'bottom': 0.16, 'top': 0.93, 'wspace': 0.05
                }
                # self.ax_cbar = [0.18, 0.09, 0.7, 0.015]
            elif self.n_row ==4:
                self.panels_adjust = {
                    'left': 0.1, 'right': 0.93, 'hspace': 0.15,
                    'bottom': 0.1, 'top': 0.93, 'wspace': 0.05
                }
                # self.ax_cbar = [0.2, 0.07, 0.7, 0.02]
            else:
                self.panels_adjust = {
                    'left': 0.08, 'right': 0.95, 'hspace': 0.2,
                    'bottom': 0.1, 'top': 0.88, 'wspace': 0.05
                }
                # self.ax_cbar = [0.37, 0.07, 0.3, 0.03]
        elif self.n_col == 2:
            if self.n_row == 1:
                self.panels_adjust = {
                    'left': 0.1, 'right': 0.95, 'hspace': 0.15,
                    'bottom': 0.25, 'top': 0.8, 'wspace': 0.1
                }
                # self.ax_cbar = [0.325, 0.11, 0.4, 0.04]
            elif self.n_row == 2:
                self.panels_adjust = {
                    'left': 0.1, 'right': 0.95, 'hspace': 0.2,
                    'bottom': 0.15, 'top': 0.88, 'wspace': 0.1
                }
                # self.ax_cbar = [0.325, 0.08, 0.4, 0.02]
            elif self.n_row == 3:
                self.panels_adjust = {
                    'left': 0.1, 'right': 0.95, 'hspace': 0.15,
                    'bottom': 0.11, 'top': 0.92, 'wspace': 0.1
                }
                # self.ax_cbar = [0.325, 0.052, 0.4, 0.02]
            elif self.n_row == 4:
                self.panels_adjust = {
                    'left': 0.1, 'right': 0.95, 'hspace': 0.1,
                    'bottom': 0.1, 'top': 0.94, 'wspace': 0.1
                }
                # self.ax_cbar = [0.27, 0.05, 0.5, 0.013]
            else:
                self.panels_adjust = {
                    'left': 0.1, 'right': 0.95, 'hspace': 0.2,
                    'bottom': 0.1, 'top': 0.88, 'wspace': 0.1
                }
                # self.ax_cbar = [0.37, 0.07, 0.3, 0.03]
        elif self.n_col == 3:
            if self.n_row == 1:
                self.panels_adjust = {
                    'left': 0.06, 'right': 0.95, 'hspace': 0.2,
                    'bottom': 0.25, 'top': 0.85, 'wspace': 0.05
                }
                # self.ax_cbar = [0.35, 0.14, 0.3, 0.04]
            elif self.n_row == 2:
                self.panels_adjust = {
                    'left': 0.08, 'right': 0.95, 'hspace': 0.2,
                    'bottom': 0.15, 'top': 0.88, 'wspace': 0.05
                }
                # self.ax_cbar = [0.37, 0.07, 0.3, 0.03]
            elif self.n_row == 3:
                self.panels_adjust = {
                    'left': 0.06, 'right': 0.95, 'hspace': 0.15,
                    'bottom': 0.09, 'top': 0.92, 'wspace': 0.1
                }
                # self.ax_cbar = [0.37, 0.07, 0.3, 0.03]
            elif self.n_row == 4:
                self.panels_adjust = {
                    'left': 0.06, 'right': 0.95, 'hspace': 0.15,
                    'bottom': 0.09, 'top': 0.94, 'wspace': 0.1
                }
                # self.ax_cbar = [0.355, 0.06, 0.3, 0.013]
            else:
                self.panels_adjust = {
                    'left': 0.08, 'right': 0.95, 'hspace': 0.2,
                    'bottom': 0.1, 'top': 0.88, 'wspace': 0.05
                }
                # self.ax_cbar = [0.37, 0.07, 0.3, 0.03]
        elif self.n_col == 4:
            if self.n_row == 1:
                self.panels_adjust = {
                    'left': 0.07, 'right': 0.98, 'hspace': 0.08,
                    'bottom': 0.15, 'top': 0.8, 'wspace': 0.1
                }
            else:
                self.panels_adjust = {
                    'left': 0.06, 'right': 0.95, 'hspace': 0.2,
                    'bottom': 0.25, 'top': 0.85, 'wspace': 0.05
                }
        else:
            self.panels_adjust = {
                'left': 0.08, 'right': 0.95, 'hspace': 0.2,
                'bottom': 0.1, 'top': 0.88, 'wspace': 0.05
            }
            # self.ax_cbar = [0.37, 0.07, 0.3, 0.03]

    def _make_colorbar(self, contours, label, for_all_panels=True, 
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
            cbar = plt.gcf().colorbar(
                contours, ax=plt.gca(), drawedges=False,
                orientation='vertical', spacing='proportional', extend=extend,
                fraction=0.1, aspect=16, pad=0.03
            )
        # Set tick and label properties.
        if np.any(ticks):
            cbar.set_ticks(ticks)
        if ticklabels not in (None, False):
            cbar.set_ticklabels(ticklabels)
        cbar.ax.tick_params(labelsize='x-small')
        cbar.set_label(label, fontsize='x-small', labelpad=0.5)

    def create_fig(self):
        """Create the figure and set up the subplots."""
        from matplotlib import pyplot as plt

        self.fig = plt.figure(figsize=(self.n_col*self.col_size,
                                       self.n_row*self.row_size))
        self.fig.subplots_adjust(**self.panels_adjust)
        if self.fig_title:
            self.fig.suptitle(self.fig_title, fontsize=11)

        for n in range(self.n_ax):
            self.Ax[n].ax = self.fig.add_subplot(self.n_row, self.n_col, n+1)
            self.Ax[n]._set_axes_props()

    def make_plots(self):
        """Render the plots in every Ax."""
        for n in range(self.n_ax):
            self.Ax[n].make_plots()
        if self.do_colorbar:
            pass
            # self._make_colorbar()

    def savefig(self, *args, **kwargs):
        """Save the Figure using matplotlib's built-in 'savefig' method."""
        self.fig.savefig(*args, **kwargs)

class Ax(object):
    # Which labels to include based on position in the figure.
    labels = {
        'left': {'xticklabels': ' ', 'yticklabels': True,
                 'xlabel': False, 'ylabel': True, 'do_colorbar': False},
        'interior': {'xticklabels': ' ', 'yticklabels': ' ',
                     'xlabel': False, 'ylabel': False, 'do_colorbar': False},
        'bottomleft': {'xticklabels': True, 'yticklabels': True,
                       'xlabel': True, 'ylabel': True, 'do_colorbar': True},
        'bottom': {'xticklabels': True, 'yticklabels': ' ',
                   'xlabel': True, 'ylabel': False, 'do_colorbar': True}
    }

    def __init__(self, Figure, ax_num, ax_loc):
        self.Figure = Figure
        self.ax_specs = self.Figure.ax_specs
        self.plot_specs = self.Figure.plot_specs
        self.data_specs = self.Figure.data_specs
        self.ax_num = ax_num
        self.ax_loc = ax_loc
        self.n_plot = Figure.n_plot[ax_num]
        self.n_data = Figure.n_data[ax_num]
        self.Plot = []

        self._copy_attrs_from_fig()
        self._set_ax_loc_specs()
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
        """Copy the attrs of the parent Figure that correspond to this Ax."""
        for attr in self.ax_specs:
            value = getattr(self.Figure, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'ax'))

        for attr in self.plot_specs:
            value = getattr(self.Figure, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'plot'))

        for attr in self.data_specs:
            value = getattr(self.Figure, attr)[self.ax_num]
            setattr(self, attr, self._traverse_child_tree(value, 'data'))

    def _set_ax_loc_specs(self):
        """Set attrs that depend on Ax location within the Figure."""
        # Take the Figure's attr value if it's neeeded; otherwise set False.
        for key, val in self.labels[self.ax_loc].iteritems():
            if val:
                if val == ' ':
                    new_val = ' '
                else:
                    new_val = getattr(self.Figure, key, False)

            else:
                new_val = False
            setattr(self, key, new_val)

    def _set_axes_props(self):
        """Set the properties of the matplotlib Axes instance."""
        from matplotlib import pyplot as plt
        if self.ylim:
            self.ax.set_ylim(self.ylim)
        if self.yticks:
            self.ax.set_yticks(self.yticks)
        if self.yticklabels:
            self.ax.set_yticklabels(self.yticklabels, fontsize='small')
        if self.ylabel:
            # if self.ylabel == 'units': ylabel = var.plot_units
            self.ax.set_ylabel(self.ylabel, fontsize='small', labelpad=-2)
        if self.xlim:
            if self.xlim == 'ann_cycle':
                self.xlim = (1,12)
                self.xticks = range(1,13)
                self.xticklabels = tuple('JFMAMJJASOND')
            self.ax.set_xlim(self.xlim)
        if self.xticks:
            self.ax.set_xticks(self.xticks)
        if self.xticklabels:
            self.ax.set_xticklabels(self.xticklabels, fontsize='x-small')
        if self.xlabel:
            # if self.xlabel == 'units': self.xlabel = var.plot_units
            self.ax.set_xlabel(self.xlabel, fontsize='small', labelpad=1)

        plt.tick_params(labelsize='x-small')
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.yaxis.set_ticks_position('left')

        # Axis title.
        if self.ax_title:
            self.ax.set_title(self.ax_title, fontsize='small')
        # Axis panel labels, i.e. (a), (b), (c), etc.
        if self.do_ax_label:
            self.panel_label = self.ax.text(
                0.04, 0.9, '(%s)' % tuple('abcdefghijklmnopqrs')[self.ax_num],
                fontsize='small', transform=self.ax.transAxes
                )

    def _make_plot_objs(self):
        """Create the Plot object for each plotted element."""
        plots = {'scatter': Scatter, 'map': Map}
        self.Plot = []
        for n in range(self.n_plot):
            try:
                p = plots[self.plot_type[n]](self, n)
            except KeyError:
                raise TypeError("Plot type '%s' not recognized."
                                % self.plot_type[n])
            else:
                self.Plot.append(p)

    def make_plots(self):
        """Call the matplotlib plotting command for each Plot."""
        for n in range(self.n_plot):
            self.Plot[n].plot()

class Plot(object):
    def __init__(self, Ax, plot_num):
        self.Ax = Ax
        self.Figure = Ax.Figure
        self.plot_specs = self.Ax.plot_specs
        self.data_specs = self.Ax.data_specs
        self.plot_num = plot_num
        self.n_data = Ax.n_data[plot_num]

        self._copy_attrs_from_ax()
        self._conv_to_aospy_obj()
        self._load_data()
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

    def _conv_to_aospy_obj(self):
        """Convert any string labels of aospy objects to those objects."""
        from aospy.io import _aospy_inst
        self.proj, self.model, self.run, self.var = _aospy_inst(
            proj=self.proj, model=self.model, run=self.run, var=self.var)

    def _load_data(self):
        if self.n_data == 1:
            self._load_2d_data()
        elif self.n_data == 2:
            self._load_xy_data()

    def _load_2d_data(self):
        """Load 2D-data and set the corresponding attrs."""
        import numpy as np
        from aospy.io import load_plot_data
        self.data = np.squeeze(load_plot_data(
            self.proj, self.model, self.run, self.ens_mem[0],
            self.var, self.level[0], self.intvl[0], self.dtype[0],
            self.yr_range[0], region=self.region[0]
        ))

    def _load_xy_data(self):
        """Load x- and y-data and set the corresponding attrs."""
        from aospy.io import load_plot_data
        for i, data in enumerate(('xdata', 'ydata')):
            setattr(self, data, load_plot_data(
                self.proj[i], self.model[i], self.run[i], self.ens_mem[i],
                self.var[i], self.level[i], self.intvl[i], self.dtype[i],
                self.yr_range[i], region=self.region[i]
            ))

    def _subtract_mean(self, data):
        import numpy as np
        return np.subtract(data, np.mean(data))

    def _apply_data_transforms(self):
        """Apply any specified transformations to the data once loaded."""
        import numpy as np
        transforms = {'do_subtract_mean': self._subtract_mean}
        for attr, method in transforms.iteritems():
            for data, do_method in zip(['xdata', 'ydata'], getattr(self, attr)):
                if do_method:
                    setattr(self, data, method(getattr(self, data)))
        return data

class Scatter(Plot):
    def __init__(self, Ax, plot_num):
        Plot.__init__(self, Ax, plot_num)
        assert self.n_data == 2

    def plot(self):
        self.Ax.ax.scatter(
            self.xdata, self.ydata, s=self.marker_size, c=self.marker_color,
            marker=self.marker_shape
        )
        if self.do_best_fit_line:
            self.best_fit_line(print_slope=self.print_best_fit_slope)

        if self.print_corr_coeff:
            self.corr_coeff()

    def corr_coeff(self):
        """Compute the Pearson correlation coefficient and plot it."""
        import scipy.stats
        pearsonr, p_val = scipy.stats.pearsonr(self.xdata, self.ydata)
        self.Ax.ax.text(0.3, 0.89, r'$r=$ %.2f' % pearsonr,
                        transform=self.Ax.ax.transAxes, fontsize='x-small')

    def best_fit_line(self, print_slope=True):
        """Plot the best fit line to the scatterplot data."""
        import numpy as np
        best_fit = np.polyfit(self.xdata, self.ydata, 1)
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

    def _make_basemap(self):
        import numpy as np
        from mpl_toolkits.basemap import Basemap

        if self.Ax.xlim:
            llcrnrlon, urcrnrlon = self.Ax.xlim
        else:
            llcrnrlon, urcrnrlon = -180, 180
        if self.Ax.ylim:
            llcrnrlat, urcrnrlat = self.Ax.ylim
        else:
            llcrnrlat, urcrnrlat = -90, 90

        self.basemap = Basemap(
            projection=self.Ax.map_proj, resolution=self.Ax.map_res,
            ax=self.Ax.ax, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
            urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat
        )

    def _shiftgrid(self, lon0, datain, lonsin, start=True, cyclic=360.0):
        from mpl_toolkits.basemap import shiftgrid
        return shiftgrid(lon0, datain, lonsin, start=start, cyclic=cyclic)

    def _set_cntr_lvls(self):
        import numpy as np
        self.cntr_lvls = np.linspace(
            self.min_cntr, self.max_cntr, self.num_cntr + 1
        )

    def plot(self):
        import numpy as np
        self._make_basemap()
        self._set_cntr_lvls()
        if self.Ax.shiftgrid_start:
            lon0 = 180
        else:
            lon0 = 180
        plot_data, lons = self._shiftgrid(
            lon0, self.data, self.model.lon,
            start=self.Ax.shiftgrid_start, cyclic=self.Ax.shiftgrid_cyclic
        )
        x, y = self.basemap(*np.meshgrid(lons, self.model.lat))
        self.handle = self.basemap.contourf(
            x, y, plot_data, self.cntr_lvls, extend=self.contourf_extend
        )
        self.handle.set_cmap(self.col_map)

        self.basemap.drawcoastlines(linewidth=0.3)
        self.basemap.drawmapboundary(linewidth=0.3)
