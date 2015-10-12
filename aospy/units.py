"""units.py: Units class for representing physical units, e.g. meters."""


class Units(object):
    vert_int_str = r'kg m$^{-2}$'

    def __init__(self, units='', plot_units=False, plot_units_conv=1.,
                 vert_int_units=False, vert_int_plot_units=False,
                 vert_int_plot_units_conv=False):
        """String representation of physical units and conversion methods."""
        self.units = units
        if plot_units:
            self.plot_units = plot_units
        else:
            self.plot_units = units
        self.plot_units_conv = plot_units_conv
        if vert_int_units:
            self.vert_int_units = vert_int_units
        else:
            self.vert_int_units = ' '.join(
                [Units.vert_int_str, units]).replace('  ', ' ')
        if vert_int_plot_units:
            self.vert_int_plot_units = vert_int_plot_units
        else:
            self.vert_int_plot_units = ' '.join(
                [Units.vert_int_str, self.plot_units]).replace('  ', ' ')
        if vert_int_plot_units_conv:
            self.vert_int_plot_units_conv = vert_int_plot_units_conv
        else:
            self.vert_int_plot_units_conv = plot_units_conv
