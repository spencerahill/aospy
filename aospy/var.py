"""var.py: Var class for representing a physical variable in aospy."""
import numpy as np

from .units import Units


class Var(object):
    """Physical variables."""
    def __init__(self, name, alt_names=False, func=False, variables=False,
                 units=False, plot_units='', plot_units_conv=1, domain='atmos',
                 description='', def_time=False, def_vert=False, def_lat=False,
                 def_lon=False, in_nc_grid=False, math_str=False,
                 colormap='RdBu_r', valid_range=False):
        self.name = name
        if alt_names:
            self.alt_names = alt_names

        # Identity transform if no function specified.
        if not func:
            self.func = lambda x: x
            self.variables = False
        else:
            self.func = func
            self.variables = variables

        # `units` kwarg can be `Units` object or string
        if isinstance(units, Units):
            self._Units = units
            for var_attr, units_attr in zip(
                    ('units', 'plot_units', 'plot_units_conv',
                     'vert_int_units', 'vert_int_plot_units',
                     'vert_int_plot_units_conv'),
                    ('units', 'plot', 'plot_conv', 'vert_int', 'vert_int_plot',
                     'vert_int_plot_conv')
            ):
                setattr(self, var_attr, getattr(units, units_attr))
        else:
            self.units = units
            self.plot_units = plot_units
            self.plot_units_conv = plot_units_conv
            self.vert_int_plot_conv = 1
        self.domain = domain
        self.description = description
        self.def_time = def_time
        self.def_vert = def_vert
        self.def_lat = def_lat
        self.def_lon = def_lon
        self.in_nc_grid = in_nc_grid
        self.math_str = math_str
        self.colormap = colormap
        self.valid_range = valid_range

    def __str__(self):
        return 'Var instance "' + self.name + '"'

    def to_plot_units(self, data, vert_int=False):
        """
        Multiply the given data by the plotting units conversion if it exists.
        """
        if vert_int:
            conv_factor = self.vert_int_plot_conv
        else:
            conv_factor = self.plot_units_conv
        return data*conv_factor

    def mask_unphysical(self, data):
        """Mask data array where values are outside physically valid range."""
        if not self.valid_range:
            return data
        else:
            return np.ma.masked_outside(data, np.min(self.valid_range),
                                        np.max(self.valid_range))
