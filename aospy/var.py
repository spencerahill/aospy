"""var.py: Var class for representing a physical variable in aospy."""
import numpy as np

from .units import Units


class Var(object):
    """Physical variables."""
    def __init__(self, name, alt_names=False, func=False, variables=False,
                 units=False, plot_units='', plot_units_conv=1, domain='atmos',
                 description='', def_time=False, def_vert=False, def_lat=False,
                 def_lon=False, in_nc_grid=False, math_str=False,
                 colormap=False, valid_range=False):
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

        if type(units) is Units:
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
        self.domain = domain
        self.description = description
        self.def_time = def_time
        self.def_vert = def_vert
        self.def_lat = def_lat
        self.def_lon = def_lon
        self.in_nc_grid = in_nc_grid
        if math_str:
            self.math_str = math_str
        if colormap:
            self.colormap = colormap
        if valid_range:
            self.valid_range = valid_range

    def __str__(self):
        return 'Var instance "' + self.name + '"'

    def convert_to_plot_units(self, data):
        """
        Multiply the given data by the plotting units conversion if it exists.
        """
        try:
            return data*self.plot_units_conv
        except:
            pass

    def mask_unphysical(self, data):
        """Mask data array where values are outside physically valid range."""
        try:
            return np.ma.masked_outside(data, self.valid_range[0],
                                        self.valid_range[1])
        except AttributeError:
            return data
