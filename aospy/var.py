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
        """Create Var object."""
        self.name = name
        if alt_names:
            self.alt_names = alt_names
            self.names = tuple([name] + list(alt_names))
        else:
            self.names = tuple([name])

        if not func:
            self.func = lambda x: x
            self.variables = False
        else:
            self.func = func
            self.variables = variables

        if not isinstance(units, Units):
            self.units = Units(units=units)
        else:
            self.units = units

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

    __repr__ = __str__

    def to_plot_units(self, data, vert_int=False):
        """
        Multiply the given data by the plotting units conversion if it exists.
        """
        if vert_int:
            return data*self.units.vert_int_plot_units_conv
        else:
            return data*self.units.plot_units_conv

    def mask_unphysical(self, data):
        """Mask data array where values are outside physically valid range."""
        if not self.valid_range:
            return data
        else:
            return np.ma.masked_outside(data, np.min(self.valid_range),
                                        np.max(self.valid_range))
