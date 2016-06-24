"""var.py: Var class for representing a physical variable in aospy."""
import numpy as np

from .units import Units


class Var(object):
    """Physical variables."""
    def __init__(self, name, alt_names=False, func=False, variables=False,
                 units='', plot_units='', plot_units_conv=1, domain='atmos',
                 description='', def_time=False, def_vert=False, def_lat=False,
                 def_lon=False, in_nc_grid=False, math_str=False,
                 colormap='RdBu_r', valid_range=False,
                 func_input_dtype='DataArray'):
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
            self.func_input_dtype = None
        else:
            self.func = func
            self.variables = variables
        self.func_input_dtype = func_input_dtype

        if not isinstance(units, Units):
            self.units = Units(units=units)
        else:
            self.units = units

        if not description:
            try:
                self.description = self.func.func_doc
            except AttributeError:
                self.description = description
        else:
            self.description = description
        self.__doc__ = self.description

        self.domain = domain
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

    def to_plot_units(self, data, dtype_vert=False):
        """
        Multiply the given data by the plotting units conversion if it exists.
        """
        if dtype_vert == 'vert_av' or not dtype_vert:
            conv_factor = self.units.plot_units_conv
        elif dtype_vert == ('vert_int'):
            conv_factor = self.units.vert_int_plot_units_conv
        else:
            raise ValueError("dtype_vert value `{0}` not recognized.  Only "
                             "bool(dtype_vert) = False, 'vert_av', and "
                             "'vert_int' supported.".format(dtype_vert))
        if isinstance(data, dict):
            return {key: val*conv_factor for key, val in data.items()}
        return data*conv_factor

    def mask_unphysical(self, data):
        """Mask data array where values are outside physically valid range."""
        if not self.valid_range:
            return data
        else:
            return np.ma.masked_outside(data, np.min(self.valid_range),
                                        np.max(self.valid_range))
