import imp
from numpy import ma

from . import user_path
from .units import Units

class Var(object):
    def __init__(self, name, alt_names=False, func=False, vars=False,
                 units=False, plot_units='', plot_units_conv=1, domain='atmos',
                 description='', def_time=False, def_vert=False, def_lat=False,
                 def_lon=False, in_nc_grid=False, math_str=False, cmap=False,
                 valid_range=False):
        """Physical variables."""
        self.name = name
        if alt_names:
            self.alt_names = alt_names
        # Identity transform if no function specified.
        if not func:
            self.func = lambda x: x
        else:
            self.func = func
        # `units` kwarg can be `Units` object or string
        
        if type(units) is Units:
            self._Units = units
            for var_attr, units_attr in zip(
                    ('units', 'plot_units', 'plot_units_conv', 'vert_int_units',
                     'vert_int_plot_units', 'vert_int_plot_units_conv'),
                    ('units', 'plot', 'plot_conv', 'vert_int',
                     'vert_int_plot', 'vert_int_plot_conv')
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
        if cmap:
            self.cmap = cmap
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
            return ma.masked_outside(data, self.valid_range[0],
                                     self.valid_range[1])
        except AttributeError:
            return data

variables = imp.load_source(
    'variables', (user_path + '/variables/__init__.py').replace('//','/')
)

def var_inst(var):
    """Convert string of an aospy.var name to an aospy.var instance."""
    if type(var) is Var:
        var_out = var
    elif type(var) is str:
        try:
            var_out = getattr(variables, var)
        except AttributeError:
            raise AttributeError('Not a recognized aospy.Var name: %s'
                                 % var)
    elif type(var) in (list, tuple):
        var_out = [var_inst(v) for v in var]
        if type(var) is tuple:
            var_out = tuple(var_out)
    return var_out
