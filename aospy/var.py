"""Funcionality for representing a physical variable in aospy."""
import numpy as np


class Var(object):
    """An object representing a physical quantity to be computed.

    Attributes
    ----------
    name : str
        The variable's name
    alt_names : tuple of strings
        All other names that the variable may be referred to in the input data
    names : tuple of strings
        The combination of `name` and `alt_names`
    description : str
        A description of the variable
    func : function
        The function with which to compute the variable
    variables : sequence of aospy.Var objects
        The variables passed to `func` to compute it
    units : str
        The variable's physical units
    domain : str
        The physical domain of the variable, e.g. 'atmos', 'ocean', or 'land'
    def_time, def_vert, def_lat, def_lon : bool
        Whether the variable is defined, respectively, in time, vertically, in
        latitude, and in longitude
    math_str : str
        The mathematical representation of the variable
    colormap : str
        The name of the default colormap to be used in plots of this variable
    valid_range : length-2 tuple
        The range of values outside which to flag as unphysical/erroneous

    """

    def __init__(self, name, alt_names=None, func=None, variables=None,
                 units='', plot_units='', plot_units_conv=1, domain='atmos',
                 description='', def_time=False, def_vert=False, def_lat=False,
                 def_lon=False, math_str=False, colormap='RdBu_r',
                 valid_range=None):
        """Instantiate a Var object.

        Parameters
        ----------
        name : str
            The variable's name
        alt_names : tuple of strings
            All other names that the variable might be referred to in any input
            data.  Each of these should be unique to this variable in order to
            avoid loading the wrong quantity.
        description : str
            A description of the variable
        func : function
            The function with which to compute the variable
        variables : sequence of aospy.Var objects
            The variables passed to `func` to compute it.  Order matters:
            whenever calculations are performed to generate data corresponding
            to this Var, the data corresponding to the elements of `variables`
            will be passed to `self.function` in the same order.
        units : str
            The variable's physical units
        domain : str
            The physical domain of the variable, e.g. 'atmos', 'ocean', or
            'land'.  This is only used by aospy by some types of `DataLoader`,
            including `GFDLDataLoader`.
        def_time, def_vert, def_lat, def_lon : bool
            Whether the variable is defined, respectively, in time, vertically,
            in latitude, and in longitude
        math_str : str
            The mathematical representation of the variable.  This is typically
            a raw string of LaTeX math-mode, e.g. r'$T_\mathrm{sfc}$' for
            surface temperature.
        colormap : str
            (Currently not used by aospy) The name of the default colormap to
            be used in plots of this variable.
        valid_range : length-2 tuple
            The range of values outside which to flag as unphysical/erroneous
        """  # noqa: W605
        self.name = name
        if alt_names is None:
            self.names = tuple([name])
        else:
            self.alt_names = alt_names
            self.names = tuple([name] + list(alt_names))

        if func is None:
            self.func = lambda x: x
            self.variables = None
        else:
            self.func = func
            self.variables = variables

        self.units = units

        if not description:
            if self.func.__doc__ is None:
                self.description = ''
            else:
                self.description = self.func.__doc__
        else:
            self.description = description

        self.domain = domain
        self.def_time = def_time
        self.def_vert = def_vert
        self.def_lat = def_lat
        self.def_lon = def_lon
        self.math_str = math_str
        self.colormap = colormap
        self.valid_range = valid_range

    def __str__(self):
        return 'Var instance "' + self.name + '"'

    __repr__ = __str__

    def to_plot_units(self, data, dtype_vert=False):
        """Convert the given data to plotting units."""
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
