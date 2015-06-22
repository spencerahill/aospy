import imp
from .__config__ import user_path

class Units(object):
    vint_str = r'kg m$^{-2}$'
    def __init__(self, units='', plot=False, plot_conv=1., vert_int=False,
                 vert_int_plot=False, vert_int_plot_conv=False):
        """String representation of physical units and conversion methods.""" 
        self.units = units
        self.plot = plot if plot else units
        self.plot_conv = plot_conv
        if vert_int:
            self.vert_int = vert_int
        else:
            self.vert_int = ' '.join(
                [Units.vint_str, units]).replace('  ', ' ')
        if vert_int_plot:
            self.vert_int_plot = vert_int_plot
        else:
            self.vert_int_plot = ' '.join(
                [Units.vint_str, self.plot]).replace('  ', ' ')
        if vert_int_plot_conv:
            self.vert_int_plot_conv = vert_int_plot_conv
        else:
            self.vert_int_plot_conv = plot_conv

units = imp.load_source(
    'units', (user_path + '/units/__init__.py').replace('//','/')
)
