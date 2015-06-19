from .__config__ import user_path
from . import constants
from . import utils
from . import io
from . import proj
from .proj import Proj
from . import model
from .model import Model
from . import run
from .run import Run
from . import var
from .var import Var
from . import region
from .region import Region
from . import units
from .units import Units
from . import calc
from .calc import Calc
from . import plotting
from .plotting import Fig, Ax, Plot

# import imp as _imp
# calcs = _imp.load_source('calcs', '/home/s1h/aospy/calcs/__init__.py')
# import calcs as cc
# del _imp

__all__ = ['Proj', 'Model', 'Run', 'Var', 'Region', 'Fig', 'Ax',
           'Plot','core', 'units', 'calc', 'constants', 'utils',
           'io', 'plotting']
# from os import getenv
# aospy_path = '/'.join([getenv('HOME'), 'aospy']).replace('//','/')

# def _load_user_data(name):
#     """Load user data from aospy_path for given module name.

#     File must be located in the `aospy_path` directory and be the same name
#     as the desired aospy module subpackage, namely one of `regions`, `calcs`,
#     `variables`, and `projects`.
#     """
#     import imp
#     return imp.load_source(
#         name, '/'.join([aospy_path, name, '__init__.py']).replace('//','/')
#     )

# calcs = _load_user_data('calcs')
# regions = _load_user_data('regions')
# variables = _load_user_data('variables')
# from main import main
