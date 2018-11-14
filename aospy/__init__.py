"""aospy: management, and analysis of gridded climate data."""
from .__config__ import user_path
from . import _constants
from . import utils
from . import var
from .var import Var
from . import region
from .region import Region
from . import run
from .run import Run
from . import model
from .model import Model
from . import proj
from .proj import Proj
from . import calc
from .calc import Calc
from .automate import submit_mult_calcs
from . import examples

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = ['user_path', '_constants', 'utils', 'var', 'Var', 'region',
           'Region', 'run', 'Run', 'model', 'Model', 'proj', 'Proj',
           'calc', 'Calc', 'submit_mult_calcs', 'examples']
