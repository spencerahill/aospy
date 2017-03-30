"""aospy: management, and analysis of gridded climate data."""
from .__config__ import user_path
from . import constants
from .constants import Constant
from . import utils
from . import units
from .units import Units
from . import operator
from .operator import Operator
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
from .calc import CalcInterface, Calc
from .automate import submit_mult_calcs
from . import examples

__all__ = ['Proj', 'Model', 'Run', 'Var', 'Units', 'Constant', 'Region',
           'units', 'calc', 'constants', 'utils']
