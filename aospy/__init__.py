"""aospy: management, analysis, and plotting of gridded climate data."""
from .__config__ import (user_path, LAT_STR, LON_STR, PFULL_STR, PHALF_STR,
                         PLEVEL_STR, TIME_STR, TIME_STR_IDEALIZED)
from . import constants
from .constants import Constant
from . import numerics
from .numerics import FiniteDiff
from . import utils
from . import io
from . import timedate
from .timedate import TimeManager
from . import units
from .units import Units
from . import operator
from .operator import Operator
#from . import spharm_interface # On hold in python3
#from .spharm_interface import SpharmInterface
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
from . import plotting
from .plotting import Fig, Ax, Plot

__all__ = ['Proj', 'Model', 'Run', 'Var', 'Units', 'Constant', 'Region',
           'Fig', 'Ax', 'Plot', 'units', 'calc', 'constants', 'utils', 'io',
           'plotting']
