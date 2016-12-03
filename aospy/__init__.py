"""aospy: management, and analysis of gridded climate data."""
from .__config__ import (user_path, LAT_STR, LON_STR, PFULL_STR, PHALF_STR,
                         PLEVEL_STR, TIME_STR, TIME_STR_IDEALIZED, YEAR_STR)
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
from . import find_obj
from .find_obj import to_iterable, to_proj, to_model, to_run, to_var, to_region

__all__ = ['Proj', 'Model', 'Run', 'Var', 'Units', 'Constant', 'Region',
           'units', 'calc', 'constants', 'utils']
