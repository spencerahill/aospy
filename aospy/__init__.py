"""aospy: management, analysis, and plotting of gridded climate data."""
from .__config__ import user_path
from . import constants
from .constants import Constant
from . import utils
from . import io
from . import units
from .units import Units
from . import numerics
from .numerics import FiniteDiff
from . import windspharm_interface
from .windspharm_interface import WindspharmInterface
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
from . import plotting
from .plotting import Fig, Ax, Plot

__all__ = ['Proj', 'Model', 'Run', 'Var', 'Units', 'Constant', 'Region',
           'Fig', 'Ax', 'Plot', 'units', 'calc', 'constants', 'utils', 'io',
           'plotting']
