from aospy.proj import Proj

from . import models


aospy_test = Proj(
    name='aospy_test',
    models=(models.am2, models.idealized_moist, models.am3,
            models.idealized_moist_rad)
)
