from datetime import datetime
import os

from aospy import Proj, Model, Run, Var, Region
from aospy.data_loader import NestedDictDataLoader

ROOT_PATH = os.path.dirname(__file__)


def total_precipitation(convection_rain, condensation_rain):
    return convection_rain + condensation_rain

precip_files = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                            '000[4-6]0101.precip_monthly.nc')
sphum_files = os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                           '00060101.sphum_monthly.nc')
file_map = {'monthly': {'condensation_rain': precip_files,
                        'convection_rain': precip_files,
                        'sphum': sphum_files,
                        'ps': sphum_files}}
example_run = Run(
    name='example_run',
    description=(
        'Control simulation of the idealized moist model'
    ),
    data_loader=NestedDictDataLoader(file_map),
    default_start_date=datetime(4, 1, 1),
    default_end_date=datetime(6, 12, 31)
)

example_model = Model(
    name='example_model',
    grid_file_paths=(
        (os.path.join(os.path.split(ROOT_PATH)[0],
                      'netcdf', '00060101.sphum_monthly.nc'),
         os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                      'im.landmask.nc')),
    ),
    runs=[example_run]
)

example_proj = Proj(
    'example_proj',
    direc_out=os.path.join(os.path.dirname(__file__), 'test-files'),
    tar_direc_out=os.path.join(os.path.dirname(__file__), 'test-tar-files'),
    models=[example_model]
)

condensation_rain = Var(
    name='condensation_rain',
    alt_names=('prec_ls',),
    def_time=True,
    description=('condensation rain'),
)

convection_rain = Var(
    name='convection_rain',
    alt_names=('prec_conv',),
    def_time=True,
    description=('convection rain'),
)

precip = Var(
    name='total_precipitation',
    def_time=True,
    description=('total precipitation rate'),
    func=total_precipitation,
    variables=(convection_rain, condensation_rain)
)

ps = Var(
    name='ps',
    def_time=True,
    description=('surface pressure'),
    def_vert=False
)

sphum = Var(
    name='sphum',
    def_time=True,
    description=('specific humidity'),
    def_vert=True
)

globe = Region(
    name='globe',
    description='Entire globe',
    lat_bounds=(-90, 90),
    lon_bounds=(0, 360),
    do_land_mask=False
)

sahel = Region(
    name='sahel',
    description='African Sahel',
    mask_bounds=[((10, 20), (0, 40)), ((10, 20), (342, 360))],
    do_land_mask=True
)
