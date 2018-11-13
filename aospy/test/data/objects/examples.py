import os

from cftime import DatetimeNoLeap

from aospy import Proj, Model, Run, Var, Region
from aospy.data_loader import NestedDictDataLoader
from aospy.internal_names import LAND_MASK_STR, LON_STR


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
    default_start_date=DatetimeNoLeap(4, 1, 1),
    default_end_date=DatetimeNoLeap(6, 12, 31)
)

example_model = Model(
    name='example_model',
    grid_file_paths=(
        (os.path.join(os.path.split(ROOT_PATH)[0],
                      'netcdf', '00060101.sphum_monthly.nc'),
         os.path.join(os.path.split(ROOT_PATH)[0], 'netcdf',
                      'im.landmask.nc')),
    ),
    runs=[example_run],
    load_grid_data=True,
    grid_attrs={LAND_MASK_STR: 'custom_land_mask', LON_STR: 'custom_lon'}
)

example_proj = Proj(
    'example_proj',
    direc_out=os.path.join(os.path.dirname(__file__), 'test-files'),
    tar_direc_out=os.path.join(os.path.dirname(__file__), 'test-tar-files'),
    models=[example_model]
)

var_not_time_defined = Var(
    name='var_no_time_def',
    def_time=False,
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
    west_bound=0,
    east_bound=360,
    south_bound=-90,
    north_bound=90,
    do_land_mask=False
)

sahel = Region(
    name='sahel',
    description='African Sahel',
    mask_bounds=[(0, 40, 10, 20),
                 (342, 360, 10, 20)],
    do_land_mask=True
)

bk = Var(
    name='bk',
    description=('hybrid coefficients'),
    def_time=False,
    def_vert=True,
    def_lon=False,
    def_lat=False
)

p = Var(
    name='p',
    description=('pressure'),
    def_time=True,
    def_vert=True,
    def_lon=True,
    def_lat=True
)

dp = Var(
    name='dp',
    description=('pressure thicknesses'),
    def_time=True,
    def_vert=True,
    def_lon=True,
    def_lat=True
)
