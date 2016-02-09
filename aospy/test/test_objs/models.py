import os

from aospy.model import Model
import runs


am2 = Model(
    name='am2',
    grid_file_paths=(
        ('/archive/Yi.Ming/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie_rerun6.YIM/'
         'pp/atmos/atmos.static.nc'),
        ('/archive/Yi.Ming/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie_rerun6.YIM/'
         'pp/atmos_level/atmos_level.static.nc'),
        ('/archive/Yi.Ming/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie_rerun6.YIM/'
         'pp/atmos_level/ts/monthly/5yr/atmos_level.011601-012012.vcomp.nc'),
        ('/archive/Yi.Ming/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie_rerun6.YIM/'
         'pp/atmos/ts/monthly/100yr/atmos.000101-010012.temp.nc')
    ),
    runs=[runs.test_am2],
    default_runs=[runs.test_am2]
)


idealized_moist = Model(
    name='idealized_moist',
    grid_file_paths=(
        ('/archive/skc/idealized_moist_alb_T42/control_gaussian_T42/'
         'gfdl.ncrc2-default-prod/1x0m720d_32pe/history/00000.1x20days.nc'),
        '/home/skc/inputdata/aquaplanet.land_mask.nc',
    ),
    runs=[runs.test_idealized_moist],
    default_runs=[runs.test_idealized_moist],
)


# Lifted from Spencer Hill's obj library
am3 = Model(
    name='am3',
    grid_file_paths=(
        ('/archive/Spencer.Hill/am3/am3clim_hurrell/gfdl.ncrc2-intel-prod-'
         'openmp/pp/atmos/atmos.static.nc'),
        ('/archive/Spencer.Hill/am3/am3clim_hurrell/gfdl.ncrc2-intel-prod-'
         'openmp/pp/atmos/ts/monthly/1yr/atmos.198101-198112.ucomp.nc'),
        ('/archive/Spencer.Hill/am3/am3clim_hurrell/gfdl.ncrc2-intel-prod-'
         'openmp/pp/atmos_level/ts/monthly/1yr/'
         'atmos_level.198101-198112.ucomp.nc')
    ),
    runs=[runs.test_am3],
    default_runs=[runs.test_am3]
)


hiram = Model(
    name='hiram',
    grid_file_paths=(
        '/archive/Yi.Ming/siena_201211/c180_hiram_clim/'
        'gfdl.ncrc2-default-prod/'
        'pp/atmos/atmos.static.nc',
        '/archive/Yi.Ming/siena_201211/c180_hiram_clim/'
        'gfdl.ncrc2-default-prod/'
        'pp/atmos/ts/monthly/17yr/atmos.197901-199512.ucomp.nc'
    ),
    runs=[runs.test_hiram],
)


root_dir = '/archive/pcmdi/repo/CMIP5/output/'
cesm1_cam5 = Model(
    name='cesm1-cam5',
    description='',
    data_in_dir_struc='gfdl_repo',
    data_in_direc=os.path.realpath(os.path.join(root_dir,
                                                'NSF-DOE-NCAR/CESM1-CAM5')),
    data_in_dur=30,
    data_in_start_date=1979,
    data_in_end_date=2008,
    default_date_range=(1979, 2008),
    runs=[runs.test_amip],
    default_runs=False
)
