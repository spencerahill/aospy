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
