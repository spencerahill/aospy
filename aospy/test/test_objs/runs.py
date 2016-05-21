from aospy.run import Run


test_am2 = Run(
    name='test_am2',
    description=(
        'Preindustrial control simulation.'
    ),
    data_in_direc=('/archive/Yi.Ming/sm2.1_fixed/'
                   'SM2.1U_Control-1860_lm2_aie_rerun6.YIM/pp'),
    data_in_dur=5,
    data_in_start_date='0001-01-01',
    data_in_end_date='0080-12-31',
    default_date_range=('0021-01-01', '0080-12-31')
)


test_idealized_moist = Run(
    name='test_idealized_moist',
    description=(
        'Control case at T42 spectral resolution'
    ),
    data_in_direc='/archive/skc/idealized_moist_alb_T42/control_alb_T42/'
                  'gfdl.ncrc2-default-prod/1x0m720d_32pe/history',
    data_in_dir_struc='one_dir',
    data_in_files={'20-day': {v: '00000.1x20days.nc'
                              for v in ['t_surf', 'temp', 'ps']},
                   '3-hourly': {v: ['{}.8xday.nc'.format(v)]
                                for v in ['ps', 'temp']}},
)


test_idealized_moist_rad = Run(
    name='test_idealized_moist_rad',
    description=(
        'Control simulation of idealized moist run with realistic'
        'radiative transfer.'
    ),
    data_in_direc=('/home/skc/archive/imr_skc/'
                   'control/gfdl.ncrc3-default-repro/1/'
                   'history'),
    default_date_range=('0003-01-01', '0006-12-31'),
    data_in_dir_struc='one_dir',
    data_in_files={'monthly': {v: ['00010101.atmos_month.nc',
                                   '00020101.atmos_month.nc',
                                   '00030101.atmos_month.nc',
                                   '00040101.atmos_month.nc',
                                   '00050101.atmos_month.nc',
                                   '00060101.atmos_month.nc']
                               for v in ['t_surf', 'temp', 'ps']},
                   '3-hourly': {v: ['000{}0101.atmos_8xday.{}.nc'.format(y, v)
                                    for y in range(1, 7)]
                                for v in ['temp', 'ps']}}
)


test_am3 = Run(
    name='hurrell_cont',
    description='am3_hc from Spencer Hills obj library',
    data_in_direc=('/archive/Spencer.Hill/am3/am3clim_hurrell/'
                   'gfdl.ncrc2-intel-prod-openmp/pp'),
    data_in_dur=1,
    data_in_start_date=1980,
    data_in_end_date=2010,
    default_date_range=(1981, 2010)
)


test_hiram = Run(
    name='cont',
    description=(
        '1981-2000 HadISST climatological annual cycle of SSTs '
        'and sea ice repeated annually, with PD atmospheric composition.'
        'hiram_cont from Spencer Hills obj library'
    ),
    data_in_dur=17,
    data_in_start_date=1979,
    data_in_end_date=1995,
    default_date_range=(1979, 1995),
    data_in_direc=('/archive/Yi.Ming/siena_201211/c180_hiram_clim/'
                   'gfdl.ncrc2-default-prod/pp')
)


test_amip = Run(
    name='amip',
    description=('Atmosphere only'
                 'amip from Spencer Hills obj library'),
    data_in_direc='mon/atmos/Amon/r1i1p1',
    data_in_dir_struc='gfdl_repo',
    default_date_range=(1979, 2008),
)
