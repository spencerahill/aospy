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
    data_in_direc='/archive/skc/idealized_moist_alb_T42/control_gaussian_T42/'
                  'gfdl.ncrc2-default-prod/1x0m720d_32pe/history',
    data_in_dir_struc='one_dir',
    data_in_files={'20-day': {v: '00000.1x20days.nc'
                              for v in ['olr', 'temp']}},
)
