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
