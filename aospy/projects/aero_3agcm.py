def aero_3agcm():
    from aospy import Proj, Model, Run, _set_named_attr_dict
    from aospy.regions import *
    from aospy.variables import master_vars_list
    
    ### Proj ###
    a3gcm = Proj(
        name='aero_3agcm',
        vars=master_vars_list,
        direc_out='/archive/s1h/aero_3agcm/',
        nc_dir_struc='gfdl',
    )
    ### Runs ###
    ## AM2 runs.
    # AM2 cont
    am2_cont = Run(
        proj=a3gcm, 
        name='cont',
        description=(
            '1981-2000 HadISST climatological annual cycle of SSTs and sea '
            'ice repeated annually, with PD atmospheric composition.'
        ),
        direc_nc=('/archive/yim/siena_201203/m45_am2p14_1990/'
                  'gfdl.ncrc2-intel-prod/pp'),
        nc_dur=16,
        nc_start_yr=1983,
        nc_end_yr=1998,
        default_yr_range=(1983, 1998)
    )
    # AM2 aero
    # am2_aero = Run(
    #     proj=a3gcm,
    #     name='aero',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_aero/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 aero_tm
    # am2_atm = Run(
    #     proj=a3gcm,
    #     name='aero_tm',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual tropical mean equilibrium SST anomaly from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_aero_trop_mean/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 aero_mtm
    # am2_amtm = Run(
    #     proj=a3gcm,
    #     name='aero_mtm',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, subtracting annual tropical mean equilibrium SST anomaly from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_aero_m_trop_mean_fixed2/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 aero_pac
    # am2_apac = Run(
    #     proj=a3gcm,
    #     name='aero_pac',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Pacific Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_aero_pac/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 aero_atl
    # am2_aatl = Run(
    #     proj=a3gcm,
    #     name='aero_atl',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Atlantic Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_aero_atl/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 aero_ind
    # am2_aind = Run(
    #     proj=a3gcm,
    #     name='aero_ind',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Indian Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_aero_ind/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 gas
    # am2_gas = Run(
    #     proj=a3gcm,
    #     name='gas',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies from a PI-to-PD WMGG and ozone simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_gas/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 gas_tm
    # am2_gtm = Run(
    #     proj=a3gcm,
    #     name='gas_tm',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual tropical mean equilibrium SST anomaly from a PI-to-PD WMGG and ozone simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_gas_trop_mean/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 gas_mtm
    # am2_gmtm = Run(
    #     proj=a3gcm,
    #     name='gas_mtm',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies minus their annual tropical mean from a PI-to-PD WMGG and ozone simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201203/m45_am2p14_1990_clim_gas_m_trop_mean/gfdl.ncrc2-intel-prod/pp',
    #     nc_dur=16,
    #     nc_start_yr=1983,
    #     nc_end_yr=1998,
    #     default_yr_range=(1983, 1998)
    # )
    # # AM2 noT
    # am2_noT = Run(
    #     proj=a3gcm,
    #     name='noTok',
    #     description='',
    #     direc_nc='/archive/miz/GCM/miz_cess_noT/cess/am2_cess/pp/',
    #     nc_dur=5,
    #     nc_start_yr=1983,
    #     nc_end_yr=1987,
    #     default_yr_range=(1983, 1987)
    # )
    # # AM2 noT_p2K
    # am2_noT_p2K = Run(
    #     proj=a3gcm,
    #     name='noTok_p2K',
    #     description='',
    #     direc_nc='/archive/miz/GCM/miz_cess_noT/cess+2/am2_cess+2/pp/',
    #     nc_dur=5,
    #     nc_start_yr=1983,
    #     nc_end_yr=1987,
    #     default_yr_range=(1983, 1987)
    # )
    # AM2 amip
    am2_amip = Run(
        proj=a3gcm,
        name='amip',
        description='',
        direc_nc='/archive/fjz/AM2.1_1870-2004/AM2.1_1870-2004-HGlob-SST' + \
                 '-ICE-AllForc_B1-_B10_ens/pp/',
        nc_dur=130,
        nc_start_yr=1870,
        nc_end_yr=1999,
        default_yr_range=(1870, 1999)
    )
    # AM2 extended control
    am2_reyoi_cont=Run(
        proj=a3gcm,
        name='reyoi_cont',
        tags=['reyoi', 'cont'],
        description='PI atmos and Reynolds OI climatological SSTs',
        direc_nc='/archive/s1h/am2/am2clim_reyoi/gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -0.25K
    am2_reyoi_m0p25 = Run(
        proj=a3gcm,
        name='reyoi-0.25K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-0p25K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )    
    # AM2 -0.5K
    am2_reyoi_m0p5=Run(
        proj=a3gcm,
        name='reyoi-0.5K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-0p5K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=15,
        nc_start_yr=1983,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -1K
    am2_reyoi_m1 = Run(
        proj=a3gcm,
        name='reyoi-1K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-1K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=15,
        nc_start_yr=1983,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -1.5K
    am2_reyoi_m1p5 = Run(
        proj=a3gcm,
        name='reyoi-1.5K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-1p5K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )    
    # AM2 -2K
    am2_reyoi_m2 = Run(
        proj=a3gcm,
        name='reyoi-2K',
        tags=['reyoi', '-2K'],
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-2K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=15,
        nc_start_yr=1983,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -3K
    am2_reyoi_m3 = Run(
        proj=a3gcm,
        name='reyoi-3K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-3K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -4K
    am2_reyoi_m4 = Run(
        proj=a3gcm,
        name='reyoi-4K',
        description='',
        direc_nc=('/archive/s1h/am2/am2clim_reyoi-4K/'
                  'gfdl.ncrc2-default-prod/pp/'),
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -6K
    am2_reyoi_m6 = Run(
        proj=a3gcm,
        name='reyoi-6K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-6K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -8K
    am2_reyoi_m8 = Run(
        proj=a3gcm,
        name='reyoi-8K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-8K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 -10K
    am2_reyoi_m10 = Run(
        proj=a3gcm,
        name='reyoi-10K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi-10K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +0.25K
    am2_reyoi_p0p25 = Run(
        proj=a3gcm,
        name='reyoi+0.25K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+0p25K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )    
    # AM2 +0.5K
    am2_reyoi_p0p5 = Run(
        proj=a3gcm,
        name='reyoi+0.5K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+0p5K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=15,
        nc_start_yr=1983,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +1K
    am2_reyoi_p1 = Run(
        proj=a3gcm,
        name='reyoi+1K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+1K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=15,
        nc_start_yr=1983,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +1.5K
    am2_reyoi_p1p5 = Run(
        proj=a3gcm,
        name='reyoi+1.5K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+1p5K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +2K
    am2_reyoi_p2 = Run(
        proj=a3gcm,
        name='reyoi+2K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+2K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +3K
    am2_reyoi_p3 = Run(
        proj=a3gcm,
        name='reyoi+3K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+3K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +4K
    am2_reyoi_p4 = Run(
        proj=a3gcm,
        name='reyoi+4K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+4K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +6K
    am2_reyoi_p6 = Run(
        proj=a3gcm,
        name='reyoi+6K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+6K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +8K
    am2_reyoi_p8 = Run(
        proj=a3gcm,
        name='reyoi+8K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+8K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 +10K
    am2_reyoi_p10 = Run(
        proj=a3gcm,
        name='reyoi+10K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi+10K/' + \
                 'gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 WPWP+2K
    am2_reyoi_wpwp_p2 = Run(
        proj=a3gcm,
        name='reyoi_wpwp+2K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi_wpwp+2K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 WPWP-2K
    am2_reyoi_wpwp_m2=Run(
        proj=a3gcm,
        name='reyoi_wpwp-2K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi_wpwp-2K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=30,
        nc_start_yr=1983,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 WPWP-2K
    am2_reyoi_wpwp_m2=Run(
        proj=a3gcm,
        name='reyoi_wpwp-2K',
        description='',
        direc_nc='/archive/s1h/am2/am2clim_reyoi_wpwp-2K/gfdl.ncrc2-default-prod/pp/',
        nc_dur=30,
        nc_start_yr=1983,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012)
    )
    # AM2 locked clouds control
    am2_cld_lock_cont=Run(
        proj=a3gcm,
        name='cld_lock_cont',
        description='',
        direc_nc=('/archive/yim/quickstart/m45_am2p14_1990_nocre_1995/'
                  'gfdl.ncrc2-default-prod/pp/'),
        nc_dur=16,
        nc_start_yr=1983,
        nc_end_yr=1998,
        default_yr_range=(1983, 1998)
    )
    # AM2 locked clouds +2K clouds control SST
    am2_cld_lock_cld=Run(
        proj=a3gcm,
        name='cld_lock+2Kcld',
        description='',
        direc_nc=('/archive/yim/quickstart/m45_am2p14_1990_nocre_1995_p2K_fix2/'
                  'gfdl.ncrc2-default-prod/pp/'),
        nc_dur=16,
        nc_start_yr=1983,
        nc_end_yr=1998,
        default_yr_range=(1983, 1998)
    )
    # AM2 locked clouds +2K SST control clouds
    am2_cld_lock_sst=Run(
        proj=a3gcm,
        name='cld_lock+2Ksst',
        description='',
        direc_nc=('/archive/yim/quickstart/m45_am2p14_1990_nocre_1995_p2K_fix1/'
                  'gfdl.ncrc2-default-prod/pp/'),
        nc_dur=16,
        nc_start_yr=1983,
        nc_end_yr=1998,
        default_yr_range=(1983, 1998)
    )
    # AM2 locked clouds +2K both
    am2_cld_lock_p2=Run(
        proj=a3gcm,
        name='cld_lock+2K',
        description='',
        direc_nc=('/archive/yim/quickstart/m45_am2p14_1990_nocre_1995_p2K/'
                  'gfdl.ncrc2-default-prod/pp/'),
        nc_dur=16,
        nc_start_yr=1983,
        nc_end_yr=1998,
        default_yr_range=(1983, 1998)
    )
    # # AM2 Hurrell SSTs control
    # am2_hurrell_cont = Run(
    #     proj=a3gcm,
    #     name='hurrell_cont',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_hurrell/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 Hurrell SSTs +2 K
    # am2_hurrell_p2 = Run(
    #     proj=a3gcm,
    #     name='hurrell+2K',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_hurrell+2K/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 ReynoldsEOF SSTs control
    # am2_reynolds = Run(
    #     proj=a3gcm,
    #     name='reynolds_cont',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_reynoldsEOF/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 ReynoldsEOF SSTs +2K
    # am2_reynolds_p2 = Run(
    #     proj=a3gcm,
    #     name='reynolds+2K',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_reynoldsEOF+2K/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 AMIP1 SSTs control
    # am2_amip1 = Run(
    #     proj=a3gcm,
    #     name='amip1_cont',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_amip1/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 AMIP1 SSTs +2K
    # am2_amip1_p2 = Run(
    #     proj=a3gcm,
    #     name='amip1+2K',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_amip1+2K/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 cld_seed_all +2K
    # am2_cld_seed_all_p2 = Run(
    #     proj=a3gcm,
    #     name='cld_seed_all+2K',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_reyoi_cld_seed_all+2K/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 cld_seed_np +2K
    # am2_cld_seed_np_p2 = Run(
    #     proj=a3gcm,
    #     name='cld_seed_np+2K',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_reyoi_cld_seed_np+2K/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 cld_seed_all +2K
    # am2_cld_seed_sp_p2 = Run(
    #     proj=a3gcm,
    #     name='cld_seed_sp+2K',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_reyoi_cld_seed_sp+2K/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )
    # # AM2 cld_seed_sa +2K
    # am2_cld_seed_sa_p2 = Run(
    #     proj=a3gcm,
    #     name='cld_seed_sa+2K',
    #     description='',
    #     direc_nc='/archive/s1h/am2/am2clim_reyoi_cld_seed_sa+2K/' + \
    #              'gfdl.ncrc2-default-prod/pp/',
    #     nc_dur=1,
    #     nc_start_yr=1983,
    #     nc_end_yr=2012,
    #     default_yr_range=(1983, 2012)
    # )

    # ## AM3 runs
    # # AM3 cont
    # am3_cont=Run(proj=a3gcm,
    #                name='cont',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, with PD atmospheric composition.',
    #                direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10/gfdl.ncrc2-intel-prod/pp'
    #                )
    # # AM3 aero
    # am3_aero=Run(proj=a3gcm,
    #                name='aero',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_aero/gfdl.ncrc2-intel-prod/pp'
    #                )
    # # AM3 aero_tm
    # am3_atm=Run(proj=a3gcm,
    #               name='aero_tm',
    #               description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual tropical mean equilibrium SST anomaly from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #               direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_aero_trop_mean/gfdl.ncrc2-intel-prod/pp'
    #               )
    # # AM3 aero_mtm
    # am3_amtm=Run(proj=a3gcm,
    #                name='aero_mtm',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, subtracting annual tropical mean equilibrium SST anomaly from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_aero_m_trop_mean_fixed/gfdl.ncrc2-intel-prod-openmp/pp'
    #                )
    # # AM3 aero_pac
    # am3_apac=Run(proj=a3gcm,
    #                name='aero_pac',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Pacific Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_aero_pac/gfdl.ncrc2-intel-prod/pp'
    #                )
    # # AM3 aero_atl
    # am3_aatl=Run(proj=a3gcm,
    #                name='aero_atl',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Atlantic Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_aero_atl/gfdl.ncrc2-intel-prod/pp'
    #                )
    # # AM3 aero_ind
    # am3_aind=Run(proj=a3gcm,
    #                name='aero_ind',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Indian Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_aero_ind/gfdl.ncrc2-intel-prod/pp'
    #                )
    # # AM3 gas
    # am3_gas=Run(
    #     proj=a3gcm,
    #     name='gas',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies from a PI-to-PD WMGG and ozone simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_gas/gfdl.ncrc2-intel-prod/pp'
    # )
    # # AM3 gas_tm
    # am3_gtm=Run(
    #     proj=a3gcm,
    #     name='gas_tm',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual tropical mean equilibrium SST anomaly from a PI-to-PD WMGG and ozone simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_gas_trop_mean/gfdl.ncrc2-intel-prod/pp'
    # )
    # # AM3 gas_mtm
    # am3_gmtm=Run(
    #     proj=a3gcm,
    #     name='gas_mtm',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies minus their annual tropical mean from a PI-to-PD WMGG and ozone simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/fms/siena_201211/c48L48_am3p10_gas_m_trop_mean/gfdl.ncrc2-intel-prod/pp'
    # )
    # AM3 AMIP
    am3_amip=Run(
        proj=a3gcm,
        name='amip',
        description='',
        ens_mem_prefix='/archive/lwh/fms/riga_201104/c48L48_am3p9_',
        ens_mem_ext=['ext', 'ext2', 'ext3'],
        ens_mem_suffix='/gfdl.intel-prod/pp',
        nc_dur=136,
        nc_start_yr=1870,
        nc_end_yr=2005,
        default_yr_range=(1870,2005)
    )
    # AM3 Hurrell control re-run
    am3_hc = Run(
        proj=a3gcm,
        name='hurrell_cont',
        description='',
        direc_nc='/archive/s1h/am3/am3clim_hurrell/' + \
                 'gfdl.ncrc2-intel-prod-openmp/pp',
        nc_dur=1,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell +1K
    am3_hp1k = Run(
        proj=a3gcm,
        name='hurrell+1K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell+1K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell +2K
    am3_hp2k = Run(
        proj=a3gcm,
        name='hurrell+2K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell+2K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=1,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell +4K
    am3_hp4k = Run(
        proj=a3gcm,
        name='hurrell+4K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell+4K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=1,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell +6K
    am3_hp6k = Run(
        proj=a3gcm,
        name='hurrell+6K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell+6K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell +8K
    am3_hp8k = Run(
        proj=a3gcm,
        name='hurrell+8K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell+8K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell +10K
    am3_hp10k = Run(
        proj=a3gcm,
        name='hurrell+10K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell+10K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell -1K
    am3_hm1k = Run(
        proj=a3gcm,
        name='hurrell-1K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell-1K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=1,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell -2K
    am3_hm2k = Run(
        proj=a3gcm,
        name='hurrell-2K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell-2K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=1,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell -4K
    am3_hm4k = Run(
        proj=a3gcm,
        name='hurrell-4K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell-4K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell -6K
    am3_hm6k = Run(
        proj=a3gcm,
        name='hurrell-6K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell-6K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell -8K
    am3_hm8k = Run(
        proj=a3gcm,
        name='hurrell-8K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell-8K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell -10K
    am3_hm10k = Run(
        proj=a3gcm,
        name='hurrell-10K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell-10K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell -15K
    am3_hm15k = Run(
        proj=a3gcm,
        name='hurrell-15K',
        description='',
        direc_nc=('/archive/s1h/am3/am3clim_hurrell-15K/'
                  'gfdl.ncrc2-intel-prod-openmp/pp'),
        nc_dur=31,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    # AM3 Hurrell WPWP +2K
    am3_hwpwp_p2k = Run(
        proj=a3gcm,
        name='hurrell_wpwp+2K',
        description='',
        direc_nc='/archive/s1h/am3/am3clim_hurrell_wpwp+2K/' + \
                 'gfdl.ncrc2-intel-prod-openmp/pp',
        nc_dur=1,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981, 2010)
    )
    ## HiRAM runs
    hiram_cont=Run(
        proj=a3gcm,
        name='cont',
        description=(
            '1981-2000 HadISST climatological annual cycle of SSTs '
            'and sea ice repeated annually, with PD atmospheric composition.'
        ),
        direc_nc=('/archive/yim/siena_201211/c180_hiram_clim/'
                  'gfdl.ncrc2-default-prod/pp')
    )
    # # HiRAM aero
    # hiram_aero=Run(proj=a3gcm,
    #                name='aero',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/siena_201211/c180_hiram_clim_aero/gfdl.ncrc2-default-prod/pp'
    #                )
    # # HiRAM aero_tm
    # hiram_atm=Run(proj=a3gcm,
    #               name='aero_tm',
    #               description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual tropical mean equilibrium SST anomaly from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #               direc_nc='/archive/yim/siena_201211/c180_hiram_clim_aero_trop_mean/gfdl.ncrc2-default-prod/pp'
    #               )
    # # HiRAM aero_mtm
    # hiram_amtm=Run(proj=a3gcm,
    #                name='aero_mtm',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, subtracting annual tropical mean equilibrium SST anomaly from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/siena_201211/c180_hiram_clim_aero_m_trop_mean/gfdl.ncrc2-default-prod/pp'
    #                )
    # # HiRAM aero_pac
    # hiram_apac=Run(proj=a3gcm,
    #                name='aero_pac',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Pacific Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/siena_201211/c180_hiram_clim_aero_pac/gfdl.ncrc2-default-prod/pp'
    #                )
    # # HiRAM aero_atl
    # hiram_aatl=Run(proj=a3gcm,
    #                name='aero_atl',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Atlantic Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/siena_201211/c180_hiram_clim_aero_atl/gfdl.ncrc2-default-prod/pp'
    #                )
    # # HiRAM aero_ind
    # hiram_aind=Run(proj=a3gcm,
    #                name='aero_ind',
    #                description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid in Indian Ocean only with annual cycle of equilibrium SST anomalies from a PI-to-PD aerosols simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #                direc_nc='/archive/yim/siena_201211/c180_hiram_clim_aero_ind/gfdl.ncrc2-default-prod/pp'
    #                )
    # # HiRAM gas
    # hiram_gas=Run(
    #     proj=a3gcm,
    #     name='gas',
    #     description='1981-2000 HadISST climatological annual cycle of SSTs and sea ice repeated annually, overlaid with annual cycle of equilibrium SST anomalies from a PI-to-PD WMGG and ozone simulation of AM2.1 with a mixed layer ocean. PD atmospheric composition.',
    #     direc_nc='/archive/yim/siena_201211/c180_hiram_clim_gas_rerun2/gfdl.ncrc2-default-prod/pp'
    # )
    # HiRAM gas_tm
    hiram_gtm=Run(
        proj=a3gcm,
        name='gas_tm',
        description=(
            '1981-2000 HadISST climatological annual cycle of SSTs and sea '
            'ice repeated annually, overlaid with annual tropical mean '
            'equilibrium SST anomaly from a PI-to-PD WMGG and ozone simulation '
            'of AM2.1 with a mixed layer ocean. PD atmospheric composition.'
            ),
        direc_nc=('/archive/yim/siena_201211/c180_hiram_clim_gas_trop_mean/'
                  'gfdl.ncrc2-default-prod/pp')
    )
    # # HiRAM gas_mtm
    # hiram_gmtm=Run(
    #     proj=a3gcm,
    #     name='gas_mtm',
    #     description=(
    #         '1981-2000 HadISST climatological annual cycle of SSTs and sea ice'
    #         'repeated annually, overlaid with annual cycle of equilibrium SST'
    #         'anomalies minus their annual tropical mean from a PI-to-PD WMGG &'
    #         'ozone simulation of AM2.1 with a mixed layer ocean. PD atmos'
    #         'composition.'
    #     ),
    #     direc_nc=('/archive/yim/siena_201211/c180_hiram_clim_gas_m_trop_mean'
    #               '/gfdl.ncrc2-default-prod/pp')
    # )
    hiram_amip = Run(
        proj=a3gcm,
        name='amip',
        description='',
        ens_mem_prefix='/archive/hrao/ornl/cmip5/c180_hiram_',
        ens_mem_ext=['H1', 'H3'],
        ens_mem_suffix='/pp',
        nc_dur=5,
        nc_start_yr=1979,
        nc_end_yr=2008,
        default_yr_range=(1979, 2008),
        nc_dir_struc='gfdl'
    )
    # ## SM2.1 runs.
    # # SM2.1 cont
    # sm2_cont=Run(
    #     proj=a3gcm,
    #     name='cont',
    #     description='',
    #     direc_nc=('/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie'
    #               '_rerun6.YIM/pp'),
    #     nc_dur=20,
    #     nc_start_yr=1,
    #     nc_end_yr=120
    # )
    # # SM2.1 aero
    # sm2_aero=Run(
    #     proj=a3gcm,
    #     name='aero',
    #     description='',
    #     direc_nc=('/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie2'
    #               '_rerun6.YIM/pp'),
    #     nc_dur=100,
    #     nc_start_yr=1,
    #     nc_end_yr=100
    # )
    # # SM2.1 gas
    # sm2_gas=Run(
    #     proj=a3gcm,
    #     name='gas',
    #     description='',
    #     direc_nc=('/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie3'
    #               '_rerun8.YIM/pp'),
    #     nc_dur=5,
    #     nc_start_yr=1,
    #     nc_end_yr=80
    # )
    # # SM2.1 both
    # sm2_both=Run(
    #     proj=a3gcm,
    #     name='both',
    #     description='',
    #     direc_nc=('/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie4'
    #               '_rerun6.YIM/pp'),
    #     nc_dur=100,
    #     nc_start_yr=1,
    #     nc_end_yr=100
    # )

    ## HiRAM_ming runs
    # HiRAM ming0
    hiram_c48_0=Run(
        proj=a3gcm,
        name='ming0',
        description='',
        direc_nc=('/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0/'
                  'gfdl.ncrc2-intel-prod/pp')
    )
    # HiRAM ming0_p2K
    hiram_c48_0_p2K=Run(
        proj=a3gcm,
        name='ming0_p2K',
        description='',
        direc_nc=('/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0_p2K/'
                  'gfdl.ncrc2-intel-prod/pp')
    )
    # # HiRAM ming1
    # hiram_c48_1=Run(
    #     proj=a3gcm,
    #     name='ming1',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0b/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming1_p2K
    # hiram_c48_1_p2K=Run(
    #     proj=a3gcm,
    #     name='ming1_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0b_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming2
    # hiram_c48_2=Run(
    #     proj=a3gcm,
    #     name='ming2',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0e/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming2_p2K
    # hiram_c48_2_p2K=Run(
    #     proj=a3gcm,
    #     name='ming2_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0e_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming3
    # hiram_c48_3=Run(
    #     proj=a3gcm,
    #     name='ming3',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0f/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming3_p2K
    # hiram_c48_3_p2K=Run(
    #     proj=a3gcm,
    #     name='ming3_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0f_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming4
    # hiram_c48_4=Run(
    #     proj=a3gcm,
    #     name='ming4',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0c/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming4_p2K
    # hiram_c48_4_p2K=Run(
    #     proj=a3gcm,
    #     name='ming4_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0c_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming5
    # hiram_c48_5=Run(
    #     proj=a3gcm,
    #     name='ming5',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X01/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming5_p2K
    # hiram_c48_5_p2K=Run(
    #     proj=a3gcm,
    #     name='ming5_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X01_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming6
    # hiram_c48_6=Run(
    #     proj=a3gcm,
    #     name='ming6',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X02/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming6_p2K
    # hiram_c48_6_p2K = Run(
    #     proj=a3gcm,
    #     name='ming6_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X02_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming7
    # hiram_c48_7 = Run(
    #     proj=a3gcm,
    #     name='ming7',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X03/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming7_p2K
    # hiram_c48_7_p2K = Run(
    #     proj=a3gcm,
    #     name='ming7_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X03_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming8
    # hiram_c48_8 = Run(
    #     proj=a3gcm,
    #     name='ming8',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X04/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # HiRAM ming8_p2K
    # hiram_c48_8_p2K = Run(
    #     proj=a3gcm,
    #     name='ming8_p2K',
    #     description='',
    #     direc_nc='/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X04_p2K/'+\
    #              'gfdl.ncrc2-intel-prod/pp'
    # )

    # ## AM3_c90 runs.
    # # AM3_c90 cont
    # am3c90_cont = Run(
    #     proj=a3gcm,
    #     name='cont',
    #     description='',
    #     direc_nc='/archive/h1g/FMS/siena_201203/c90L48_am3p10_v6_clim/' + \
    #              'gfdl.ncrc2-intel-prod-openmp/pp'
    # )
    # # AM3_c90 p2K
    # am3c90_p2K = Run(
    #     proj=a3gcm,
    #     name='p2K',
    #     description='',
    #     direc_nc='/archive/h1g/FMS/siena_201203/c90L48_am3p10_v6_clim_p2k/' + \
    #              'gfdl.ncrc2-intel-prod-openmp/pp'
    # )

    # ## AM2.5 runs
    # # AM2.5 cont
    # am2p5_cont = Run(
    #     proj=a3gcm,
    #     name='cont',
    #     description='',
    #     direc_nc='/archive/miz/hiramdp/siena_201204/c180l32_am2_C0/' + \
    #              'gfdl.ncrc2-intel-prod/pp'
    # )
    # # AM2.5 p2K
    # am2p5_p2K = Run(
    #     proj=a3gcm,
    #     name='p2K',
    #     description='',
    #     direc_nc='/archive/miz/hiramdp/siena_201204/c180l32_am2_C0_p2K/gfdl.ncrc2-intel-prod/pp'
    # )

    # ## AM4 prototypes runs
    # # AM4a1 control
    am4_a1c = Run(
        proj=a3gcm,
        name='cont',
        description='',
        direc_nc='/archive/Ming.Zhao/awg/tikal_201403/c96L48_am4a1_' + \
                 '2000climo_highsen1/gfdl.ncrc2-intel-prod-openmp/pp'
    )
    # AM4a1 +2K
    am4_a1p2k = Run(
        proj=a3gcm,
        name='+2K',
        description='',
        direc_nc='/archive/Ming.Zhao/awg/tikal_201403/c96L48_am4a1_' + \
                 '2000climo_highsen1_p2K/gfdl.ncrc2-intel-prod-openmp/pp'
    )
    # # AM4a2 control
    # am4_a2c = Run(
    #     proj=a3gcm,
    #     name='cont',
    #     description='',
    #     direc_nc='/archive/cjg/awg/tikal_201403/c96L48_am4a2r1_' + \
    #              '2000climo/gfdl.ncrc2-intel-prod-openmp/pp'
    # )
    # # AM4a2 +2K
    # am4_a2p2k = Run(
    #     proj=a3gcm,
    #     name='+2K',
    #     description='',
    #     direc_nc='/archive/cjg/awg/tikal_201403/c96L48_am4a2r1_' + \
    #              '2000climo_p2K/gfdl.ncrc2-intel-prod-openmp/pp'
    # )
    # AM4c1 control
    am4_c1c = Run(
        proj=a3gcm,
        name='cont',
        description='',
        direc_nc='/archive/miz/tikal_201409_awgUpdates_mom6_2014.08.29/' + \
                 'c96L48_am4c1r2_2000climo/gfdl.ncrc2-intel-prod-openmp/pp'
    )
    # AM4c1 +2K
    am4_c1p2k = Run(
        proj=a3gcm,
        name='+2K',
        description='',
        direc_nc='/archive/miz/tikal_201409_awgUpdates_mom6_2014.08.29/' + \
                 'c96L48_am4c1r2_2000climo_p2K/gfdl.ncrc2-intel-prod-openmp/pp'
    )

    ### Models ###
    am2 = Model(
        name='am2',
        proj=a3gcm,
        nc_grid_paths=(
            ('/archive/s1h/am2/am2clim_reyoi/gfdl.ncrc2-default-prod/'
             'pp/atmos/atmos.static.nc'),
            # ['/archive/s1h/am2/am2clim_reyoi/gfdl.ncrc2-default-prod/pp/'
             # 'atmos/ts/monthly/1yr/atmos.' + str(n) + '01-' + str(n) +
             # '12.ucomp.nc' for n in range(1982, 2013)],
            ('/archive/yim/siena_201203/m45_am2p14_1990/gfdl.ncrc2-intel-prod/'
             'pp/atmos/ts/monthly/16yr/atmos.198301-199812.temp.nc'),
            ('/archive/s1h/am2/am2clim_reyoi/gfdl.ncrc2-default-prod/'
             'pp/atmos_level/atmos_level.static.nc'),
            ('/archive/s1h/am2/am2clim_reyoi/gfdl.ncrc2-default-prod/'
             'pp/atmos_level/ts/monthly/1yr/atmos_level.198201-198212.temp.nc'),
            ('/archive/s1h/am2/am2clim_reyoi/gfdl.ncrc2-default-prod/'
             'pp/atmos_level/ts/monthly/1yr/atmos_level.198201-198212.hght.nc')
        ),
        nc_dur=1,
        nc_start_yr=1982,
        nc_end_yr=2012,
        default_yr_range=(1983, 2012),
        runs=[
            # am2_cont, am2_aero, am2_atm, am2_amtm, am2_gas, am2_gtm, am2_gmtm,
            # am2_aatl, am2_aind, am2_apac, am2_noT, am2_noT_p2K,
            am2_amip,
            am2_reyoi_cont, am2_reyoi_m0p25, am2_reyoi_m0p5, am2_reyoi_m1,
            am2_reyoi_m1p5, am2_reyoi_m2, am2_reyoi_m3, am2_reyoi_m4,
            am2_reyoi_p0p25, am2_reyoi_p0p5, am2_reyoi_p1, am2_reyoi_p1p5,
            am2_reyoi_p2, am2_reyoi_p3, am2_reyoi_p4, am2_reyoi_p6,
            am2_reyoi_p8, am2_reyoi_m6, am2_reyoi_m8, am2_reyoi_m10,
            am2_reyoi_p10, am2_reyoi_wpwp_p2, am2_reyoi_wpwp_m2,
            am2_cld_lock_cont, am2_cld_lock_p2, am2_cld_lock_sst,
            am2_cld_lock_cld
            # am2_amip1, am2_amip1_p2, am2_reynolds,
            # am2_reynolds_p2, am2_hurrell_cont, am2_hurrell_p2,
            # am2_cld_seed_all_p2, am2_cld_seed_np_p2
            # am2_cld_seed_sp_p2, am2_cld_seed_sa_p2, 
        ],
        default_runs=[am2_reyoi_cont, am2_reyoi_p2]
    )
    am3 = Model(
        name='am3',
        proj=a3gcm,
        nc_grid_paths=(
            '/archive/s1h/am3/am3clim_hurrell/gfdl.ncrc2-intel-prod-openmp/'
            'pp/atmos/atmos.static.nc',
            ['/archive/s1h/am3/am3clim_hurrell/gfdl.ncrc2-intel-prod-openmp/pp/'
             'atmos/ts/monthly/1yr/atmos.' + str(n) + '01-' + str(n) +
             '12.ucomp.nc' for n in range(1980, 2011)]
        ),
        nc_dur=1,
        nc_start_yr=1980,
        nc_end_yr=2010,
        default_yr_range=(1981,2010),
        runs=[
            # am3_cont, am3_aero, am3_atm, am3_amtm, am3_gas, am3_gtm, am3_gmtm,
            # am3_aatl, am3_aind, am3_apac,
            am3_hc, am3_hp1k, am3_hp2k, am3_hp4k, am3_hp6k, am3_hp8k,
            am3_hp10k, am3_hm1k, am3_hm2k, am3_hm4k, am3_hm6k, am3_hm8k,
            am3_hm10k, am3_hm15k, am3_amip, am3_hwpwp_p2k
        ],
        default_runs=[am3_hc, am3_hp2k]
    )
    hiram = Model(
        name='hiram',
        proj=a3gcm,
        nc_grid_paths=(
            '/archive/yim/siena_201211/c180_hiram_clim/gfdl.ncrc2-default-prod/'
            'pp/atmos/atmos.static.nc',
            '/archive/yim/siena_201211/c180_hiram_clim/gfdl.ncrc2-default-prod/'
            'pp/atmos/ts/monthly/17yr/atmos.197901-199512.ucomp.nc'
        ),
        nc_dur=17,
        nc_start_yr=1979,
        nc_end_yr=1995,
        default_yr_range=(1979,1995),
        runs=[hiram_amip,
            hiram_cont, #hiram_aero, hiram_atm, hiram_amtm, hiram_gas,
            hiram_gtm#, hiram_gmtm, hiram_aatl, hiram_aind, hiram_apac
        ],
        default_runs=[hiram_cont, hiram_gtm]
    )
    # sm2 = Model(
    #     name='sm2',
    #     proj=a3gcm,
    #     description='AM2.1 atmosphere coupled to mixed-layer ocean.',
    #     nc_grid_path='/archive/s1h/nc_grid/a3gcm.am2.nc',
    #     lat_bounds_name='latb',
    #     lon_bounds_name='lonb',
    #     land_mask_path='/archive/s1h/nc_grid/walker.cm3.atmos_av.nc',
    #     direc_nc='/archive/yim/sm2.1_fixed/',
    #     nc_dur=20,
    #     nc_start_yr=1,
    #     default_yr_range=(61, 80),
    #     runs=[sm2_cont, sm2_aero, sm2_gas, sm2_both]
    # )
    hiram_c48 = Model(
        name='hiram_mz',
        proj=a3gcm,
        description=('Low resolution version of HiRAM used by Zhao 2014,'
                     ' J. Climate.'),
        nc_grid_paths=(
            '/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0/'
            'gfdl.ncrc2-intel-prod/pp/atmos/atmos.static.nc',
            '/archive/Ming.Zhao/hiramdp/siena_201204/c48l32_him_X0/'
            'gfdl.ncrc2-intel-prod/pp/atmos/ts/monthly/15yr/'
            'atmos.198101-199512.ucomp.nc'
        ),
        nc_dur=15,
        nc_start_yr=1981,
        nc_end_yr=1995,
        default_yr_range=(1981,1995),
        runs=[
            hiram_c48_0, hiram_c48_0_p2K#, hiram_c48_1, hiram_c48_1_p2K,
            # hiram_c48_2, hiram_c48_2_p2K, hiram_c48_3, hiram_c48_3_p2K,
            # hiram_c48_4, hiram_c48_4_p2K, hiram_c48_5, hiram_c48_5_p2K,
            # hiram_c48_6, hiram_c48_6_p2K, hiram_c48_7, hiram_c48_7_p2K,
            # hiram_c48_8, hiram_c48_8_p2K
        ],
        default_runs=[hiram_c48_0, hiram_c48_0_p2K]
    )
    # am3c90 = Model(
    #     name='am3c90',
    #     proj=a3gcm,
    #     nc_grid_path='/archive/s1h/nc_grid/am3c90.atmos_ts.nc',
    #     land_mask_path='/archive/s1h/nc_grid/am3c90.atmos_static.nc',
    #     nc_dur=10,
    #     nc_start_yr=1981,
    #     nc_end_yr=1990,
    #     default_yr_range=(1981,1990),
    #     runs=[am3c90_cont, am3c90_p2K]
    # )
    # am2p5 = Model(
    #     name='am2p5',
    #     proj=a3gcm,
    #     nc_grid_path='/archive/s1h/nc_grid/am2p5.atmos_ts.nc',
    #     # land_mask_path='/archive/s1h/nc_grid/am2p5.atmos_static.nc',
    #     land_mask_path=None,
    #     nc_dur=10,
    #     nc_start_yr=1981,
    #     nc_end_yr=2000,
    #     default_yr_range=(1981,2000),
    #     runs=[am2p5_cont, am2p5_p2K]
    # )
    am4a1 = Model(
        name='am4a1',
        proj=a3gcm,
        nc_grid_paths=(
            '/archive/Ming.Zhao/awg/tikal_201403/c96L48_am4a1_'
            '2000climo_highsen1/gfdl.ncrc2-intel-prod-openmp/pp/atmos/'
            'atmos.static.nc',
            ['/archive/Ming.Zhao/awg/tikal_201403/c96L48_am4a1_'
            '2000climo_highsen1/gfdl.ncrc2-intel-prod-openmp/pp/atmos/'
            'ts/monthly/1yr/atmos.00%02d01-00%02d12.temp.nc' % (n, n)
             for n in range(2,12)]
        ),
        nc_dur=1,
        nc_start_yr=2,
        nc_end_yr=11,
        default_yr_range=(2,11),
        runs=[am4_a1c, am4_a1p2k]
    )
    # am4a2 = Model(
    #     name='am4a2',
    #     proj=a3gcm,
    #     nc_grid_path='/archive/s1h/nc_grid/am4.ts.nc',
    #     land_mask_path='/archive/s1h/nc_grid/am4.static.nc',
    #     nc_dur=5,
    #     nc_start_yr=2,
    #     nc_end_yr=11,
    #     default_yr_range=(2,11),
    #     runs=[am4_a2c, am4_a2p2k]
    # )
    am4c1 = Model(
        name='am4c1',
        proj=a3gcm,
        nc_grid_paths=(
            '/archive/miz/tikal_201409_awgUpdates_mom6_2014.08.29/' 
            'c96L48_am4c1r2_2000climo/gfdl.ncrc2-intel-prod-openmp/pp/'
            'atmos/atmos.static.nc',
            '/archive/miz/tikal_201409_awgUpdates_mom6_2014.08.29/' 
            'c96L48_am4c1r2_2000climo/gfdl.ncrc2-intel-prod-openmp/pp/'
            'atmos/ts/monthly/10yr/atmos.000101-001012.temp.nc'
        ),
        nc_dur=10,
        nc_start_yr=1,
        nc_end_yr=10,
        default_yr_range=(1,10),
        runs=[am4_c1c, am4_c1p2k]
    )
    ### Set Models and Regions of the Proj.
    _set_named_attr_dict(
        a3gcm, 'models', [am2, am3, hiram, hiram_c48, am4a1, am4c1]
                          #sm2, am3c90, am2p5, am4a2
    )
    _set_named_attr_dict(a3gcm, 'default_models', [am2, am3, hiram, hiram_c48])
    a3gcm.regions = [
        globe, nh, sh, tropics, wpwp, epac, sahel, sahel2, sahel3, sahara, 
        ind_monsoon, land, ocean, trop_land, trop_ocean, sahel_south,
        sahel_north, sahel_east, sahel_west
    ]
    return a3gcm
