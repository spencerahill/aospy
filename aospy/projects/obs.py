def obs():
    from aospy import Proj, Model, Run, _set_named_attr_dict
    from aospy.regions import *
    from aospy.variables import master_vars_list

    ### Proj ###
    obs_obj = Proj(
        name='obs', 
        vars=master_vars_list,
        direc_out='/archive/s1h/obs/',
    )
    ### Runs ###
    ## CRU
    # v3.22
    cru_v322 = Run(
        proj=obs_obj,
        name='v3.22',
        description='',
        direc_nc='/archive/s1h/obs/HadCRU/3.22',
        nc_dir_structure='one_dir',
        nc_dur=113,
        nc_start_yr=1901,
        nc_start_month=1,
        nc_end_yr=2013,
        nc_end_month=12,
        default_yr_range=(1901, 2013),
        nc_files={'precip': 'cru_ts3.22.1901.2013.pre.dat.nc',
                  'cld_amt': 'cru_ts3.22.1901.2013.cld.dat.nc',
                  'diurnal_temp_range': 'cru_ts3.22.1901.2013.dtr.dat.nc',
                  'ground_frost_freq': 'cru_ts3.22.1901.2013.frs.dat.nc',
                  'pet': 'cru_ts3.22.1901.2013.pet.dat.nc',
                  't_surf_min': 'cru_ts3.22.1901.2013.tmn.dat.nc',
                  't_surf_max': 'cru_ts3.22.1901.2013.tmx.dat.nc',
                  't_surf': 'cru_ts3.22.1901.2013.tmp.dat.nc',
                  'vap_pres': 'cru_ts3.22.1901.2013.vap.dat.nc',
                  'wet_day_freq': 'cru_ts3.22.1901.2013.wet.dat.nc'}
    )
    ## PREC/L
    # 0.5x0.5 degree resolution
    prec_l_0p5deg = Run(
        proj=obs_obj,
        name='0.5deg',
        description='',
        direc_nc='/archive/s1h/obs/PREC_L/20150212',
        nc_dir_structure='one_dir',
        nc_dur=64,
        nc_start_yr=1948,
        nc_start_month=1,
        nc_end_yr=2011,
        nc_end_month=12,
        default_yr_range=(1948, 2011),
        nc_files={'precip': 'precip.mon.mean.0.5x0.5.nc'}
    )
    # 1x1 degree resolution
    prec_l_1deg = Run(
        proj=obs_obj,
        name='1deg',
        description='',
        direc_nc='/archive/s1h/obs/PREC_L/20150212',
        nc_dir_structure='one_dir',
        nc_dur=67,
        nc_start_yr=1948,
        nc_start_month=1,
        nc_end_yr=2014,
        nc_end_month=12,
        default_yr_range=(1948, 2014),
        nc_files={'precip': 'precip.mon.mean.1x1.nc'}
    )
    # 2.5x2.5 degree resolution
    prec_l_2p5deg = Run(
        proj=obs_obj,
        name='2.5deg',
        description='',
        direc_nc='/archive/s1h/obs/PREC_L/20150212',
        nc_dir_structure='one_dir',
        nc_dur=67,
        nc_start_yr=1948,
        nc_start_month=1,
        nc_end_yr=2014,
        nc_end_month=12,
        default_yr_range=(1948, 2014),
        nc_files={'precip': 'precip.mon.mean.2.5x2.5.nc'}
    )
    ## CERES runs
    # CERES EBAF
    ceres_ebaf = Run(
        proj=obs_obj,
        name='ebaf',
        description='',
        direc_nc=('/archive/pcmdi/repo/obs4MIPs/NASA-LaRC/CERES-EBAF/'
                  'atmos/mon/v20140402'),
        nc_dir_structure='one_dir',
        nc_dur=14,
        nc_start_yr=2000,
        nc_start_month=3,
        nc_end_yr=2013,
        nc_end_month=10,
        default_yr_range=(2000, 2013),
        nc_suffix='_CERES-EBAF_L3B_Ed2-8_200003-201310.nc',
        nc_files={
            'swdn_toa': 'rsdt_CERES-EBAF_L3B_Ed2-8_200003-201310.nc',
            'swup_toa': 'rsut_CERES-EBAF_L3B_Ed2-8_200003-201310.nc',
            'swup_toa_clr': 'rsutcs_CERES-EBAF_L3B_Ed2-8_200003-201310.nc',
            'olr': 'rlut_CERES-EBAF_L3B_Ed2-8_200003-201310.nc',
            'olr_clr': 'rlutcs_CERES-EBAF_L3B_Ed2-8_200003-201310.nc'
        }
    )
    ## GPCP runs
    # GPCP v2p2
    gpcp_v2p2=Run(
        proj=obs_obj, 
        name='v2p2',
        description='GPCP v2.2 gridded precipitation, from blend of ' + \
                    'satellite and station gauge data.',
        direc_nc='/archive/pcmdi/repo/obs4MIPs/NASA-GSFC/GPCP/atmos/',
        nc_dur=10,
        nc_start_yr=1979,
        nc_end_yr=2013,
        default_yr_range=(1979, 2013),
        nc_dir_structure='one_dir',
        nc_files={'monthly':
                  ['mon/v20130401/pr_GPCP-SG_L3_v2.2_' + yrs + '.nc' for yrs in
                   ('197901-197912', '198001-198912', '199001-199912',
                    '200001-200912', '201001-201312')],
                  'pentad': 'day/v20121003/'}
    )
    ## TRMM runs
    # TRMM v7
    trmm_v7a = Run(
        proj=obs_obj,
        name='v7a',
        description='TRMM v7 gridded precipitation, from satellite data',
        direc_nc='/archive/pcmdi/repo/obs4MIPs/NASA-GSFC/TRMM/atmos/',
        nc_dur=2,
        nc_start_yr=2000,
        nc_start_month=1,
        nc_end_yr=2010,
        nc_end_month=9,
        default_yr_range=(2000, 2009),
        nc_dir_structure='one_dir',
        nc_files={'monthly': ['mon/v20130204/pr_TRMM-L3_v7A_' + yrs + '.nc'
                              for yrs in ('200001-200912', '201001-201009')]}

    )
    ## CMAP runs
    cmap_standard = Run(
        proj=obs_obj,
        name='standard',
        description=('CMAP standard version, which does not include NCEP '
                     'reanalysis data to fill in gaps.'),
        direc_nc='/archive/s1h/obs/CMAP/standard',
        nc_dir_structure='one_dir',
        nc_dur=36,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2014,
        nc_end_month=12,
        default_yr_range=(1979, 2014),
        nc_files={'monthly': 'precip.mon.mean.nc',
                  'pentad': 'precip.pentad.mean.nc'}
    )
    cmap_enhanced = Run(
        proj=obs_obj, 
        name='enhanced',
        description=('CMAP enhanced version, which includes NCEP reanalysis '
                     'data to fill in gaps.'),
        direc_nc='/archive/s1h/obs/CMAP/enhanced',
        nc_dur=36,
        nc_start_yr=1979,
        nc_end_yr=2014,
        default_yr_range=(1979, 2014),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'precip.mon.mean.nc',
                  'pentad': 'precip.pentad.mean.nc'}
    )
    ## U_Del runs
    # v201
    udel_v201 = Run(
        proj=obs_obj, 
        name='v201',
        description='',
        direc_nc='/archive/s1h/obs/U_Del',
        nc_dur=109,
        nc_start_yr=1900,
        nc_end_yr=2008,
        default_yr_range=(1900, 2008),
        nc_dir_structure='one_dir',
        nc_files={'precip': 'precip.mon.total.v201.nc',
                  't_surf': 'air.mon.total.v201.nc'}
        )
    # v301
    udel_v301 = Run(
        proj=obs_obj, 
        name='v301',
        description='',
        direc_nc='/archive/s1h/obs/U_Del',
        nc_dur=111,
        nc_start_yr=1900,
        nc_end_yr=2010,
        default_yr_range=(1900, 2010),
        nc_dir_structure='one_dir',
        nc_files={'precip': 'precip.mon.total.v301.nc',
                  't_surf': 'air.mon.total.v301.nc'}
        )
    ## ERA-Interim
    era_i = Run(
        proj=obs_obj, 
        name='interim',
        description='',
        direc_nc=('/archive/pcmdi/repo/ana4MIPs/ECMWF/ERA-Interim/atmos/'
                  'mon/v20140416'),
        nc_dur=1,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2013,
        nc_end_month=12,
        default_yr_range=(1979, 2013),
        nc_dir_structure='one_dir',
        # 2015-02-10: HACK: have to manually change variable name in the
        # nc_files dict, which is highly prone to errors and prevents
        # full functionality.  Need to resolve ASAP.
        nc_files={
                  'cld_amt': ['cl_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'evap': ['evspsbl_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'hght': ['zg_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'omega': ['wap_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                            for yrs in [str(yr) + '01-' + str(yr) + '12'
                                        for yr in range(1979, 2014)]],
                  'precip': ['pr_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'ps': ['ps_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'rh': ['hur_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'shflx': ['hfss_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'slp': ['psl_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                            for yrs in [str(yr) + '01-' + str(yr) + '12'
                                        for yr in range(1979, 2014)]],
                  'sphum': ['hus_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'temp': ['ta_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'ucomp': ['ua_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'vcomp': ['va_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'wvp': ['prw_Amon_reanalysis_IFS-Cy31r2_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]]
                  }
        )
    ## MERRA
    merra_run = Run(
        proj=obs_obj, 
        name='merra',
        description='',
        direc_nc=('/archive/pcmdi/repo/ana4MIPs/NASA-GMAO/MERRA/atmos/mon/'
                  'v20140624'),
        nc_dur=1,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2011,
        nc_end_month=12,
        default_yr_range=(1979, 2011),
        nc_dir_structure='one_dir',
        # 2015-02-10: HACK: have to manually change variable name in the
        # nc_files dict, which is highly prone to errors and prevents
        # full functionality.  Need to resolve ASAP.
        nc_files={'cld_amt': ['cl_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'evap': ['evspsbl_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'hght': ['zg_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'omega': ['wap_Amon_reanalysis_MERRA_' + yrs + '.nc'
                            for yrs in [str(yr) + '01-' + str(yr) + '12'
                                        for yr in range(1979, 2012)]],
                  'precip': ['pr_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'ps': ['ps_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'rh': ['hur_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'shflx': ['hfss_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'slp': ['psl_Amon_reanalysis_MERRA_' + yrs + '.nc'
                            for yrs in [str(yr) + '01-' + str(yr) + '12'
                                        for yr in range(1979, 2012)]],
                  'sphum': ['hus_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'temp': ['ta_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'ucomp': ['ua_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'vcomp': ['va_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]],
                  'wvp': ['prw_Amon_reanalysis_MERRA_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2012)]]
                  }
        )
    ## NCEP CFSR
    cfsr_run = Run(
        proj=obs_obj, 
        name='cfsr',
        description='',
        direc_nc=('/archive/pcmdi/repo/ana4MIPs/NOAA-NCEP/CFSR/atmos/'
                  'mon/v20140822'),
        nc_dur=1,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2013,
        nc_end_month=12,
        default_yr_range=(1979, 2013),
        nc_dir_structure='one_dir',
        # 2015-02-10: HACK: have to manually change variable name in the
        # nc_files dict, which is highly prone to errors and prevents
        # full functionality.  Need to resolve ASAP.
        nc_files={'cld_amt': ['cl_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'evap': ['evspsbl_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'hght': ['zg_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'omega': ['wap_Amon_reanalysis_CFSR_' + yrs + '.nc'
                            for yrs in [str(yr) + '01-' + str(yr) + '12'
                                        for yr in range(1979, 2014)]],
                  'precip': ['pr_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'ps': ['ps_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'rh': ['hur_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'shflx': ['hfss_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'slp': ['psl_Amon_reanalysis_CFSR_' + yrs + '.nc'
                            for yrs in [str(yr) + '01-' + str(yr) + '12'
                                        for yr in range(1979, 2014)]],
                  'sphum': ['hus_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'temp': ['ta_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'ucomp': ['ua_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'vcomp': ['va_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]],
                  'wvp': ['prw_Amon_reanalysis_CFSR_' + yrs + '.nc'
                           for yrs in [str(yr) + '01-' + str(yr) + '12'
                                       for yr in range(1979, 2014)]]
                  }

        )
    # ## JMA JRA-25
    # jra25 = Run(
    #     proj=obs_obj, 
    #     name='jra-25',
    #     description='',
    #     direc_nc='/archive/pcmdi/repo/ana4MIPs/JMA/JRA-25/atmos/mon/v20140408',
    #     nc_dur=1,
    #     nc_start_yr=1979,
    #     nc_start_month=1,
    #     nc_end_yr=2013,
    #     nc_end_month=12,
    #     default_yr_range=(1979, 2013),
    #     nc_dir_structure='one_dir',
    #     # 2015-02-10: HACK: have to manually change variable name in the
    #     # nc_files dict, which is highly prone to errors and prevents
    #     # full functionality.  Need to resolve ASAP.
    #     nc_files={'monthly': ['va_Amon_reanalysis_JRA-25_' + yrs + '.nc'
    #                           for yrs in [str(yr) + '01-' + str(yr) + '12'
    #                                       for yr in range(1979, 2014)]]}
    #     )
    ## LandFlux-EVAL 1989-2005 evapotranspiration runs
    # Using all products
    lfe_all = Run(
        proj=obs_obj, 
        name='all',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=2005,
        nc_end_month=12,
        default_yr_range=(1989, 2005),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-05.monthly.all.nc',
                  'annual': 'LandFluxEVAL.merged.89-05.yearly.all.nc'}
        )
    # Using diagnostic-based products only
    lfe_diag = Run(
        proj=obs_obj, 
        name='diagnostic',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=2005,
        nc_end_month=12,
        default_yr_range=(1989, 2005),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-05.monthly.diagnostic.nc',
                  'annual': 'LandFluxEVAL.merged.89-05.yearly.diagnostic.nc'}
        )
    # Using land surface model-based products only
    lfe_lsm = Run(
        proj=obs_obj, 
        name='lsm',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=2005,
        nc_end_month=12,
        default_yr_range=(1989, 2005),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-05.monthly.lsm.nc',
                  'annual': 'LandFluxEVAL.merged.89-05.yearly.lsm.nc'}
        )
    # Using reanalyses-based products only
    lfe_rean = Run(
        proj=obs_obj, 
        name='reanalyses',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=2005,
        nc_end_month=12,
        default_yr_range=(1989, 2005),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-05.monthly.reanalyses.nc',
                  'annual': 'LandFluxEVAL.merged.89-05.yearly.reanlayses.nc'}
        )
    ## LandFlux-EVAL 1989-1995 evapotranspiration runs
    # Using all products
    lfe95_all = Run(
        proj=obs_obj, 
        name='all',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=1995,
        nc_end_month=12,
        default_yr_range=(1989, 1995),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-95.monthly.all.nc',
                  'annual': 'LandFluxEVAL.merged.89-95.yearly.all.nc'}
        )
    # Using diagnostic-based products only
    lfe95_diag = Run(
        proj=obs_obj, 
        name='diagnostic',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=1995,
        nc_end_month=12,
        default_yr_range=(1989, 1995),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-95.monthly.diagnostic.nc',
                  'annual': 'LandFluxEVAL.merged.89-95.yearly.diagnostic.nc'}
        )
    # Using land surface model-based products only
    lfe95_lsm = Run(
        proj=obs_obj, 
        name='lsm',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=1995,
        nc_end_month=12,
        default_yr_range=(1989, 1995),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-95.monthly.lsm.nc',
                  'annual': 'LandFluxEVAL.merged.89-95.yearly.lsm.nc'}
        )
    # Using reanalyses-based products only
    lfe95_rean = Run(
        proj=obs_obj, 
        name='reanalyses',
        description='',
        direc_nc='/archive/s1h/obs/LandFlux-EVAL',
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=1995,
        nc_end_month=12,
        default_yr_range=(1989, 1995),
        nc_dir_structure='one_dir',
        nc_files={'monthly': 'LandFluxEVAL.merged.89-95.monthly.reanalyses.nc',
                  'annual': 'LandFluxEVAL.merged.89-95.yearly.reanlayses.nc'}
        )

    ### Models ###
    # HadCRU
    cru = Model(
        name='cru', 
        description='HadCRU',
        nc_grid_paths=('/archive/s1h/obs/HadCRU/3.22/'
                       'cru_ts3.22.1901.2013.pre.dat.nc',),
        nc_dur=113,
        nc_start_yr=1901,
        nc_end_yr=2013,
        default_yr_range=(1901, 2013),
        runs=[cru_v322],
        default_runs=[cru_v322]        
    )
    # NOAA PREC/L
    prec_l = Model(
        name='prec_l', 
        description='NOAA PRECipitation REConstruction over Land (PREC/L)',
        nc_grid_paths=('/archive/s1h/obs/PREC_L/20150212/'
                       'precip.mon.mean.1x1.nc',),
        nc_dur=64,
        nc_start_yr=1948,
        nc_start_month=1,
        nc_end_yr=2012,
        nc_end_month=12,
        default_yr_range=(1948, 2012),
        # runs=[prec_l_0p5deg, prec_l_1deg, prec_l_2p5deg]
        runs=[prec_l_1deg],
        default_runs=[prec_l_1deg]        
    )
    # GPCP
    gpcp = Model(
        name='gpcp', 
        description='Global Precipitation Climatology Project: ' + \
                    'http://www.gewex.org/gpcp.html',
        nc_grid_paths=([
            '/archive/pcmdi/repo/obs4MIPs/NASA-GSFC/GPCP/atmos/' + \
            'mon/v20130401/pr_GPCP-SG_L3_v2.2_' + yrs + '.nc' for yrs in
            ('197901-197912', '198001-198912', '199001-199912',
             '200001-200912', '201001-201312')
        ],),
        nc_dur=10,
        nc_start_yr=1979,
        nc_end_yr=2013,
        default_yr_range=(1979, 2013),
        runs=[gpcp_v2p2],
        default_runs=[gpcp_v2p2]
    )
    # CERES
    ceres = Model(
        name='ceres', 
        proj=obs_obj,
        nc_grid_paths=('/archive/pcmdi/repo/obs4MIPs/NASA-LaRC/'
                       'CERES-EBAF/atmos/mon/v20140402/'
                       'rsut_CERES-EBAF_L3B_Ed2-8_200003-201310.nc',),
        nc_dur=14,
        nc_start_yr=2000,
        nc_start_month=3,
        nc_end_yr=2013,
        default_yr_range=(2000,2013),
        runs=[ceres_ebaf],
        default_runs=[ceres_ebaf]        
    )
    # CMAP
    cmap = Model(
        name='cmap', 
        proj=obs_obj,
        description='CPC Merged Analysis of Precipitation',
        nc_grid_paths=('/archive/s1h/obs/CMAP/standard/'
                       'precip.mon.mean.nc',),
        nc_dur=36,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2014,
        nc_end_month=12,
        default_yr_range=(1979, 2014),
        runs=[cmap_standard, cmap_enhanced],
        default_runs=[cmap_standard]        
    )
    # TRMM
    trmm = Model(
        name='trmm', 
        description=('Tropical Rainfall Measuring Mission: '
                     'http://trmm.gsfc.nasa.gov/'),
        nc_grid_paths=(['/archive/pcmdi/repo/obs4MIPs/NASA-GSFC/TRMM/atmos/'
                        'mon/v20130204/pr_TRMM-L3_v7A_' + yrs + '.nc' for yrs
                        in ('200001-200912', '201001-201009')],),
        nc_dur=10,
        nc_start_yr=2000,
        nc_start_month=1,
        nc_end_yr=2010,
        nc_end_month=9,
        runs=[trmm_v7a],
        default_runs=[trmm_v7a]
    )
    # U_Del
    udel = Model(
        name='udel', 
        description='U. Delaware gridded land data from station obs',
        nc_grid_paths=('/archive/s1h/obs/U_Del/precip.mon.total.v301.nc',),
        nc_dur=111,
        nc_start_yr=1900,
        runs=[udel_v201, udel_v301],
        default_runs=[udel_v301]        
    )
    # ERA reanalyses
    era = Model(
        name='era', 
        proj=obs_obj,
        nc_grid_paths=(
            ['/archive/pcmdi/repo/ana4MIPs/ECMWF/ERA-Interim/atmos/'
             'mon/v20140416/wap_Amon_reanalysis_IFS-Cy31r2_' + yrs +
             '.nc' for yrs in [str(yr) + '01-' + str(yr) + '12'
                               for yr in range(1979, 2014)]],
        ),
        nc_dur=1,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2013,
        nc_end_month=12,
        default_yr_range=(1979, 2013),
        runs=[era_i],
        default_runs=[era_i]
    )
    # MERRA reanalyses
    merra = Model(
        name='merra', 
        proj=obs_obj,
        nc_grid_paths=(
            # ['/archive/pcmdi/repo/ana4MIPs/NASA-GMAO/MERRA/atmos/'
             # 'mon/v20140624/hfss_Amon_reanalysis_MERRA_' + yrs +
             # '.nc' for yrs in [str(yr) + '01-' + str(yr) + '12'
                               # for yr in range(1979, 2012)]],
            ['/archive/pcmdi/repo/ana4MIPs/NASA-GMAO/MERRA/atmos/'
             'mon/v20140624/wap_Amon_reanalysis_MERRA_' + yrs +
             '.nc' for yrs in [str(yr) + '01-' + str(yr) + '12'
                               for yr in range(1979, 2012)]],
        ),
        nc_dur=1,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2011,
        nc_end_month=12,
        default_yr_range=(1979, 2011),
        runs=[merra_run],
        default_runs=[merra_run]
    )
    # NCEP CFSR reanalyses
    cfsr = Model(
        name='cfsr', 
        proj=obs_obj,
        nc_grid_paths=(
            ['/archive/pcmdi/repo/ana4MIPs/NOAA-NCEP/CFSR/atmos/'
             'mon/v20140822/zg_Amon_reanalysis_CFSR_' + yrs +
             '.nc' for yrs in [str(yr) + '01-' + str(yr) + '12'
                               for yr in range(1979, 2014)]],
        ),
        nc_dur=1,
        nc_start_yr=1979,
        nc_start_month=1,
        nc_end_yr=2013,
        default_yr_range=(1979, 2013),
        runs=[cfsr_run],
        default_runs=[cfsr_run]
    )
    # # JRA-25 reanalyses
    # jra = Model(
    #     name='jra', 
    #     proj=obs_obj,
    #     nc_grid_paths=('/archive/pcmdi/repo/ana4MIPs/JMA/JRA-25/atmos/mon/'
    #                    'v20140408/va_Amon_reanalysis_JRA-25_197901-201312.nc',),
    #     nc_dur=1,
    #     nc_start_yr=1979,
    #     nc_start_month=1,
    #     nc_end_yr=2013,
    #     default_yr_range=(1979, 2013),
    #     runs=[jra25]
    # )
    # LandFlux-EVAL evapotranspiration data, 1989-2005
    landflux = Model(
        name='landflux-eval', 
        proj=obs_obj,
        nc_grid_paths=('/archive/s1h/obs/LandFlux-EVAL/'
                       'LandFluxEVAL.merged.89-05.monthly.all.nc',),
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=2005,
        nc_end_month=12,
        default_yr_range=(1989, 2005),
        runs=[lfe_all, lfe_diag, lfe_lsm, lfe_rean],
        default_runs=[lfe_all]
    )
    # LandFlux-EVAL evapotranspiration data, 1989-1995
    landflux95 = Model(
        name='landflux-eval95', 
        proj=obs_obj,
        nc_grid_paths=('/archive/s1h/obs/LandFlux-EVAL/'
                       'LandFluxEVAL.merged.89-95.monthly.all.nc',),
        nc_dur=17,
        nc_start_yr=1989,
        nc_start_month=1,
        nc_end_yr=1995,
        nc_end_month=12,
        default_yr_range=(1989, 1995),
        runs=[lfe95_all, lfe95_diag, lfe95_lsm, lfe95_rean],
        default_runs=[lfe95_all]
    )

    ### Set Models and Regions of the Proj.
    _set_named_attr_dict(obs_obj, 'models', 
                        [cru, prec_l, gpcp, ceres, trmm, cmap, udel, era, 
                         merra, cfsr, landflux, landflux95]
    )
    _set_named_attr_dict(obs_obj, 'default_models', 
                        [cru, prec_l, gpcp, trmm, cmap, udel]
    )
    obs_obj.regions = [
        globe, nh, sh, tropics, wpwp, epac, sahel, sahel2, sahel3, sahara, 
        ind_monsoon, land, ocean, trop_land, trop_ocean, sahel_south,
        sahel_north, sahel_east, sahel_west
    ]
    return obs_obj
