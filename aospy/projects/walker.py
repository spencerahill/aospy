from aospy import Proj, Model, Run
from aospy.av_stat import grid_sfc_area
from aospy.regions import globe, nh, sh, tropics, wpwp, epac, sahel
from aospy.regions import burls_wpac, burls_epac, burls_trop_pac, burls_ext_pac
### Runs ###
## AM3 runs. 'wa' prefix is for 'Walker AM3'
# wa_pi
wa_pi = Run()
wa_pi.name = 'preind'
wa_pi.description = 'Observed historical SSTs with preindustrial forcing.'
wa_pi.ens_mem_prefix = '/archive/lwh/fms/riga_201104/c48L48_am3p9_1860_'
wa_pi.ens_mem_ext = ['ext', 'ext2', 'ext3']
wa_pi.ens_mem_suffix = '/gfdl.intel-prod/pp'
wa_pi.set_direc()
wa_pi.nc_dur = 5
wa_pi.nc_start_yr = 1870
wa_pi.nc_end_yr = 2004
# wa_allf
wa_allf = Run()
wa_allf.name = 'all_forc'
wa_allf.description = 'Observed historical SSTs with all historical forcings, anthropogenic and natural.'
wa_allf.ens_mem_prefix = '/archive/lwh/fms/riga_201104/c48L48_am3p9_'
wa_allf.ens_mem_ext = ['ext', 'ext2', 'ext3']
wa_allf.ens_mem_suffix = '/gfdl.intel-prod/pp'
wa_allf.set_direc()
wa_allf.nc_dur = 5
wa_allf.nc_start_yr = 1870
wa_allf.nc_end_yr = 2004
# wa_aero
wa_aero = Run()
wa_aero.name = 'aero'
wa_aero.description = 'Observed historical SSTs with historical aerosol emissions and otherwise PI atmospheric composition.'
wa_aero.ens_mem_prefix = '/archive/lwh/fms/riga_201104/c48L48_am3p9_aeroOnly_'
wa_aero.ens_mem_ext = ['ext', 'ext2', 'ext3']
wa_aero.ens_mem_suffix = '/gfdl.intel-prod/pp'
wa_aero.set_direc()
wa_aero.nc_dur = 5
wa_aero.nc_start_yr = 1870
wa_aero.nc_end_yr = 2004
# wa_bc
wa_bc = Run()
wa_bc.name = 'aero'
wa_bc.description = 'Observed historical SSTs with historical black carbon aerosol emissions and otherwise PI atmospheric composition.'
wa_bc.ens_mem_prefix = '/archive/lwh/fms/riga_201104/c48L48_am3p9_BCOnly_'
wa_bc.ens_mem_ext = ['ext', 'ext2', 'ext3']
wa_bc.ens_mem_suffix = '/gfdl.intel-prod/pp'
wa_bc.set_direc()
wa_bc.nc_dur = 5
wa_bc.nc_start_yr = 1870
wa_bc.nc_end_yr = 2004
# wa_sulf
wa_sulf = Run()
wa_sulf.name = 'aero'
wa_sulf.description = 'Observed historical SSTs with historical sulfate aerosol emissions and otherwise PI atmospheric composition.'
wa_sulf.ens_mem_prefix = '/archive/lwh/fms/riga_201104/c48L48_am3p9_sulfateOnly_'
wa_sulf.ens_mem_ext = ['ext', 'ext2', 'ext3']
wa_sulf.ens_mem_suffix = '/gfdl.intel-prod/pp'
wa_sulf.set_direc()
wa_sulf.nc_dur = 5
wa_sulf.nc_start_yr = 1870
wa_sulf.nc_end_yr = 2004
# wa_volc
wa_volc = Run()
wa_volc.name = 'aero'
wa_volc.description = 'Observed historical SSTs with historical volcano emissions and otherwise PI atmospheric composition.'
wa_volc.ens_mem_prefix = '/archive/lwh/fms/riga_201104/c48L48_am3p9_volcOnly_'
wa_volc.ens_mem_ext = ['ext', 'ext2', 'ext3']
wa_volc.ens_mem_suffix = '/gfdl.intel-prod/pp'
wa_volc.set_direc()
wa_volc.nc_dur = 5
wa_volc.nc_start_yr = 1870
wa_volc.nc_end_yr = 2004
# wa_wmggO3
wa_wmggO3 = Run()
wa_wmggO3.name = 'wmggO3'
wa_wmggO3.description = 'Observed historical SSTs with historical well-mixed greenhouse gas and ozone emissions and otherwise PI atmospheric composition.'
wa_wmggO3.ens_mem_prefix = '/archive/lwh/fms/riga_201104/c48L48_am3p9_WMGGO3Only_'
wa_wmggO3.ens_mem_ext = ['ext', 'ext2', 'ext3']
wa_wmggO3.ens_mem_suffix = '/gfdl.intel-prod/pp'
wa_wmggO3.set_direc()
wa_wmggO3.nc_dur = 5
wa_wmggO3.nc_start_yr = 1870
wa_wmggO3.nc_end_yr = 2004
# AM3 climatological SSTs, PI forcing
wa_clim = Run()
wa_clim.name = 'climSST'
wa_clim.description = 'Prescribed observed climatological SSTs with PI atmospheric composition.'
wa_clim.direc_nc = '/archive/lwh/fms/preR/c48L48_am3p9/pp'
wa_clim.nc_dur = 20
wa_clim.nc_start_yr = 1981
wa_clim.nc_end_yr = 2000
# AM3 climatological SSTs, PD aerosol forcing
wa_clim_aero = Run()
wa_clim_aero.name = 'climSSTaero'
wa_clim_aero.description = 'Prescribed historical climatological SSTs with PD aerosols and otherwise PI atmospheric composition'
wa_clim_aero.direc_nc = '/archive/lwh/fms/preR/c48L48_am3p9_1860aero/pp'
wa_clim_aero.nc_dur = 20
wa_clim_aero.nc_start_yr = 1981
wa_clim_aero.nc_end_yr = 2000

## CM3 runs
# CM3 aero
wc_aero = Run()
wc_aero.name = 'aero'
wc_aero.description = 'Preindustrial to present-day aerosol emissions; otherwise preindustrial atmospheric composition.'
wc_aero.ens_mem_prefix = '/archive/cm3/ipcc_ar5/CM3Z_D1_1860-2005_Aerosol_'
wc_aero.ens_mem_ext = ['P1', 'P3', 'P5']
wc_aero.ens_mem_suffix = '/pp'
wc_aero.set_direc()
wc_aero.nc_dur = 146
wc_aero.nc_start_yr = 1860
wc_aero.nc_end_yr = 2005
# CM3 AIE
wc_aie = Run()
wc_aie.name = 'aero_ind'
wc_aie.description = ''
wc_aie.ens_mem_prefix = '/archive/cm3/ipcc_ar5/CM3Z_D1_1860-2005_AerosolIndirect_'
wc_aie.ens_mem_ext = ['PI1', 'PI3', 'PI5']
wc_aie.ens_mem_suffix = '/pp'
wc_aie.set_direc()
wc_aie.nc_dur = 146
wc_aie.nc_start_yr = 1860
wc_aie.nc_end_yr = 2005
# CM3 all forcing
wc_allf = Run()
wc_allf.name = 'all_forc'
wc_allf.description = 'Historical all forcings.'
wc_allf.ens_mem_prefix = '/archive/cm3/ipcc_ar5/CM3Z_D1_1860-2005_AllForc_'
wc_allf.ens_mem_ext = ['H1', 'H2', 'H3', 'H4', 'H5']
wc_allf.ens_mem_suffix = '/pp'
wc_allf.set_direc()
wc_allf.nc_dur = 146
wc_allf.nc_start_yr = 1860
wc_allf.nc_end_yr = 2005
# CM3 anthro
wc_anth = Run()
wc_anth.name = 'anthro'
wc_anth.description = 'All historical anthropogenic emissions.'
wc_anth.ens_mem_prefix = '/archive/cm3/ipcc_ar5/CM3Z_D1_1860-2005_Anthro_'
wc_anth.ens_mem_ext = ['A1', 'A3', 'A5']
wc_anth.ens_mem_suffix = '/pp'
wc_anth.set_direc()
wc_anth.nc_dur = 146
wc_anth.nc_start_yr = 1860
wc_anth.nc_end_yr = 2005
# CM3 black carbon
wc_bc = Run()
wc_bc.name = 'bc'
wc_bc.description = 'Black carbon only.'
wc_bc.ens_mem_prefix = '/archive/lwh/fms/siena/CM3Z_D1_1860-2005_BCOnly_'
wc_bc.ens_mem_ext = ['P1BC', 'P3BC', 'P5BC']
wc_bc.ens_mem_suffix = '/gfdl.ncrc2-intel-prod-openmp/pp'
wc_bc.set_direc()
wc_bc.nc_dur = 146
wc_bc.nc_start_yr = 1860
wc_bc.nc_end_yr = 2005
# CM3 cont
wc_cont = Run()
wc_cont.name = 'cont'
wc_cont.description = 'Preindustrial control.'
wc_cont.direc_nc = '/archive/cm3/ipcc_ar5/CM3Z_Control-1860_D1/pp'
wc_cont.nc_dur = 100
wc_cont.nc_start_yr = 1
wc_cont.nc_end_yr = 3800
# CM3 natural
wc_nat = Run()
wc_nat.name = 'natural'
wc_nat.description = 'Natural forcings only; includes volcanoes.'
wc_nat.ens_mem_prefix = '/archive/cm3/ipcc_ar5/CM3Z_D1_1860-2005_Natural_'
wc_nat.ens_mem_ext = ['N1', 'N3', 'N5']
wc_nat.ens_mem_suffix = '/pp'
wc_nat.set_direc()
wc_nat.nc_dur = 146
wc_nat.nc_start_yr = 1860
wc_nat.nc_end_yr = 2005
# CM3 sulfate
wc_sulf = Run()
wc_sulf.name = 'sulfate'
wc_sulf.description = 'Sulfate only.'
wc_sulf.ens_mem_prefix = '/archive/lwh/fms/siena/CM3Z_D1_1860-2005_sulfateOnly_'
wc_sulf.ens_mem_ext = ['P1SO4', 'P3SO4', 'P5SO4']
wc_sulf.ens_mem_suffix = '/gfdl.ncrc2-intel-prod-openmp/pp'
wc_sulf.set_direc()
wc_sulf.nc_dur = 146
wc_sulf.nc_start_yr = 1860
wc_sulf.nc_end_yr = 2005
# CM3 sulfate and OC
wc_sulf_oc = Run()
wc_sulf_oc.name = 'sulf_oc'
wc_sulf_oc.description = 'Sulfate and organic carbon only.'
wc_sulf_oc.ens_mem_prefix = '/archive/lwh/fms/siena/CM3Z_D1_1860-2005_sulfateOCOnly_'
wc_sulf_oc.ens_mem_ext = ['P1SO4OC', 'P3SO4OC', 'P5SO4OC']
wc_sulf_oc.ens_mem_suffix = '/gfdl.ncrc2-intel-prod-openmp/pp'
wc_sulf_oc.set_direc()
wc_sulf.nc_dur = 146
wc_sulf.nc_start_yr = 1860
wc_sulf.nc_end_yr = 2005
# CM3 WMGG and O3
wc_wmggO3 = Run()
wc_wmggO3.name = 'wmggO3'
wc_wmggO3.description = 'Preindustrial to present-day well-mixed greenhouse gas and ozone emissions; otherwise preindustrial atmospheric composition.'
wc_wmggO3.ens_mem_prefix = '/archive/cm3/ipcc_ar5/CM3Z_D1_1860-2005_WMGGO3_'
wc_wmggO3.ens_mem_ext = ['G1', 'G3', 'G5']
wc_wmggO3.ens_mem_suffix = '/pp'
wc_wmggO3.set_direc()
wc_wmggO3.nc_dur = 146
wc_wmggO3.nc_start_yr = 1860
wc_wmggO3.nc_end_yr = 2005
# CM3 Chinese aerosols
wc_china = Run()
wc_china.name = 'china'
wc_china.description = 'Historical aerosol emissions over China only; otherwise preindustrial.'
wc_china.ens_mem_prefix = '/archive/ds/fms/siena/CM3Z_D1_1860-2005_AllForc_1860AerChn_'
wc_china.ens_mem_ext = ['Ch1', 'Ch3', 'Ch5']
wc_china.ens_mem_suffix = '/gfdl.ncrc2-default-prod/pp'
wc_china.set_direc()
wc_china.nc_dur = 146
wc_china.nc_start_yr = 1860
wc_china.nc_end_yr = 2005
# CM3 historical aerosols until 1980 then flat
wc_aero1980 = Run()
wc_aero1980.name = 'aero1980'
wc_aero1980.description = 'Historical aerosol emissions up to 1980, then held fixed at 1980 levels.'
wc_aero1980.ens_mem_prefix = '/archive/lwh/fms/siena/CM3Z_D1_1980-2005_1980Aerosol_'
wc_aero1980.ens_mem_ext = ['P1', 'P3', 'P5']
wc_aero1980.ens_mem_suffix = '/gfdl.ncrc2-intel-prod-openmp/pp'
wc_aero1980.set_direc()
wc_aero1980.nc_dur = 26
wc_aero1980.nc_start_yr = 1980
wc_aero1980.nc_end_yr = 2005

### Models ###
# am3
wam3 = Model()
wam3.name = 'am3'
wam3.description = 'AM3 with prescribed observed SSTs 1860-2005'
wam3.set_nc_grid('/home/s1h/nc_grid/walker.cm3.atmos_ts.nc')
wam3.lat_name = 'lat'
wam3.lat = wam3.nc_grid.variables[wam3.lat_name][:]
wam3.lat_bounds_name = 'lat_bnds'
wam3.lat_bounds = wam3.nc_grid.variables[wam3.lat_bounds_name][:]
wam3.lon_name = 'lon'
wam3.lon = wam3.nc_grid.variables[wam3.lon_name][:]
wam3.lon_bounds_name = 'lon_bnds'
wam3.lon_bounds = wam3.nc_grid.variables[wam3.lon_bounds_name][:]
wam3.sfc_area = grid_sfc_area(wam3.nc_grid)
wam3.level_name = 'level'
wam3.level = wam3.nc_grid.variables[wam3.level_name][:]
wam3.set_land_mask('/home/s1h/nc_grid/walker.cm3.atmos_av.nc')
wam3.dt_name = 'average_DT'
wam3.direc_nc = '/archive/lwh/fms/riga_201104/'
wam3.nc_dur = 136
wam3.nc_start_yr = 1870
wam3.default_yr_range = (1870, 2005)
wam3.set_runs([wa_pi, wa_aero, wa_wmggO3, wa_allf, wa_clim, wa_clim_aero])
# cm3
wcm3 = Model()
wcm3.name = 'cm3'
wcm3.description = 'Fully coupled CM3 AOGCM'
wcm3.set_nc_grid('/home/s1h/nc_grid/walker.cm3.atmos_ts.nc')
wcm3.lat_name = 'lat'
wcm3.lat = wcm3.nc_grid.variables[wcm3.lat_name][:]
wcm3.lat_bounds_name = 'lat_bnds'
wcm3.lat_bounds = wcm3.nc_grid.variables[wcm3.lat_bounds_name][:]
wcm3.lon_name = 'lon'
wcm3.lon = wcm3.nc_grid.variables[wcm3.lon_name][:]
wcm3.lon_bounds_name = 'lon_bnds'
wcm3.lon_bounds = wcm3.nc_grid.variables[wcm3.lon_bounds_name][:]
wcm3.sfc_area = grid_sfc_area(wcm3.nc_grid)
wcm3.level_name = 'level'
wcm3.level = wcm3.nc_grid.variables[wcm3.level_name][:]
wcm3.set_land_mask('/home/s1h/nc_grid/walker.cm3.atmos_av.nc')
wcm3.dt_name = 'average_DT'
wcm3.direc_nc = '/archive/cm3/ipcc_ar5/'
wcm3.nc_dur = 146
wcm3.nc_start_yr = 1860
wcm3.default_yr_range = (1860, 2005)
wcm3.set_runs([wc_aero, wc_aero1980, wc_aie, wc_allf, wc_anth, wc_cont, wc_nat, 
               wc_sulf, wc_wmggO3, wc_china])

### Proj ###
walker = Proj()
walker.name = 'walker'
walker.direc_nc = '/archive/cm3/ipcc_ar5/'
walker.direc_out = '/archive/s1h/walker/'
walker.set_models ([wam3, wcm3])
walker.regions = [globe, nh, sh, tropics,  wpwp, epac, 
                  burls_wpac, burls_epac, burls_trop_pac, burls_ext_pac]
