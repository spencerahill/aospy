from aospy import Proj, Model, Run
from aospy.av_stat import grid_sfc_area
from aospy.regions import globe, nh, sh, tropics, wpwp, epac, sahel
from aospy.regions import burls_wpac, burls_epac, burls_trop_pac, burls_ext_pac
### Runs ###
# 1pctCO2
c5_1pctCO2 = Run()
c5_1pctCO2.name = ''
c5_1pctCO2.description = ''
c5_1pctCO2.ens_mem_prefix = ''
c5_1pctCO2.ens_mem_ext = []
c5_1pctCO2.ens_mem_suffix = ''
c5_1pctCO2.set_direc()
c5_1pctCO2.nc_dur = 0
c5_1pctCO2.nc_start_yr = 0
c5_1pctCO2.nc_end_yr = 0
# 1pctCO2
c5_1pctCO2 = Run()
c5_1pctCO2.name = ''
c5_1pctCO2.description = ''
c5_1pctCO2.ens_mem_prefix = ''
c5_1pctCO2.ens_mem_ext = []
c5_1pctCO2.ens_mem_suffix = ''
c5_1pctCO2.set_direc()
c5_1pctCO2.nc_dur = 0
c5_1pctCO2.nc_start_yr = 0
c5_1pctCO2.nc_end_yr = 0
# 1pctCO2
c5_1pctCO2 = Run()
c5_1pctCO2.name = ''
c5_1pctCO2.description = ''
c5_1pctCO2.ens_mem_prefix = ''
c5_1pctCO2.ens_mem_ext = []
c5_1pctCO2.ens_mem_suffix = ''
c5_1pctCO2.set_direc()
c5_1pctCO2.nc_dur = 0
c5_1pctCO2.nc_start_yr = 0
c5_1pctCO2.nc_end_yr = 0
# 1pctCO2
c5_1pctCO2 = Run()
c5_1pctCO2.name = ''
c5_1pctCO2.description = ''
c5_1pctCO2.ens_mem_prefix = ''
c5_1pctCO2.ens_mem_ext = []
c5_1pctCO2.ens_mem_suffix = ''
c5_1pctCO2.set_direc()
c5_1pctCO2.nc_dur = 0
c5_1pctCO2.nc_start_yr = 0
c5_1pctCO2.nc_end_yr = 0
# 1pctCO2
c5_1pctCO2 = Run()
c5_1pctCO2.name = ''
c5_1pctCO2.description = ''
c5_1pctCO2.ens_mem_prefix = ''
c5_1pctCO2.ens_mem_ext = []
c5_1pctCO2.ens_mem_suffix = ''
c5_1pctCO2.set_direc()
c5_1pctCO2.nc_dur = 0
c5_1pctCO2.nc_start_yr = 0
c5_1pctCO2.nc_end_yr = 0
# 1pctCO2
c5_1pctCO2 = Run()
c5_1pctCO2.name = ''
c5_1pctCO2.description = ''
c5_1pctCO2.ens_mem_prefix = ''
c5_1pctCO2.ens_mem_ext = []
c5_1pctCO2.ens_mem_suffix = ''
c5_1pctCO2.set_direc()
c5_1pctCO2.nc_dur = 0
c5_1pctCO2.nc_start_yr = 0
c5_1pctCO2.nc_end_yr = 0

### Models ###
# GFDL CM3
c5_cm3 = Model()
c5_cm3.name = 'cm3'
c5_cm3.description = 'NOAA-GFDL CM3'
c5_cm3.set_nc_grid('/home/s1h/nc_grid/cmip5.cm3.atmos_ts.nc')
c5_cm3.lat_name = 'lat'
c5_cm3.lat = c5_cm3.nc_grid.variables[c5_cm3.lat_name][:]
c5_cm3.lat_bounds_name = 'lat_bnds'
c5_cm3.lat_bounds = c5_cm3.nc_grid.variables[c5_cm3.lat_bounds_name][:]
c5_cm3.lon_name = 'lon'
c5_cm3.lon = c5_cm3.nc_grid.variables[c5_cm3.lon_name][:]
c5_cm3.lon_bounds_name = 'lon_bnds'
c5_cm3.lon_bounds = c5_cm3.nc_grid.variables[c5_cm3.lon_bounds_name][:]
c5_cm3.sfc_area = grid_sfc_area(c5_cm3.nc_grid)
c5_cm3.level_name = 'level'
c5_cm3.level = c5_cm3.nc_grid.variables[c5_cm3.level_name][:]
c5_cm3.set_land_mask('/home/s1h/nc_grid/cmip5.cm3.atmos_av.nc')
c5_cm3.dt_name = 'average_DT'
c5_cm3.direc_nc = '/archive/cm3/ipcc_ar5/'
c5_cm3.nc_dur = 146
c5_cm3.nc_start_yr = 1860
c5_cm3.runs = {run.name: run for run in 
            [wc_aero, wc_aero1980, wc_aie, wc_allf, wc_anth, wc_cont, wc_nat, 
             wc_sulf, wc_wmggO3, wc_china]}

### Proj ###
cmip5 = Proj()
cmip5.name = 'cmip5'
cmip5.direc_nc = '/archive/pcmdi/repo/CMIP5/'
cmip5.direc_out = '/archive/s1h/cmip5/'
cmip5.models = {model.name: model for model in [wam3, wcm3]}
cmip5.regions = [globe, nh, sh, tropics,  wpwp, epac, 
                 burls_wpac, burls_epac, burls_trop_pac, burls_ext_pac]
