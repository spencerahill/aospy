from aospy import Proj, Model, Run
from aospy.av_stat import grid_sfc_area

#execfile('/home/s1h/py/vars.py')
#execfile('/home/s1h/py/calcs.py')
#execfile('/home/s1h/py/regions.py')

### Vars ###
# Set Var attributes to values specific to this project, as necessary.
# Ideally this will be handled more intelligently eventually.

### Runs ###
## AM2 runs. 'bl' prefix is for 'bc_layers' project; 'a' is for 'AM2.1'
# AM2.1 cont
bla_cont = Run()
bla_cont.name = 'cont'
bla_cont.description = ''
bla_cont.direc_nc = '/archive/yim/fms/lima_mix/PI2/pp'
bla_cont.nc_dur = 17
bla_cont.nc_start_yr = 1983
bla_cont.nc_end_yr = 1998

# AM2.1 layer 11
bla11 = Run()
bla11.name = 'bc_layer11'
bla11.description = ''
bla11.direc_nc = '/archive/yim/fms/lima_mix/PD2_bc_layer11.YIM/pp'
bla11.nc_dur = 17
bla11.nc_start_yr = 1983
bla11.nc_end_yr = 1998

# AM2.1 layer 14
bla14 = Run()
bla14.name = 'bc_layer14'
bla14.description = ''
bla14.direc_nc = '/archive/yim/fms/lima_mix/PD2_bc_layer14.YIM/pp'
bla14.nc_dur = 17
bla14.nc_start_yr = 1983
bla14.nc_end_yr = 1998

# AM2.1 layer 16
bla16 = Run()
bla16.name = 'bc_layer16'
bla16.description = ''
bla16.direc_nc = '/archive/yim/fms/lima_mix/PD2_bc_layer16.YIM/pp'
bla16.nc_dur = 17
bla16.nc_start_yr = 1983
bla16.nc_end_yr = 1998

# AM2.1 layer 18
bla18 = Run()
bla18.name = 'bc_layer18'
bla18.description = ''
bla18.direc_nc = '/archive/yim/fms/lima_mix/PD2_bc_layer18.YIM/pp'
bla18.nc_dur = 17
bla18.nc_start_yr = 1983
bla18.nc_end_yr = 1998

# AM2.1 layer 20
bla20 = Run()
bla20.name = 'bc_layer20'
bla20.description = ''
bla20.direc_nc = '/archive/yim/fms/lima_mix/PD2_bc_layer20.YIM/pp'
bla20.nc_dur = 17
bla20.nc_start_yr = 1983
bla20.nc_end_yr = 1998

# AM2.1 layer 22
bla22 = Run()
bla22.name = 'bc_layer22'
bla22.description = ''
bla22.direc_nc = '/archive/yim/fms/lima_mix/PD2_bc_layer22.YIM/pp'
bla22.nc_dur = 17
bla22.nc_start_yr = 1983
bla22.nc_end_yr = 1998

# AM2.1 layer 24
bla24 = Run()
bla24.name = 'bc_layer24'
bla24.description = ''
bla24.direc_nc = '/archive/yim/fms/lima_mix/PD2_bc_layer24.YIM/pp'
bla24.nc_dur = 17
bla24.nc_start_yr = 1983
bla24.nc_end_yr = 1998

## SM2.1 runs. 's' after 'bl' is for 'SM2.1'
# SM2.1 cont
bls_cont = Run()
bls_cont.name = 'cont'
bls_cont.description = ''
bls_cont.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie_rerun6.YIM/pp'
bls_cont.nc_dur = 20
bls_cont.nc_start_yr = 1
bls_cont.nc_end_yr = 120

# SM2.1 layer 11
bls11 = Run()
bls11.name = 'bc_layer11'
bls11.description = ''
bls11.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie10_layer11_rerun7.YIM/pp'
bls11.nc_dur = 20
bls11.nc_start_yr = 1
bls11.nc_end_yr = 120

# SM2.1 layer 14
bls14 = Run()
bls14.name = 'bc_layer14'
bls14.description = ''
bls14.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie10_layer14_rerun7.YIM/pp'
bls14.nc_dur = 20
bls14.nc_start_yr = 1
bls14.nc_end_yr = 120

# SM2.1 layer 16
bls16 = Run()
bls16.name = 'bc_layer16'
bls16.description = ''
bls16.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie10_layer16_rerun7.YIM/pp'
bls16.nc_dur = 20
bls16.nc_start_yr = 1
bls16.nc_end_yr = 120

# SM2.1 layer 18
bls18 = Run()
bls18.name = 'bc_layer18'
bls18.description = ''
bls18.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie10_layer18_rerun7.YIM/pp'
bls18.nc_dur = 20
bls18.nc_start_yr = 1
bls18.nc_end_yr = 120

# SM2.1 layer 20
bls20 = Run()
bls20.name = 'bc_layer20'
bls20.description = ''
bls20.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie10_layer20_rerun7.YIM/pp'
bls20.nc_dur = 20
bls20.nc_start_yr = 1
bls20.nc_end_yr = 120

# SM2.1 layer 22
bls22 = Run()
bls22.name = 'bc_layer22'
bls22.description = ''
bls22.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie10_layer22_rerun7.YIM/pp'
bls22.nc_dur = 20
bls22.nc_start_yr = 1
bls22.nc_end_yr = 120

# SM2.1 layer 24
bls24 = Run()
bls24.name = 'bc_layer24'
bls24.description = ''
bls24.direc_nc = '/archive/yim/sm2.1_fixed/SM2.1U_Control-1860_lm2_aie10_layer24_rerun7.YIM/pp'
bls24.nc_dur = 20
bls24.nc_start_yr = 1
bls24.nc_end_yr = 120


### Models ###
# bl_am2
bl_am2 = Model()
bl_am2.name = 'am2'
bl_am2.set_nc_grid('/home/s1h/nc_grid/a3gcm.am2.nc')
bl_am2.lat_name = 'lat'
bl_am2.lat = bl_am2.nc_grid.variables[bl_am2.lat_name][:]
bl_am2.lat_bounds_name = 'latb'
bl_am2.lat_bounds = bl_am2.nc_grid.variables[bl_am2.lat_bounds_name][:]
bl_am2.lon_name = 'lon'
bl_am2.lon = bl_am2.nc_grid.variables[bl_am2.lon_name][:]
bl_am2.lon_bounds_name = 'lonb'
bl_am2.lon_bounds = bl_am2.nc_grid.variables[bl_am2.lon_bounds_name][:]
bl_am2.sfc_area = grid_sfc_area(bl_am2.nc_grid)
bl_am2.level_name = 'level'
bl_am2.level = bl_am2.nc_grid.variables[bl_am2.level_name][:]
bl_am2.dt_name = 'average_DT'
bl_am2.direc_nc = '/archive/yim/am2.1_fixed/'
bl_am2.nc_dur = 20
bl_am2.nc_start_yr = 1
bl_am2.set_runs = ([bla_cont, bla11, bla14, bla16, bla18, bla20, bla22, bla24])


# bl_sm2
bl_sm2 = Model()
bl_sm2.name = 'sm2'
bl_sm2.set_nc_grid('/home/s1h/nc_grid/a3gcm.am2.nc')
bl_sm2.lat_name = 'lat'
bl_sm2.lat = bl_sm2.nc_grid.variables[bl_sm2.lat_name][:]
bl_sm2.lat_bounds_name = 'latb'
bl_sm2.lat_bounds = bl_sm2.nc_grid.variables[bl_sm2.lat_bounds_name][:]
bl_sm2.lon_name = 'lon'
bl_sm2.lon = bl_sm2.nc_grid.variables[bl_sm2.lon_name][:]
bl_sm2.lon_bounds_name = 'lonb'
bl_sm2.lon_bounds = bl_sm2.nc_grid.variables[bl_sm2.lon_bounds_name][:]
bl_sm2.sfc_area = grid_sfc_area(bl_sm2.nc_grid)
bl_sm2.level_name = 'level'
bl_sm2.level = bl_sm2.nc_grid.variables[bl_sm2.level_name][:]
bl_sm2.dt_name = 'average_DT'
bl_sm2.direc_nc = '/archive/yim/sm2.1_fixed/'
bl_sm2.nc_dur = 20
bl_sm2.nc_start_yr = 1
bl_sm2.set_runs = ([bls_cont, bls11, bls14, bls16, bls18, bls20, bls22, bls24])

### Proj ###
bc_layers = Proj()
bc_layers.name = 'bc_layers'
bc_layers.direc_out = '/archive/s1h/bc_layers/'
bc_layers.set_models([bl_am2, bl_sm2])
bc_layers.regions = [globe, nh, sh, tropics, wpwp, epac]
