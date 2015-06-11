from aospy import Var, Units
from aospy.constants import c_p, r_e, day2sec, sec2day
import aospy.calcs
            
unitless = Units(units='')
K = Units(units='K')
s=Units(
    units=r's',
    plot=r'day',
    plot_conv=sec2day
)
s1=Units(
    units=r's$^{-1}$',
    plot=r'day$^{-1}$',
    plot_conv=day2sec,
)
K_s1=Units(
    units=r'K s$^{-1}$',
    vert_int_plot=r'W m$^{-2}$',
    vert_int_plot_conv=c_p
)
m_s=Units(
    units=r'm s$^{-1}$',
    vert_int=r'kg m$^{-1}$ s$^{-1}$'
)
m_s2 = Units(
    units=r'm s$^{-2}$',
    vert_int=r'kg m$^{-1}$ s$^{-2}$'
)
kg_m2_s = Units(
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot=r'mm day$^{-1}$',
    plot_conv=day2sec
)
kg_m2_s_mass = Units(           # For vertical integrals of divergence.
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot=r'10$^{-2}$ kg m$^{-2}$ s',    
    plot_conv=1e2
)
W_m2 = Units(
    units=r'W m$^{-2}$',
    vert_int=r''
    )
J_kg1 = Units(
    units=r'J kg$^{-1}$',
    plot='K',
    plot_conv=1/c_p,
    vert_int='J m$^{-2}$',
    vert_int_plot='10$^6$ J m$^{-2}$',
    vert_int_plot_conv=1e-6
)
J_kg1_s1 = Units(
    units=r'J kg$^{-1}$ s$^{-1}$',
    plot='K day$^{-1}',
    plot_conv=day2sec/c_p,
    vert_int='W m$^{-2}$',
    vert_int_plot='W m$^{-2}$',
    vert_int_plot_conv=1
)
Pa=Units(
    units=r'Pa',
    plot=r'hPa',
    plot_conv=1e-2
)
hPa=Units(
    units=r'hPa',
)

alb_sfc = Var(
    name='alb_sfc',
    units='Surface albedo',
    plot_units='Surface albedo',
    plot_units_conv=1.,
    domain='atmos',
    description='Surface albedo.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
cld_amt = Var(
    name='cld_amt',
    alt_names=('cl',),
    units='Cloud fraction',
    plot_units='Cloud fraction',
    plot_units_conv=100.,
    domain='atmos',
    description='Cloud fraction at each level.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
divg = Var(
    name='divg',
    math_str=r"\nabla\cdot\vec{v}",
    units=s1,
    domain='atmos',
    description='Divergence.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False,
    cmap='RdBu'
)
esf = Var(
    name='esf',
    units=r'',
    plot_units=r'',
    plot_units_conv=1e-6,
    domain='atmos',
    description='Eddy streamfunction',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
evap = Var(
    name='evap',
    alt_names=('ET_mean', 'evspsbl'),
    math_str=r"$E$",
    units=kg_m2_s,
    domain='atmos',
    description='Surface evaporation',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
high_cld_amt = Var(
    name='high_cld_amt',
    units='Cloud fraction',
    plot_units='Cloud fraction',
    plot_units_conv=1.,
    domain='atmos',
    description='High cloud fraction.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
low_cld_amt = Var(
    name='low_cld_amt',
    units='Cloud fraction',
    plot_units='Cloud fraction',
    plot_units_conv=1.,
    domain='atmos',
    description='Low cloud fraction.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
lwdn_sfc = Var(
    name='lwdn_sfc',
    alt_names=('rlds',),
    math_str="$R^{LW\downarrow_{sfc}$",
    description='All-sky downwelling longwave radiation at the surface.',
    domain='atmos',
    units=W_m2,
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
lwdn_sfc_clr = Var(
    name='lwdn_sfc_clr',
    alt_names=('rldscs',),
    units=W_m2,
    domain='atmos',
    description='Clear-sky downwelling longwave radiation at the surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
lwup_sfc = Var(
    name='lwup_sfc',
    alt_names=('rlus',),
    units=W_m2,
    domain='atmos',
    description='All-sky upwelling longwave radiation at the surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
lwup_sfc_clr = Var(
    name='lwup_sfc_clr',
    alt_names=('rluscs',),
    units=W_m2,
    domain='atmos',
    description='Clear-sky upwelling longwave radiation at the surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
hght = Var(
    name='hght',
    alt_names=('zg',),
    units='meters',
    plot_units='meters',
    plot_units_conv=1.,
    domain='atmos',
    description='Geopotential height.',
    def_time=True,
    def_vert='phalf',
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
ice_wat = Var(
    name='ice_wat',
    units='kg/kg',
    plot_units=r'10$^{-3}$ g kg$^{-1}$',
    plot_units_conv=1.e6,
    domain='atmos',
    description='Cloud ice water specific humidity.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
liq_wat = Var(
    name='liq_wat',
    units='kg/kg',
    plot_units=r'10$^{-3}$ g kg$^{-1}$',
    plot_units_conv=1.e6,
    domain='atmos',
    description='Cloud liquid water specific humidity.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
mc = Var(
    name='mc',
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot_units=r'$M_c$, 10$^{-3}$ kg m$^{-2}$ s$^{-1}$',
    plot_units_conv=1e3,
    domain='atmos',
    description='Convective mass flux.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
mc_full = Var(
    name='mc_full',
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot_units=r'10$^{-3}$ kg m$^{-2}$ s$^{-1}$',
    plot_units_conv=1e3,
    domain='atmos',
    description='Convective mass flux at full levels.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
mc_half = Var(
    name='mc_half',
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot_units=r'10$^{-3}$ kg m$^{-2}$ s$^{-1}$',
    plot_units_conv=1e3,
    domain='atmos',
    description='Convective mass flux at half levels.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
mid_cld_amt = Var(
    name='mid_cld_amt',
    units='Cloud fraction',
    plot_units='Cloud fraction',
    plot_units_conv=1.,
    domain='atmos',
    description='Mid-level cloud fraction.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
olr = Var(
    name='olr',
    alt_names=('rlut',),
    units=W_m2,
    domain='atmos',
    description='All-sky outgoing longwave radiation at TOA.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
olr_clr = Var(
    name='olr_clr',
    alt_names=('rlutcs',),
    units=W_m2,
    domain='atmos',
    description='Clear-sky outgoing longwave radiation at TOA.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
omega = Var(
    name='omega',
    alt_names=('wap',),
    units=r'Pa s$^{-1}$',
    plot_units=r'hPa day$^{-1}$',
    plot_units_conv=24.*3600./100.,
    domain='atmos',
    description='Pressure vertical velocity.',
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
precip = Var(
    name='precip',
    alt_names=('pr', 'pre'),
    units = kg_m2_s,
    domain='atmos',
    description='Liquid precipitation reaching surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
prec_conv = Var(
    name='prec_conv',
    alt_names=('prc',),
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot_units=r'mm day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Liquid precipitation reaching surface from convection scheme.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
prec_ls = Var(
    name='prec_ls',
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot_units=r'mm day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Liquid precipitation reaching surface from large scale ascent.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
ps = Var(
    name='ps',
    units='Pa',
    plot_units='hPa',
    plot_units_conv=1e-2,
    domain='atmos',
    description='Surface pressure.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
pv = Var(
    name='pv',
    units=r's$^{-1}$',
    plot_units=r'10$^{-8}$ s$^{-1}$',
    plot_units_conv=1e8,
    domain='atmos',
    description='Potential vorticity',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
rh = Var(
    name='rh',
    alt_names=('hur',),
    units='',
    plot_units='Percent',
    plot_units_conv=1.,
    domain='atmos',
    description='Relative humidity',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
rh_ref = Var(
    name='rh_ref',
    units='',
    plot_units='Percent',
    plot_units_conv=1.,
    domain='atmos',
    description='Relative humidity at 2m above surface',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
shflx = Var(
    name='shflx',
    alt_names=('hfss',),
    units=W_m2,
    domain='atmos',
    description='Surface sensible heat flux into the atmosphere.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
slp = Var(
    name='slp',
    alt_names=('psl',),
    units=hPa,
    domain='atmos',
    description='Sea level pressure.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
snow_conv = Var(
    name='snow_conv',
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot_units=r'mm day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Snow reaching surface from convection scheme.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
snow_ls = Var(
    name='snow_ls',
    units=r'kg m$^{-2}$ s$^{-1}$',
    plot_units=r'mm day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Snow reaching surface from large scale ascent.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
soil_liq = Var(
    name='soil_liq',
    units=r'kg m$^{-2}$',
    plot_units=r'kg m$^{-2}$',
    plot_units_conv=1.,
    domain='land',
    description='Soil liquid in each level of land model',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
soil_moisture = Var(
    name='soil_moisture',
    alt_names=('water', 'water_soil',),
    units=r'kg m$^{-2}$',
    plot_units=r'kg m$^{-2}$',
    plot_units_conv=1.,
    domain='land',
    description='Mass of water in land bucket',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
sphum = Var(
    name='sphum',
    alt_names=('hus',),
    units='kg/kg',
    plot_units='g/kg',
    plot_units_conv=1e3,
    domain='atmos',
    description='Specific humidity.',
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
sst = Var(
    name='sst',
    alt_names=('ts',),
    units=K,
    domain='ocean',
    description='Sea surface temperature.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swdn_sfc = Var(
    name='swdn_sfc',
    alt_names=('rsds',),
    units=W_m2,
    domain='atmos',
    description='All-sky downwelling shortwave radiation at the surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swdn_sfc_clr = Var(
    name='swdn_sfc_clr',
    alt_names=('rsdscs',),
    units=W_m2,
    domain='atmos',
    description='Clear-sky downwelling shortwave radiation at the surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swup_sfc = Var(
    name='swup_sfc',
    alt_names=('rsus',),
    units=W_m2,
    domain='atmos',
    description='All-sky upwelling shortwave radiation at the surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swup_sfc_clr = Var(
    name='swup_sfc_clr',
    alt_names=('rsuscs',),
    units=W_m2,
    domain='atmos',
    description='Clear-sky upwelling shortwave radiation at the surface.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swdn_toa = Var(
    name='swdn_toa',
    alt_names=('rsdt',),
    units=W_m2,
    domain='atmos',
    description='Downwelling shortwave radiation at TOA.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swdn_toa_clr = Var(
    name='swdn_toa_clr',
    alt_names=('rsdtcs',),
    units=W_m2,
    domain='atmos',
    description='Downwelling shortwave radiation at TOA.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swup_toa = Var(
    name='swup_toa',
    alt_names=('rsut',),
    units=W_m2,
    domain='atmos',
    description='All-sky Upwelling shortwave radiation at TOA.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
swup_toa_clr = Var(
    name='swup_toa_clr',
    alt_names=('rsutcs',),
    units=W_m2,
    domain='atmos',
    description='Clear-sky upwelling shortwave radiation at TOA.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
t_surf = Var(
    name='t_surf',
    alt_names=('tas', 'tmp'),
    units=K,
    domain='atmos',
    description='Surface air temperature.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
tdt_conv = Var(
    name='tdt_conv',
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Convective heating rate.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False,
    valid_range=(-20, 20)
)
tdt_ls = Var(
    name='tdt_ls',
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Large-scale heating rate.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
tdt_lw = Var(
    name='tdt_lw',
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='All-sky longwave heating rate.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
tdt_lw_clr = Var(
    name='tdt_lw_clr',
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Clear-sky longwave heating rate.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
tdt_sw = Var(
    name='tdt_sw',
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='All-sky shortwave heating rate.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
tdt_sw_clr = Var(
    name='tdt_sw_clr',
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Clear-sky shortwave heating rate.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
tdt_vdif = Var(
    name='tdt_vdif',
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
    domain='atmos',
    description='Temperature tendency from vertical diffusion.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
temp = Var(
    name='temp',
    alt_names=('ta',),
    units=K,
    domain='atmos',
    description='Air temperature.',
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
tot_cld_amt = Var(
    name='tot_cld_amt',
    alt_names=('clt',),
    units='Cloud fraction',
    plot_units='Cloud fraction',
    plot_units_conv=1.,
    domain='atmos',
    description='Total cloud fraction.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
ucomp = Var(
    name='ucomp',
    alt_names=('ua',),
    units=r'm s$^{-1}$',
    plot_units=r'm s$^{-1}$',
    plot_units_conv=1.,
    domain='atmos',
    description='Eastward velocity.',
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
vcomp = Var(
    name='vcomp',
    alt_names=('va',),
    units=r'm s$^{-1}$',
    plot_units=r'm s$^{-1}$',
    plot_units_conv=1.,
    domain='atmos',
    description='Northward velocity.',
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
vort = Var(
    name='vort',
    units=r's$^{-1}$',
    plot_units=r'10$^{-6}$ s$^{-1}$',
    plot_units_conv=1e6,
    domain='atmos',
    description='Vorticity',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
wvp = Var(
    name='wvp',
    alt_names=('WVP', 'prw'),
    units=r'kg m$^{-2}$',
    plot_units=r'kg m$^{-2}$',
    plot_units_conv=1,
    domain='atmos',
    description='Water vapor path',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
# Grid variables.
lat = Var(
    name='lat',
    units='Degrees',
    domain=None,
    description='Latitude.',
    def_time=False,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    in_nc_grid=True
)
lon = Var(
    name='lon',
    units='Degrees',
    domain=None,
    description='Longitude.',
    def_time=False,
    def_vert=False,
    def_lat=False,
    def_lon=True,
    in_nc_grid=True
)
level = Var(
    name='level',
    units=hPa,
    domain='atmos',
    description='Pressure level.',
    def_time=False,
    def_vert=True,
    def_lat=False,
    def_lon=False,
    in_nc_grid=True
)
pk = Var(
    name='pk',
    units=Pa,
    domain='atmos_level',
    description='Pressure part of hybrid sigma coordinate.',
    def_time=False,
    def_vert=True,
    def_lat=False,
    def_lon=False,
    in_nc_grid=True
)
bk = Var(
    name='bk',
    units=Pa,
    domain='atmos_level',
    description='Sigma part of hybrid sigma coordinate.',
    def_time=False,
    def_vert=True,
    def_lat=False,
    def_lon=False,
    in_nc_grid=True
)
sfc_area = Var(
    name='sfc_area',
    units=r'm$^2$',
    domain=None,
    description='Grid surface area.',
    def_time=False,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
### Calculations involving one or more model-native variables.
aht = Var(name='aht',
          domain='atmos',
          description='Total northward atmospheric heat transport.',
          vars=(swdn_toa, swup_toa, olr, swup_sfc, swdn_sfc, lwup_sfc, lwdn_sfc,
                shflx, evap, snow_ls, snow_conv, sfc_area),
          def_time=True,
          def_vert=False,
          def_lat=True,
          def_lon=False,
          func=aospy.calcs.aht,
          units='Watts',
          plot_units='PW',
          plot_units_conv=1e-15
)
albedo = Var(
    name='albedo',
    domain='atmos',
    description='Column albedo: fraction of incoming SW at TOA reflectedback to space.',
    vars=(swdn_toa, swup_toa),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.albedo,
    units='Column albedo',
    plot_units='Column albedo',
    plot_units_conv=1.
)
sfc_albedo = Var(
    name='sfc_albedo',
    domain='atmos',
    description='Surface albedo: Downward SW / Upward SW at surface.',
    vars=(swdn_sfc, swup_sfc),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_albedo,
    units='Surface albedo',
    plot_units='Surface albedo',
    plot_units_conv=1.
)
ang_mom = Var(
    name='ang_mom',
    domain='atmos',
    description='Angular momentum per unit mass.',
    vars=(lat, ucomp),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.ang_mom,
    units='',
    plot_units='',
    plot_units_conv=1.
)
bowen_ratio = Var(
    name='bowen_ratio',
    domain='atmos',
    description='Bowen Ratio: ratio of surface sensible to latent heat flux',
    vars=(shflx, evap),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.bowen_ratio,
    units='',
    plot_units='',
    plot_units_conv=1.
)
column_energy = Var(
    name='column_energy',
    domain='atmos',
    description='Net energy flux into atmosphere at surface and at TOA.',
    vars=(swdn_toa, swup_toa, olr, swup_sfc, swdn_sfc, lwup_sfc, lwdn_sfc,
          shflx, evap),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.column_energy,
    units=W_m2
)
cre_lw = Var(
    name='cre_lw',
    domain='atmos',
    description='Cloudy-sky outgoing longwave radiation at TOA.',
    vars=(olr, olr_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.cre_lw,
    units=W_m2
)
cre_net = Var(
    name='cre_net',
    domain='atmos',
    description='Net cloudy-sky top-of-atmosphere radiative flux, signed positive downwards.',
    # vars=(swdn_toa, swup_toa, olr, swdn_toa_clr, swup_toa_clr, olr_clr),
    vars=(swup_toa, olr, swup_toa_clr, olr_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.cre_net,
    units=W_m2
)
cre_sw = Var(
    name='cre_sw',
    domain='atmos',
    description='Net cloudy-sky top-of-atmosphere shortwave radiative flux, signed positive downwards.',
    # vars=(swdn_toa, swup_toa, swdn_toa_clr, swup_toa_clr),
    vars=(swup_toa, swup_toa_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.cre_sw,
    units=W_m2
)
descent_tot = Var(
    name='descent_tot',
    domain='atmos',
    description='Total vertical pressure velocity (i.e. signed positive moving vertically downwards) per unit surface area. Includes both convective and large scale.',
    vars=(omega, mc),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.descent_tot,
    units=r'Pa s$^{-1}$',
    plot_units=r'hPa day$^{-1}$',
    plot_units_conv=24.*3600./100.
)
divg_mass_bal = Var(
    name='divg_mass_bal',
    domain='atmos',
    description='Divergence mass balanced in each column.',
    vars=(divg, 'dp'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_vert_int_bal,
    units=r's$^{-1}$',
    plot_units=r'day$^{-1}$',
    plot_units_conv=24.*3600.
)
dry_static_stab = Var(
    name='dry_static_stab',
    domain='atmos',
    description=('Local upper minus lower level DSE for fixed lower level and '
                 'varying upper level over all levels.'),
    vars=('temp', 'hght', 'p'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.dry_static_stab,
    units=K
)
dse = Var(
    name='dse',
    domain='atmos',
    description='Dry static energy.',
    vars=(temp, hght),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.dse,
    units=J_kg1
)
dse_horiz_advec = Var(
    name='dse_horiz_advec',
    domain='atmos',
    description='Horizontal advection of dry static energy.',
    vars=(temp, hght, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.dse_horiz_advec,
    units=J_kg1_s1
)
dse_times_horiz_divg = Var(
    name='dse_times_horiz_divg',
    domain='atmos',
    description='Horizontal mass flux divergence times dry static energy.',
    vars=(temp, hght, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.dse_times_horiz_divg,
    units=J_kg1_s1
)
dse_horiz_flux_divg = Var(
    name='dse_horiz_flux_divg',
    domain='atmos',
    description='Horizontal flux divergence of dry static energy.',
    vars=(temp, hght, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.dse_horiz_flux_divg,
    units=J_kg1_s1
)
dse_horiz_advec_divg_sum = Var(
    name='dse_horiz_advec_divg_sum',
    domain='atmos',
    description='Sum of DSE horizontal divergence and advection.',
    vars=(temp, hght, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.dse_horiz_advec_divg_sum,
    units=J_kg1_s1
)
du_dx = Var(
    name='du_dx',
    domain='atmos',
    description='',
    vars=(ucomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.d_dx_from_latlon,
    units=s1
)
dv_dy = Var(
    name='dv_dy',
    domain='atmos',
    description='',
    vars=(vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.d_dy_from_latlon,
    units=s1
)
equiv_pot_temp = Var(
    name='equiv_pot_temp',
    domain='atmos',
    description='Equivalent potential temperature.',
    vars=(temp, level, sphum),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.equiv_pot_temp,
    units=K
)
esf = Var(
    name='esf',
    domain='atmos',
    description='Eddy streamfunction.',
    vars=(lat, level, vcomp),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=None,
    units=r'kg s$^{-1}$',
    plot_units=r'10$^{10}$ kg s$^{-1}$',
    plot_units_conv=1.
)
evap_frac = Var(
    name='evap_frac',
    domain='atmos',
    description='Evaporative fraction, i.e. ratio of LH to (LH+SH).',
    vars=(evap, shflx),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.evap_frac,
    units=r'kg s$^{-1}$',
    plot_units=r'10$^{10}$ kg s$^{-1}$',
    plot_units_conv=1.
)
fmse = Var(
    name='fmse',
    domain='atmos',
    description='Frozen moist static energy.',
    vars=(temp, hght, sphum, ice_wat),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.fmse,
    units=J_kg1
)
gms_change_est = Var(
    name='gms_change_est',
    domain='atmos',
    description='Gross moist stability estimated as near surface MSE at ITCZ minus at the local latitude.',
    vars=(temp, temp, sphum, precip, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_change_est,
    units=K
)
gms_change_est2 = Var(
    name='gms_change_est2',
    domain='atmos',
    description='Gross moist stability estimated as near surface MSE at ITCZ minus at the local latitude.',
    vars=(temp, temp, sphum, precip, level, lat),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_change_est2,
    units=K
)
gms_h01 = Var(
    name='gms_h01',
    domain='atmos',
    description='Gross moist stability estimated as near surface MSE at ITCZ minus at the local latitude.',
    vars=(temp, hght, sphum, precip, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_h01,
    units=K
)
gms_h01est = Var(
    name='gms_h01est',
    domain='atmos',
    description='Gross moist stability estimated as near surface MSE at ITCZ minus at the local latitude, but neglecting the geopotential term.',
    vars=(temp, sphum, precip, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_h01est,
    units=K
)
gms_h01est2 = Var(
    name='gms_h01est2',
    domain='atmos',
    description='Gross moist stability estimated as MSE aloft at ITCZ minus near surface MSE at the local latitude.',
    vars=(temp, hght, sphum, precip, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_h01est2,
    units=K
)
gms_moc = Var(
    name='gms_moc',
    domain='atmos',
    description='Gross moist stability using only MMC MSE transport',
    vars=(temp, hght, sphum, precip, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_moc,
    units=K
)
gms_msf = Var(
    name='gms_msf',
    domain='atmos',
    description='Gross moist stability using MMC plus stationary eddy MSE transport.',
    vars=(temp, hght, sphum, precip, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_msf,
    units=K
)
gms_up_low = Var(
    name='gms_up_low',
    domain='atmos',
    description='Gross moist stability estimated as uppper minus lower level MSE.',
    vars=(temp, hght, sphum, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.gms_up_low,
    units=K
)
gms_each_level = Var(
    name='gms_each_level',
    domain='atmos',
    description='Local upper minus lower level MSE for fixed lower level and varying upper level over all levels.',
    vars=(temp, hght, sphum, 'p'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.gms_each_level,
    units=K
)
gross_dry_stab = Var(
    name='gross_dry_stab',
    domain='atmos',
    description='Gross dry stability.',
    vars=(temp, hght, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.gross_dry_stab,
    units='J kg$^{-1}$',
    plot_units='Kelvin',
    plot_units_conv=1./c_p
)
gross_moist_stab = Var(
    name='gross_moist_stab',
    domain='atmos',
    description='Gross dry stability.',
    vars=(temp, hght, sphum, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.gross_moist_stab,
    units='J kg$^{-1}$',
    plot_units='Kelvin',
    plot_units_conv=1./c_p
)
gross_moist_strat = Var(
    name='gross_moist_strat',
    domain='atmos',
    description='',
    vars=(sphum, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.gross_moist_strat,
    units='J kg$^{-1}',
    plot_units='Kelvin',
    plot_units_conv=1./c_p
)
hght_horiz_advec = Var(
    name='hght_horiz_advec',
    domain='atmos',
    description='Horizontal advection of geopotential height.',
    vars=(hght, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.horiz_advec,
    units=J_kg1_s1
)
hght_times_horiz_divg = Var(
    name='hght_times_horiz_divg',
    domain='atmos',
    description='Horizontal mass flux divergence times geopotential height.',
    vars=(hght, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_times_horiz_divg_mass_bal,
    units=J_kg1_s1
)
hght_horiz_flux_divg = Var(
    name='hght_horiz_flux_divg',
    domain='atmos',
    description='Horizontal flux divergence of geopotential height.',
    vars=(hght, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_horiz_flux_divg,
    units=J_kg1_s1
)
hght_horiz_advec_divg_sum = Var(
    name='hght_horiz_advec_divg_sum',
    domain='atmos',
    description='Sum of geopotential height horizontal divergence and advection.',
    vars=(hght, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_horiz_advec_divg_sum,
    units=J_kg1_s1
)
horiz_divg = Var(
    name='horiz_divg',
    domain='atmos',
    description='',
    vars=(ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.horiz_divg,
    units=s1,
    cmap='RdBu'
)
horiz_divg_mass_bal = Var(
    name='horiz_divg_mass_bal',
    domain='atmos',
    description='',
    vars=(ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.horiz_divg_mass_bal,
    units=s1,
    cmap='RdBu'
)
horiz_divg_vert_int_max = Var(
    name='horiz_divg_vert_int_max',
    domain='atmos',
    description=('Integral of horizontal divergence from surface to level '
                 'that maximizes its magnitude'),
    vars=(ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.horiz_divg_vert_int_max,
    units=kg_m2_s_mass
)
moist_static_stab = Var(
    name='moist_static_stab',
    domain='atmos',
    description='',
    vars=('temp', 'p', 'sphum'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.moist_static_stab,
    units=K,
    cmap='RdBu_r'
)
mse = Var(
    name='mse',
    domain='atmos',
    description='Moist static energy.',
    vars=(temp, hght, sphum),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse,
    units=J_kg1
)
mse_horiz_advec = Var(
    name='mse_horiz_advec',
    domain='atmos',
    description='Horizontal advection of moist static energy.',
    vars=(temp, hght, sphum, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse_horiz_advec,
    units=J_kg1_s1
)
mse_times_horiz_divg = Var(
    name='mse_times_horiz_divg',
    domain='atmos',
    description='Horizontal mass flux divergence times moist static energy.',
    vars=(temp, hght, sphum, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse_times_horiz_divg,
    units=J_kg1_s1
)
mse_horiz_flux_divg = Var(
    name='mse_horiz_flux_divg',
    domain='atmos',
    description='Horizontal flux divergence of moist static energy.',
    vars=(temp, hght, sphum, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse_horiz_flux_divg,
    units=J_kg1_s1
)
mse_horiz_advec_divg_sum = Var(
    name='mse_horiz_advec_divg_sum',
    domain='atmos',
    description='Sum of MSE horizontal divergence and advection.',
    vars=(temp, hght, sphum, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse_horiz_advec_divg_sum,
    units=J_kg1_s1
)
mse_vert_flux_divg = Var(
    name='mse_vert_flux_divg',
    domain='atmos',
    description='Vertical flux divergence of moist static energy.',
    vars=(temp, hght, sphum, omega, 'p'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse_vert_flux_divg,
    units=J_kg1_s1
)
mse_vert_advec = Var(
    name='mse_vert_advec',
    domain='atmos',
    description='Moist static energy times vertical divergence.',
    vars=(temp, hght, sphum, omega, 'p'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse_vert_advec,
    units=J_kg1_s1
)
mse_times_vert_divg = Var(
    name='mse_times_vert_divg',
    domain='atmos',
    description='Vertical advection of moist static energy.',
    vars=(temp, hght, sphum, omega, 'p', 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.mse_times_vert_divg,
    units=J_kg1_s1
)
msf = Var(
    name='msf',
    domain='atmos',
    description='Eulerian meridional mass streamfunction.',
    vars=(lat, level, vcomp),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.msf,
    units=r'kg s$^{-1}$',
    plot_units=r'10$^{10}$ kg s$^{-1}$',
    plot_units_conv=1e-10,
)
mass_flux = Var(
    name='mass_flux',
    domain='atmos',
    description=('Mass flux: Eulerian meridional mass streamfunction '
                 'integrated to the level of its maximum magnitude.'),
    vars=(lat, level, vcomp),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.msf_max,
    units=r'kg s$^{-1}$',
    plot_units=r'10$^{10}$ kg s$^{-1}$',
    plot_units_conv=1e-10,
)
p_minus_e = Var(
    name='p-e',
    domain='atmos',
    description='Precipitation minus evaporation',
    vars=(precip, evap),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.p_minus_e,
    units=kg_m2_s
)
pot_temp = Var(
    name='pot_temp',
    domain='atmos',
    description='Potential temperature.',
    vars=(temp, 'p'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.pot_temp,
    units=K
)
prec_conv_frac = Var(
    name='prec_conv_frac',
    domain='atmos',
    description='Fraction of liquid precipitation reaching surface originating from convection scheme (as opposed to large-scale condensation.',
    vars=(prec_conv, precip),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.prec_conv_frac,
    units=r'',
    plot_units=r'',
    plot_units_conv=1.
)
q_horiz_advec = Var(
    name='q_horiz_advec',
    domain='atmos',
    description='Horizontal advection of specific humidity.',
    vars=(sphum, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.q_horiz_advec,
    units=s1
)
q_times_horiz_divg = Var(
    name='q_times_horiz_divg',
    domain='atmos',
    description='Horizontal flux divergence times specific humidity.',
    vars=(sphum, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.q_times_horiz_divg,
    units=s1
)
q_horiz_flux_divg = Var(
    name='q_horiz_flux_divg',
    domain='atmos',
    description='Horizontal flux divergence of specific humidity.',
    vars=(sphum, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.q_horiz_flux_divg,
    units=s1
)
q_horiz_advec_divg_sum = Var(
    name='q_horiz_advec_divg_sum',
    domain='atmos',
    description='Sum of Q horizontal divergence and advection.',
    vars=(sphum, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_horiz_advec_divg_sum,
    units=s1
)
q_vert_flux_divg = Var(
    name='q_vert_flux_divg',
    domain='atmos',
    description='Vertical flux divergence of specific humidity.',
    vars=(sphum, omega, 'p'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_vert_flux_divg,
    units=s1
)
q_vert_advec = Var(
    name='q_vert_advec',
    domain='atmos',
    description='Vertical advection of specific humidity.',
    vars=(sphum, omega, 'p'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.vert_advec,
    units=s1
)
q_times_vert_divg = Var(
    name='q_times_vert_divg',
    domain='atmos',
    description='Specific humidity times vertical divergence.',
    vars=(sphum, omega, 'p', 'dp'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_times_vert_divg_mass_bal,
    units=s1
)
qu = Var(
    name='qu',
    domain='atmos',
    description='Zonal specific humidity flux.',
    vars=(sphum, ucomp),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.qu,
    units=r'(kg/kg)*(m/s)',
    plot_units=r'(g kg$^{-1}$)(m s$^{-1}$)',
    plot_units_conv=1e3
)
qv = Var(
    name='qv',
    domain='atmos',
    description='Meridional specific humidity flux.',
    vars=(sphum, vcomp),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.qv,
    units=r'(kg/kg)*(m/s)',
    plot_units=r'(g kg$^{-1}$)(m s$^{-1}$)',
    plot_units_conv=1e3
)
sfc_albedo = Var(
    name='sfc_albedo',
    domain='atmos',
    description='Surface albedo, masked where downwelling SW at surface is zero.',
    vars=(swdn_sfc, swup_sfc),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_albedo,
    units='Surface albedo',
    plot_units='Surface albedo',
    plot_units_conv=1.
)
sfc_energy = Var(
    name='sfc_energy',
    domain='atmos',
    description='Net surface energy flux, signed positive upwards.',
    vars=(swup_sfc, swdn_sfc, lwup_sfc, lwdn_sfc, shflx, evap),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_energy,
    units=W_m2
)
sfc_lw = Var(
    name='sfc_lw',
    domain='atmos',
    description='Net all-sky longwave radiative flux at surface, signed positive upwards.',
    vars=(lwup_sfc, lwdn_sfc),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_lw,
    units=W_m2
)
sfc_lw_cld = Var(
    name='sfc_lw_cld',
    domain='atmos',
    description='Net cloudy-sky longwave radiative flux at surface, signed positive upwards.',
    vars=(lwup_sfc, lwup_sfc_clr, lwdn_sfc, lwdn_sfc_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_lw_cld,
    units=W_m2
)
sfc_rad = Var(
    name='sfc_rad',
    domain='atmos',
    description='Net all-sky surface radiative flux, signed positive upwards.',
    vars=(swup_sfc, swdn_sfc, lwup_sfc, lwdn_sfc),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_rad,
    units=W_m2
)
sfc_rad_cld = Var(
    name='sfc_rad_cld',
    domain='atmos',
    description='Net cloudy-sky surface radiative flux, signed positive upwards.',
    vars=(swup_sfc, swup_sfc_clr, swdn_sfc, swdn_sfc_clr,
          lwup_sfc, lwup_sfc_clr, lwdn_sfc, lwdn_sfc_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_rad_cld,
    units=W_m2
)
sfc_sw = Var(
    name='sfc_sw',
    domain='atmos',
    description='Net all-sky shortwave radiative flux at surface, signed positive upwards.',
    vars=(swup_sfc, swdn_sfc),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_sw,
    units=W_m2
)
sfc_sw_cld = Var(
    name='sfc_sw_cld',
    domain='atmos',
    description='Net cloudy-sky shortwave radiative flux at surface, signed positive upwards.',
    vars=(swup_sfc, swup_sfc_clr, swdn_sfc, swdn_sfc_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.sfc_sw_cld,
    units=W_m2
)
temp_horiz_advec = Var(
    name='temp_horiz_advec',
    domain='atmos',
    description='Horizontal advection of temperature.',
    vars=(temp, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.horiz_advec,
    units=K_s1
)
temp_times_horiz_divg = Var(
    name='temp_times_horiz_divg',
    domain='atmos',
    description='Horizontal mass flux divergence times temperature.',
    vars=(temp, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_times_horiz_divg_mass_bal,
    units=K_s1
)
temp_horiz_flux_divg = Var(
    name='temp_horiz_flux_divg',
    domain='atmos',
    description='Horizontal flux divergence of temperature.',
    vars=(temp, ucomp, vcomp, lat, lon, r_e),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_horiz_flux_divg,
    units=K_s1
)
temp_horiz_advec_divg_sum = Var(
    name='temp_horiz_advec_divg_sum',
    domain='atmos',
    description='Sum of temperature horizontal divergence and advection.',
    vars=(temp, ucomp, vcomp, lat, lon, r_e, 'dp'),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.field_horiz_advec_divg_sum,
    units=K_s1
)
temp_vert_advec = Var(
    name='temp_vert_advec',
    domain='atmos',
    description='Vertical advection of temperature.',
    vars=(temp, omega, 'p'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.vert_advec,
    units=s1
)
tdt_diab = Var(
    name='tdt_diab',
    domain='atmos',
    description='Net temperature tendency from diabatic terms: LW, SW, convective, large-scale, and vertical diffusion.',
    #tdt_diab.vars=(tdt_lw, tdt_sw, tdt_conv, tdt_ls, tdt_vdif),
    vars=(tdt_lw, tdt_sw, tdt_conv, tdt_ls),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.tdt_diab,
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
)
tdt_lw_cld = Var(
    name='tdt_lw_cld',
    domain='atmos',
    description='Cloudy-sky temperature tendency from longwave radiation.' ,
    vars=(tdt_lw, tdt_lw_clr),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.tdt_lw_cld,
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
)
tdt_sw_cld = Var(
    name='tdt_sw_cld',
    domain='atmos',
    description='Cloudy-sky temperature tendency from shortwave radiation.' ,
    vars=(tdt_sw, tdt_sw_clr),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.tdt_sw_cld,
    units=r'Kelvin s$^{-1}$ ',
    plot_units=r'Kelvin day$^{-1}$',
    plot_units_conv=24.*3600.,
)
toa_rad = Var(
    name='toa_rad',
    domain='atmos',
    description='Net top-of-atmosphere radiative flux, signed positive downwards.',
    vars=(swdn_toa, swup_toa, olr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.toa_rad,
    units=W_m2
)
toa_rad_clr = Var(
    name='toa_rad_clr',
    domain='atmos',
    description='Clear-sky TOA radiative flux, positive downwards.',
    vars=(swdn_toa_clr, swup_toa_clr, olr_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.toa_rad_clr,
    units=W_m2
)
toa_sw = Var(
    name='toa_sw',
    domain='atmos',
    description='Net all-sky top-of-atmosphere shortwave radiative flux, signed positive downwards.',
    vars=(swdn_toa, swup_toa),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.toa_sw,
    units=W_m2
)
total_gms = Var(
    name='total_gms',
    domain='atmos',
    description='Total gross moist stability, i.e. GMS using MMC plus stationary and transient eddy MSE transports.',
    vars=(temp, hght, sphum, precip, 'p'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=False,
    func=aospy.calcs.total_gms,
    units=K
)
toa_swup_cld = Var(
    name='toa_swup_cld',
    domain='atmos',
    description='Upwelling cloudy-sky top-of-atmosphere shortwave radiative flux',
    vars=(swup_toa, swup_toa_clr),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.toa_swup_cld,
    units=W_m2
)
vert_divg = Var(
    name='vert_divg',
    domain='atmos',
    description='',
    vars=(omega, 'p'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.vert_divg,
    units=s1
)
vert_divg_mass_bal = Var(
    name='vert_divg_mass_bal',
    domain='atmos',
    description='',
    vars=(omega, 'p', 'dp'),
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.vert_divg_mass_bal,
    units=s1
)
vert_divg_vert_int_max = Var(
    name='vert_divg_vert_int_max',
    domain='atmos',
    description=('Integral of vertical divergence from surface to level '
                 'that maximizes its magnitude'),
    vars=(omega, 'p', 'dp'),
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.vert_divg_vert_int_max,
    units=kg_m2_s_mass
)
virt_pot_temp = Var(
    name='virt_pot_temp',
    domain='atmos',
    description='Virtual potential temperature.',
    vars=(temp, 'p', sphum, liq_wat),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    func=aospy.calcs.virt_pot_temp,
    units=K
)

master_vars_list = [
    alb_sfc, cld_amt, divg, esf, evap, hght, high_cld_amt, ice_wat, liq_wat,
    low_cld_amt, lwdn_sfc, lwdn_sfc_clr, lwup_sfc, lwup_sfc_clr, mc, mc_full,
    mc_half, mid_cld_amt, olr, olr_clr, omega, precip, prec_conv, prec_ls, ps,
    pv, rh, rh_ref, shflx, slp, snow_conv, snow_ls, soil_liq, soil_moisture,
    sphum, sst, swdn_sfc, swdn_sfc_clr, swup_sfc, swup_sfc_clr, swdn_toa,
    swdn_toa_clr, swup_toa, swup_toa_clr, t_surf, tdt_conv, tdt_ls, tdt_lw,
    tdt_lw_clr, tdt_sw, tdt_sw_clr, tdt_vdif, temp, tot_cld_amt, ucomp, vcomp,
    vort, wvp, lat, lon, level, pk, bk, sfc_area, aht, albedo, sfc_albedo,
    ang_mom, bowen_ratio, column_energy, cre_net, cre_lw, cre_sw, descent_tot,
    equiv_pot_temp, esf, evap_frac, fmse, gms_change_est,
    gms_change_est2, gms_h01, gms_h01est, gms_h01est2, gms_moc, gms_msf,
    gms_up_low, gms_each_level, gross_dry_stab, gross_moist_stab,
    dry_static_stab, total_gms, dse, horiz_divg, moist_static_stab,
    gross_moist_strat, mse, mse_horiz_advec, mse_times_horiz_divg,
    mse_vert_advec, msf, mass_flux, p_minus_e, pot_temp, prec_conv_frac,
    sfc_albedo, sfc_energy, sfc_lw, sfc_lw_cld, sfc_rad, sfc_rad_cld, sfc_sw,
    sfc_sw_cld, tdt_diab, tdt_lw_cld, tdt_sw_cld, toa_rad, toa_rad_clr, toa_sw,
    toa_swup_cld, vert_divg, virt_pot_temp, divg_mass_bal, horiz_divg_mass_bal,
    vert_divg_mass_bal, mse_horiz_flux_divg, q_horiz_advec, q_vert_advec,
    q_times_horiz_divg, q_horiz_flux_divg, qu, qv, du_dx, dv_dy,
    dse_horiz_flux_divg, dse_times_horiz_divg, dse_horiz_advec,
    temp_horiz_flux_divg, temp_times_horiz_divg, temp_horiz_advec,
    hght_horiz_flux_divg, hght_times_horiz_divg, hght_horiz_advec,
    mse_horiz_advec_divg_sum, q_horiz_advec_divg_sum,
    temp_horiz_advec_divg_sum, hght_horiz_advec_divg_sum,
    dse_horiz_advec_divg_sum, q_times_vert_divg, q_vert_flux_divg,
    mse_times_vert_divg, mse_vert_flux_divg, horiz_divg_vert_int_max,
    vert_divg_vert_int_max, temp_vert_advec
]
