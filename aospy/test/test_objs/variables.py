from aospy.var import Var
import units
import calcs


p = Var(
    name='p',
    units=units.Pa,
    domain='atmos',
    description='Pressure of model half levels.',
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
ps = Var(
    name='ps',
    units=units.Pa,
    domain='atmos',
    description='Surface pressure.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
bk = Var(
    name='bk',
    units=units.Pa,
    domain='atmos_level',
    description='Sigma part of hybrid sigma coordinate.',
    def_time=False,
    def_vert=True,
    def_lat=False,
    def_lon=False,
    in_nc_grid=True
)
pk = Var(
    name='pk',
    units=units.Pa,
    domain='atmos_level',
    description='Pressure part of hybrid sigma coordinate.',
    def_time=False,
    def_vert=True,
    def_lat=False,
    def_lon=False,
    in_nc_grid=True
)
temp = Var(
    name='temp',
    alt_names=('ta',),
    units=units.K,
    domain='atmos',
    description='Air temperature.',
    def_time=True,
    def_vert='pfull',
    def_lat=True,
    def_lon=True,
    in_nc_grid=False,
    colormap='RdBu_r'
)
dp = Var(
    name='dp',
    domain='atmos',
    description=('Pressure thickness of model levels.  For data interpolated '
                 'to uniform pressure levels, this does not vary in time or '
                 'space.  For data on model native coordinates, this varies '
                 'in space and time due to the spatiotemporal variations in '
                 'surface pressure.'),
    variables=(ps, bk, pk, temp),
    def_time=True,
    def_vert=True,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False,
    func=calcs.dp,
    units=units.Pa,
)
olr = Var(
    name='olr',
    alt_names=('rlut',),
    units=units.W_m2,
    domain='atmos',
    description='All-sky outgoing longwave radiation at TOA.',
    def_time=True,
    def_vert=False,
    def_lat=True,
    def_lon=True,
    in_nc_grid=False
)
