import numpy as np
import netCDF4

from .constants import r_e
from .utils import _set_named_attr_dict

class Model(object):
    """Parameters of local data associated with a single climate model."""
    def __init__(self, name='', description='', proj=False, nc_grid_paths=(),
                 nc_dur=False, nc_start_yr=False, nc_end_yr=False,
                 default_yr_range=False, runs=(), default_runs=()):
        self.name = name
        self.description = description
        self.proj = proj
        self.nc_dur = nc_dur
        self.nc_grid_paths = nc_grid_paths
        self.nc_start_yr = nc_start_yr
        self.nc_end_yr = nc_end_yr
        self.default_yr_range = default_yr_range
        self.runs = runs
        self.default_runs = default_runs

        # Use the inputted names and netCDF filepath to create grid data.
        self._set_mult_nc_grid_attr()
        self._set_sfc_area()
        # Populate the runs dictionary with the specified list of runs.
        _set_named_attr_dict(self, 'runs', runs)
        _set_named_attr_dict(self, 'default_runs', default_runs)

    def __str__(self):
        return 'Model instance "' + self.name + '"'

    def _get_nc_grid(self):
        """Get the nc_grid of an aospy object."""
        nc_objs = []
        for path in self.nc_grid_paths:
            try:
                nc_obj = netCDF4.Dataset(path)
            except TypeError:
                nc_obj = netCDF4.MFDataset(path)
            except RuntimeError:
                raise RuntimeError('No netCDF file found at the specified '
                                   'directory: %s' % path)
            nc_objs.append(nc_obj)
        return tuple(nc_objs)

    def _get_nc_grid_attr(self, nc_grid, attr_name):
        """Get attribute from the nc_grid file(s)."""
        for nc in nc_grid:
            try:
                return nc.variables[attr_name][:]
            except KeyError:
                pass
        return None

    def _set_mult_nc_grid_attr(self):
        """
        Set multiple attrs from grid file given their names in the grid file.
        """
        grid_attrs = {
            'lat':         ('lat', 'latitude', 'LATITUDE', 'y'),
            'lat_bounds':  ('latb', 'lat_bnds', 'lat_bounds'),
            'lon':         ('lon', 'longitude', 'LONGITUDE', 'x'),
            'lon_bounds':  ('lonb', 'lon_bnds', 'lon_bounds'),
            'level':       ('level', 'lev', 'plev'),
            'time':        ('time', 'TIME'),
            'time_st':     ('average_T1',),
            'time_end':    ('average_T2',),
            'time_dur':    ('average_DT',),
            'time_bounds': ('time_bounds', 'time_bnds'),
            'sfc_area':    ('area', 'sfc_area'),
            'zsurf':       ('zsurf',),
            'land_mask':   ('land_mask',),
            'pk':          ('pk',),
            'bk':          ('bk',),
            'phalf':       ('phalf',),
            'pfull':       ('pfull',)
        }
        nc_grid = self._get_nc_grid()
        try:
            for name_int, names_ext in grid_attrs.iteritems():
                for name in names_ext:
                    setattr(self, name_int,
                            self._get_nc_grid_attr(nc_grid, name))
                    break
        except:
            raise
        finally:
            for nc in nc_grid:
                nc.close()

    def grid_sfc_area(self, lonb, latb, gridcenter=False):
        """Calculate surface area of each grid cell in a lon-lat grid."""
        def diff_latlon_bnds(array):
            if array.ndim == 1:
                return array[1:] - array[:-1]
            else:
                return array[:,1] - array[:,0]
        dlon = diff_latlon_bnds(lonb)
        # Must compute gridcell edges if given values at cell centers.
        if gridcenter:
            # Grid must be evenly spaced for algorithm to work.
            assert np.allclose(dlon[0], dlon)
            dlon = np.append(dlon, lonb[0] - lonb[-1] + 360.)
            dlat = diff_latlon_bnds(latb)
            assert np.allclose(dlat[0], dlat)
            # First and last array values done separately.
            lat = latb
            latb = np.append(lat[0] - 0.5*dlat[0], lat[:-1] + 0.5*dlat)
            latb = np.append(latb, lat[-1]+0.5*dlat[-1])

        dsinlat = diff_latlon_bnds(np.sin(np.deg2rad(latb)))
        X, Y = np.meshgrid(np.abs(dlon), np.abs(dsinlat))

        return X*Y*(r_e**2)*(np.pi/180.)

    def _set_sfc_area(self):
        """Set the 2D array holding the surface area of gridboxes."""
        if getattr(self, 'sfc_area', None) is not None:
            return
        else:
            try:
                sfc_area = self.grid_sfc_area(
                    self.lon_bounds, self.lat_bounds, gridcenter=False
                )
            except AttributeError:
                sfc_area = self.grid_sfc_area(
                    self.lon, self.lat, gridcenter=True
                )
            self.sfc_area = sfc_area

    def calc_levs_thick(self, levs):
        """
        Calculates the thickness, in Pa, of each pressure level.

        Assumes that the pressure values given are at the center of that model
        level, except for the lowest value (typically 1000 hPa), which is the
        bottom boundary. The uppermost level extends to 0 hPa.

        """
        # Bottom level extends from levs[0] to halfway betwen levs[0]
        # and levs[1].
        levs = to_pascal(levs)
        lev_thick = [0.5*(levs[0] - levs[1])]
        # Middle levels extend from halfway between [k-1], [k] and [k], [k+1].
        for k in range(1, levs.size-1):
            lev_thick.append(0.5*(levs[k-1] - levs[k+1]))
        # Top level extends from halfway between top two levels to 0 hPa.
        lev_thick.append(0.5*(levs[-2] + levs[-1]))
        # Convert to numpy array and from hectopascals (hPa) to Pascals (Pa).
        return np.array(lev_thick)

    def _set_levs_thick(self):
        """
        Set the 1D array holding the pressure thickness of the model levels.
        """
        if self.level:
            self.levs_thick = calc_levs_thick(self.level)
        else:
            self.levs_thick = None

def model_inst(model, parent_proj=False):
    """Convert string of an aospy.model name to an aospy.model instance."""
    if parent_proj and type(parent_proj) is not Proj:
        parent_proj = _proj_inst(parent_proj)
    if type(model) is Model:
        model_inst =  model
    elif type(model) is str:
        model_inst = parent_proj.models[model]
    elif type(model) in (list, tuple):
        model_inst = [_model_inst(mod, proj) for (mod, proj)
                      in zip(model, parent_proj)]
        if type(model) is tuple:
            model_inst = tuple(model_inst)
    if parent_proj:
        parent_proj = _proj_inst(parent_proj)
        try:
            model_inst.proj = parent_proj
        except AttributeError:
            pass
    return model_inst
