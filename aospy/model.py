"""model.py: Model class of aospy for storing attributes of a GCM."""
import numpy as np
import netCDF4

from .constants import r_e
from .utils import dict_name_keys, level_thickness


class Model(object):
    """Parameters of local data associated with a single climate model."""
    def __init__(self, name='', description='', proj=False, nc_grid_paths=(),
                 nc_dur=False, nc_start_yr=False, nc_end_yr=False,
                 nc_start_month=False, nc_end_month=False,
                 default_yr_range=False, runs={}, default_runs={},
                 load_grid_data=False):
        self.name = name
        self.description = description
        self.proj = proj
        self.nc_dur = nc_dur
        self.nc_grid_paths = nc_grid_paths
        self.nc_start_yr = nc_start_yr
        self.nc_end_yr = nc_end_yr
        self.nc_start_month = nc_start_month
        self.nc_end_month = nc_end_month
        self.default_yr_range = default_yr_range
        self.runs = dict_name_keys(runs)
        self.default_runs = dict_name_keys(default_runs)

        self.grid_data_is_set = False
        if load_grid_data:
            self.set_grid_data()
            self.grid_data_is_set = True

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
            'lat':         ('lat', 'latitude', 'LATITUDE', 'y', 'yto'),
            'lat_bounds':  ('latb', 'lat_bnds', 'lat_bounds'),
            'lon':         ('lon', 'longitude', 'LONGITUDE', 'x', 'xto'),
            'lon_bounds':  ('lonb', 'lon_bnds', 'lon_bounds'),
            'level':       ('level', 'lev', 'plev'),
            'time':        ('time', 'TIME'),
            'time_st':     ('average_T1',),
            'time_end':    ('average_T2',),
            'time_dur':    ('average_DT',),
            'time_bounds': ('time_bounds', 'time_bnds'),
            'year':        ('yr', 'year'),
            'month':       ('mo', 'month'),
            'day':         ('dy', 'day'),
            'sfc_area':    ('area', 'sfc_area'),
            'zsurf':       ('zsurf',),
            'land_mask':   ('land_mask',),
            'pk':          ('pk',),
            'bk':          ('bk',),
            'phalf':       ('phalf',),
            'pfull':       ('pfull',),
            'nrecords':    ('nrecords',),
            'idim':        ('idim',),
            'fill_value':  ('fill_value')
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

    def _set_levs_thick(self):
        """
        Set the 1D array holding the pressure thickness of the model levels.
        """
        if self.level:
            self.levs_thick = level_thickness(self.level)
        else:
            self.levs_thick = None

    def set_grid_data(self):
        """Populate the attrs that hold grid data."""
        if self.grid_data_is_set:
            return
        self._set_mult_nc_grid_attr()
        self._set_sfc_area()
        self._set_levs_thick()
        self.grid_data_is_set = True
