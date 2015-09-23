"""model.py: Model class of aospy for storing attributes of a GCM."""
import glob
import os

import numpy as np
import netCDF4
import xray

from .constants import r_e
from .io import get_nc_direc_repo
from .utils import dict_name_keys, level_thickness


class Model(object):
    """Parameters of local data associated with a single climate model."""
    def __init__(self, name='', description='', proj=False, nc_grid_paths=(),
                 nc_dur=False, nc_start_yr=False, nc_end_yr=False,
                 nc_start_month=False, nc_end_month=False,
                 default_yr_range=False, runs={}, default_runs={},
                 load_grid_data=False, nc_dir_struc=False, nc_direc=False,
                 repo_version=-1, read_mode='netcdf4'):
        self.name = name
        self.description = description
        self.proj = proj
        
        self.read_mode = read_mode

        self.nc_direc = nc_direc
        if nc_dir_struc == 'gfdl_repo':
            self.repo_version = repo_version
            self.nc_grid_paths = self.find_nc_direc_repo()
        else:
            self.nc_grid_paths = nc_grid_paths
        self.nc_dir_struc = nc_dir_struc

        self.nc_start_yr = nc_start_yr
        self.nc_end_yr = nc_end_yr
        self.nc_start_month = nc_start_month
        self.nc_end_month = nc_end_month
        self.nc_dur = nc_dur


        self.default_yr_range = default_yr_range
        self.runs = dict_name_keys(runs)
        [setattr(run, 'parent', self) for run in self.runs.values()]
        if default_runs:
            self.default_runs = dict_name_keys(default_runs)
        else:
            self.default_runs = {}

        self.grid_data_is_set = False
        if load_grid_data:
            self.set_grid_data()
            self.grid_data_is_set = True

    def __str__(self):
        return 'Model instance "' + self.name + '"'

    __repr__ = __str__

    def find_nc_direc_repo(self, run_name='amip', var_name='ta',
                           direc_sub='mon/atmos/Amon/r1i1p1'):
        """Find the netCDF files used to populate grid attrs for a GFDL repo"""
        direc = os.path.join(self.nc_direc, run_name, direc_sub)
        direc_full = get_nc_direc_repo(direc, var_name,
                                       version=self.repo_version)
        # Some repos have all variables in the same directory.
        files = glob.glob(os.path.join(direc_full, var_name + '_*.nc'))
        if len(files) == 0:
            files = glob.glob(os.path.join(direc_full, var_name,
                                           var_name + '_*.nc'))
        # Others have subdirectories for each variable.
        if len(files) > 0:
            files.sort()
            return files
        else:
            raise IOError("Specified files don't exist for var name '%s' "
                          "in directory '%s" % (var_name, direc_full))

        return glob.glob(os.path.join(direc_full, var_name + '_*.nc'))


    def _get_nc_grid(self):
        """Get the nc_grid of an aospy object."""
        nc_objs = []
        for path in self.nc_grid_paths:
            try:
                if self.read_mode == 'netcdf4':
                    nc_obj = netCDF4.MFDataset(path)
                elif self.read_mode == 'xray':
                    # I'm not sure why we read in the `time` array from the model grid.
                    # I don't think there is any harm in not decoding the CF conventions
                    # in other words not parsing the dates. Trying to decode CF can trigger
                    # a warning (that I don't think is necessary here), so we don't bother.
                    nc_obj = xray.open_dataset(path, 
                                               drop_variables=['time_bounds', 'nv'],
                                               decode_cf=False)
                else:
                    pass
            except RuntimeError:
                nc_obj = netCDF4.MFDataset([path])
            nc_objs.append(nc_obj)
        return tuple(nc_objs)

    def _get_nc_grid_attr(self, nc_grid, attr_name):
        """Get attribute from the nc_grid file(s)."""
        for nc in nc_grid:
            try:
                return nc.variables[attr_name][:]
            except KeyError:
                pass
        raise KeyError

    def _get_nc_grid_attr_xray(self, nc_grid, attr_name):
        for nc in nc_grid:
            try:
                # Return a copy of the DataArray, because otherwise you won't catch 
                # the exception (it needs to look for `attr_name` here; otherwise
                # it is lazy and waits till it's needed later)
                return nc[attr_name].copy(deep=True)
            except RuntimeError:
                pass
        raise RuntimeError

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
            for name_int, names_ext in grid_attrs.items():
                for name in names_ext:
                    try:
                        if self.read_mode == 'netcdf4':
                            setattr(self, name_int,
                                self._get_nc_grid_attr(nc_grid, name))
                        elif self.read_mode == 'xray':
                            setattr(self, name_int, 
                                    self._get_nc_grid_attr_xray(nc_grid, name))
                    except (KeyError, RuntimeError):
                        pass
                    else:
                        break
                else:
                    setattr(self, name_int, None)
        except:
            raise
        finally:
            for nc in nc_grid:
                nc.close()

    def grid_sfc_area(self, lon, lat, lonb, latb, gridcenter=False):
        """Calculate surface area of each grid cell in a lon-lat grid."""

        if self.read_mode == 'netcdf4':
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

        elif self.read_mode == 'xray':
            # The odd business of converting to radians then back to degrees
            # is to account for a bug in xray. For whatever reason you cannot
            # call the diff function on a coordinate. If we apply an identity
            # operation (e.g. changing back and forth between degrees and radians)
            # we can disconnect the coordinate from the values. This allows diff
            # to work properly.
            dlon = np.abs(np.rad2deg(np.deg2rad(lonb)).diff(dim='lonb'))
            dsinlat = np.abs(np.sin(np.pi * latb / 180.0).diff(dim='latb'))
            
            # We will need to be clever about naming, however. For now we don't
            # use the gridbox area, so we'll leave things like this for now.
            sfc_area = dlon*dsinlat*(r_e**2) * (np.pi/180.)
            # Rename the coordinates such that they match the actual lat / lon 
            # (Not the bounds)
            sfc_area = sfc_area.rename({'latb' : 'lat', 'lonb' : 'lon'})
            sfc_area['lat'] = lat
            sfc_area['lon'] = lon
            return sfc_area
        else:
            pass
        
    def _set_sfc_area(self):
        """Set the 2D array holding the surface area of gridboxes."""
        if getattr(self, 'sfc_area', None) is not None:
            return
        else:
            try:
                sfc_area = self.grid_sfc_area(
                    self.lon, self.lat,
                    self.lon_bounds, self.lat_bounds, gridcenter=False
                )
            except AttributeError:
                sfc_area = self.grid_sfc_area(
                    self.lon, self.lat, gridcenter=True
                )
            self.sfc_area = sfc_area

    def set_grid_data(self):
        """Populate the attrs that hold grid data."""
        if self.grid_data_is_set:
            return
        self._set_mult_nc_grid_attr()
        self._set_sfc_area()
        if self.level is not None:
            self.levs_thick = level_thickness(self.level)
        else:
            self.levs_thick = None
        self.grid_data_is_set = True
