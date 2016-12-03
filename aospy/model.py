"""model.py: Model class of aospy for storing attributes of a GCM."""
from __future__ import print_function
import glob
import logging
import os

import numpy as np
import xarray as xr

from .__config__ import (LAT_STR, LON_STR, PHALF_STR, PFULL_STR, PLEVEL_STR,
                         LAT_BOUNDS_STR, LON_BOUNDS_STR, TIME_BOUNDS_STR,
                         BOUNDS_STR)
from .constants import r_e
from .utils.times import datetime_or_default
from . import utils


class Model(object):
    """Parameters of local data associated with a single climate model."""
    def __init__(self, name='', description='', proj=None,
                 grid_file_paths=None, data_direc=None,
                 data_dir_struc=None, data_dur=None,
                 data_start_date=None, data_end_date=None,
                 default_start_date=None, default_end_date=None,
                 runs={}, default_runs={},
                 load_grid_data=False, repo_version=None,
                 repo_ens_mem='r1i1p1'):
        self.name = name
        self.description = description
        self.proj = proj

        grid_file_paths = [] if grid_file_paths is None else grid_file_paths
        self.grid_file_paths = grid_file_paths
        self.repo_version = repo_version
        self.repo_ens_mem = repo_ens_mem

        self.data_direc = data_direc
        self.data_dir_struc = data_dir_struc
        self.data_dur = data_dur

        self.data_start_date = datetime_or_default(data_start_date, None)
        self.data_end_date = datetime_or_default(data_end_date, None)
        self.default_start_date = datetime_or_default(default_start_date,
                                                      self.data_start_date)
        self.default_end_date = datetime_or_default(default_end_date,
                                                    self.data_end_date)

        self.runs = utils.io.dict_name_keys(runs)
        [setattr(run, 'parent', self) for run in self.runs.values()]
        if default_runs:
            self.default_runs = utils.io.dict_name_keys(default_runs)
        else:
            self.default_runs = {}

        self.grid_data_is_set = False
        if load_grid_data:
            self.set_grid_data()
            self.grid_data_is_set = True

    def __str__(self):
        return 'Model instance "' + self.name + '"'

    __repr__ = __str__

    def find_data_direc_repo(self, run_name='amip', var_name='ta',
                             direc_sub='mon/atmos/Amon', ens_mem=None):
        """Find the netCDF files used to populate grid attrs for a GFDL repo"""
        if ens_mem is None:
            ens_mem = self.repo_ens_mem
        direc = os.path.join(self.data_direc, run_name, direc_sub,
                             ens_mem)
        if not os.path.isdir(direc):
            direc = os.path.join(self.data_direc, run_name, direc_sub,
                                 'r1i1p1')
        direc_full = utils.io.get_data_direc_repo(direc, var_name,
                                            version=self.repo_version)
        files = glob.glob(os.path.join(direc_full, var_name + '_*.nc'))
        # Some repos have all variables in the same directory.
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

    def _get_grid_files_repo(self):
        """Get grid files of data in CMIP5 archive."""
        direc_sub = 'fx/atmos/fx'
        ens_mem = 'r0i0p0'
        grid = []
        for var in ('sftlf', 'areacella', 'orog'):
            try:
                path = self.find_data_direc_repo(var_name=var,
                                                    direc_sub=direc_sub,
                                                    ens_mem=ens_mem)
            except IOError as e:
                logging.debug(str(repr(e)))
            else:
                grid.append(path)
        non_grid = [self.find_data_direc_repo()]
        return grid + non_grid

    def _get_grid_files(self):
        """Get the files holding grid data for an aospy object."""
        if getattr(self.proj, 'name', None) == 'cmip5':
            try:
                grid_file_paths = self._get_grid_files_repo()
            except OSError:
                grid_file_paths = []
            if self.grid_file_paths:
                grid_file_paths += self.grid_file_paths
        else:
            grid_file_paths = self.grid_file_paths
        datasets = []
        if isinstance(grid_file_paths, str):
            grid_file_paths = [grid_file_paths]
        for path in grid_file_paths:
            try:
                ds = xr.open_dataset(path, decode_times=False)
            except TypeError:
                ds = xr.open_mfdataset(path, decode_times=False)
            except (RuntimeError, OSError) as e:
                msg = str(e) + ': {}'.format(path)
                raise RuntimeError(msg)
            datasets.append(ds)
        return tuple(datasets)

    @staticmethod
    def _get_grid_attr(grid_objs, attr_name):
        """Get attribute from the grid_objs file(s)."""
        for xds in grid_objs:
            try:
                return getattr(xds, attr_name)
            except AttributeError:
                pass

    @staticmethod
    def _rename_coords(ds):
        """
        Renames all coordinates within a Dataset or DataArray so that they
        match the internal names.
        """
        primary_attrs = {
            LAT_STR:         ('lat', 'latitude', 'LATITUDE', 'y', 'yto'),
            LAT_BOUNDS_STR:  ('latb', 'lat_bnds', 'lat_bounds'),
            LON_STR:         ('lon', 'longitude', 'LONGITUDE', 'x', 'xto'),
            LON_BOUNDS_STR:  ('lonb', 'lon_bnds', 'lon_bounds'),
            PLEVEL_STR:      ('level', 'lev', 'plev'),
            PHALF_STR:       ('phalf',),
            PFULL_STR:       ('pfull',),
            TIME_BOUNDS_STR: ('time_bounds', 'time_bnds'),
            }
        if isinstance(ds, (xr.DataArray, xr.Dataset)):
            for name_int, names_ext in primary_attrs.items():
                # Check if coord is in dataset already.
                ds_coord_name = set(names_ext).intersection(set(ds.coords))
                if ds_coord_name:
                    # Rename to the aospy internal name.
                    try:
                        ds = ds.rename({list(ds_coord_name)[0]: name_int})
                        logging.debug("Rename coord from `{0}` to `{1}` for "
                                      "Dataset `{2}`".format(ds_coord_name,
                                                             name_int, ds))
                    # xarray throws a ValueError if the name already exists
                    except ValueError:
                        ds = ds
        return ds

    def _set_mult_grid_attr(self):
        """
        Set multiple attrs from grid file given their names in the grid file.
        """
        grid_attrs = {
            LAT_STR:         ('lat', 'latitude', 'LATITUDE', 'y', 'yto'),
            LAT_BOUNDS_STR:    ('latb', 'lat_bnds', 'lat_bounds'),
            LON_STR:         ('lon', 'longitude', 'LONGITUDE', 'x', 'xto'),
            LON_BOUNDS_STR:    ('lonb', 'lon_bnds', 'lon_bounds'),
            PLEVEL_STR:      ('level', 'lev', 'plev'),
            # 'time':        ('time', 'TIME'),
            # 'time_st':     ('average_T1',),
            # 'time_end':    ('average_T2',),
            'time_dur':    ('average_DT',),
            'time_bounds': ('time_bounds', 'time_bnds'),
            # 'year':        ('yr', 'year'),
            # 'month':       ('mo', 'month'),
            # 'day':         ('dy', 'day'),
            'sfc_area':    ('area', 'sfc_area', 'areacella'),
            'zsurf':       ('zsurf', 'orog'),
            'land_mask':   ('land_mask', 'sftlf'),
            'pk':          ('pk',),
            'bk':          ('bk',),
            PHALF_STR:     ('phalf',),
            PFULL_STR:     ('pfull',),
            BOUNDS_STR:    ('bnds', 'bounds'),
            # 'nrecords':    ('nrecords',),
            # 'idim':        ('idim',),
            # 'fill_value':  ('fill_value')
        }
        grid_objs = self._get_grid_files()
        for name_int, names_ext in grid_attrs.items():
            for name in names_ext:
                grid_attr = self._get_grid_attr(grid_objs, name)
                if grid_attr is not None:
                    # Rename coordinates to aospy's internal names.
                    renamed_attr = self._rename_coords(grid_attr)
                    # Hack: MERRA data ending up with time-defined lat- and
                    # lon-bounds.  Drop the time term.
                    # TODO: improve this hack.
                    cond_to_drop = (
                        hasattr(renamed_attr, 'time') and
                        name_int not in ['time', 'time_st', 'time_end',
                                         'time_dur', 'time_bounds']
                    )
                    if cond_to_drop:
                        tmp = renamed_attr.isel(time=0)
                        try:
                            renamed_attr = tmp.squeeze('time').drop('time')
                        except KeyError:
                            renamed_attr = tmp.drop('time')
                    setattr(self, name_int, renamed_attr)
                    break

    @staticmethod
    def bounds_from_array(arr, dim_name, bounds_name):
        """Get the bounds of an array given its center values.

        E.g. if lat-lon grid center lat/lon values are known, but not the
        bounds of each grid box.  The algorithm assumes that the bounds
        are simply halfway between each pair of center values.
        """
        # TODO: don't assume needed dimension is in axis=0
        # TODO: refactor to get rid of repetitive code
        spacing = arr.diff(dim_name).values
        lower = xr.DataArray(np.empty_like(arr), dims=arr.dims,
                             coords=arr.coords)
        lower.values[:-1] = arr.values[:-1] - 0.5*spacing
        lower.values[-1] = arr.values[-1] - 0.5*spacing[-1]
        upper = xr.DataArray(np.empty_like(arr), dims=arr.dims,
                             coords=arr.coords)
        upper.values[:-1] = arr.values[:-1] + 0.5*spacing
        upper.values[-1] = arr.values[-1] + 0.5*spacing[-1]
        bounds = xr.concat([lower, upper], dim='bounds')
        return bounds.T

    @staticmethod
    def diff_bounds(bounds, coord):
        """Get grid spacing by subtracting upper and lower bounds."""
        try:
            return bounds[:, 1] - bounds[:, 0]
        except IndexError:
            diff = np.diff(bounds, axis=0)
            return xr.DataArray(diff, dims=coord.dims, coords=coord.coords)

    @classmethod
    def grid_sfc_area(cls, lon, lat, lon_bounds=None, lat_bounds=None):
        """Calculate surface area of each grid cell in a lon-lat grid."""
        # Compute the bounds if not given.
        if lon_bounds is None:
            lon_bounds = cls.bounds_from_array(lon, LON_STR, LON_BOUNDS_STR)
        if lat_bounds is None:
            lat_bounds = cls.bounds_from_array(lat, LAT_STR, LAT_BOUNDS_STR)
        # Compute the surface area.
        dlon = cls.diff_bounds(utils.vertcoord.to_radians(lon_bounds,
                                                          is_delta=True), lon)
        sinlat_bounds = np.sin(utils.vertcoord.to_radians(lat_bounds,
                                                          is_delta=True))
        dsinlat = np.abs(cls.diff_bounds(sinlat_bounds, lat))
        sfc_area = dlon*dsinlat*(r_e**2)
        # Rename the coordinates such that they match the actual lat / lon.
        try:
            sfc_area = sfc_area.rename({LAT_BOUNDS_STR: LAT_STR,
                                        LON_BOUNDS_STR: LON_STR})
        except ValueError:
            pass
        # Clean up: correct names and dimension order.
        sfc_area = sfc_area.rename('sfc_area')
        sfc_area[LAT_STR] = lat
        sfc_area[LON_STR] = lon
        return sfc_area.transpose()

    def set_grid_data(self):
        """Populate the attrs that hold grid data."""
        if self.grid_data_is_set:
            return
        self._set_mult_grid_attr()
        if not np.any(getattr(self, 'sfc_area', None)):
            try:
                sfc_area = self.grid_sfc_area(self.lon, self.lat,
                                              self.lon_bounds, self.lat_bounds)
            except AttributeError:
                sfc_area = self.grid_sfc_area(self.lon, self.lat)
            self.sfc_area = sfc_area
        try:
            self.levs_thick = utils.vertcoord.level_thickness(self.level)
        except AttributeError:
            self.level = None
            self.levs_thick = None
        self.grid_data_is_set = True
