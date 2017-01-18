"""model.py: Model class of aospy for storing attributes of a GCM."""
from __future__ import print_function
import logging

import dask
import numpy as np
import xarray as xr

from .constants import r_e
from . import internal_names
from .utils.times import datetime_or_default
from . import utils

# Dask must use its serial scheduler if computations are to be performed
# in parallel using multiprocess
dask.set_options(get=dask.async.get_sync)


class Model(object):
    """Parameters of local data associated with a single climate model."""
    def __init__(self, name='', description='', proj=None,
                 grid_file_paths=None, default_start_date=None,
                 default_end_date=None,
                 runs={}, default_runs={},
                 load_grid_data=False,
                 data_loader=None):
        self.name = name
        self.description = description
        self.proj = proj

        grid_file_paths = [] if grid_file_paths is None else grid_file_paths
        self.grid_file_paths = grid_file_paths

        self.default_start_date = datetime_or_default(
            default_start_date, getattr(data_loader, 'data_start_date', None))
        self.default_end_date = datetime_or_default(
            default_end_date, getattr(data_loader, 'data_end_date', None))

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

        self.data_loader = data_loader

    def __str__(self):
        return 'Model instance "' + self.name + '"'

    __repr__ = __str__

    def _get_grid_files(self):
        """Get the files holding grid data for an aospy object."""
        grid_file_paths = self.grid_file_paths
        datasets = []
        if isinstance(grid_file_paths, str):
            grid_file_paths = [grid_file_paths]
        for path in grid_file_paths:
            try:
                ds = xr.open_dataset(path, decode_times=False)
            except TypeError:
                ds = xr.open_mfdataset(path, decode_times=False).load()
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
        if isinstance(ds, (xr.DataArray, xr.Dataset)):
            for name_int, names_ext in internal_names.GRID_ATTRS.items():
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
        grid_objs = self._get_grid_files()
        for name_int, names_ext in internal_names.GRID_ATTRS.items():
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
            lon_bounds = cls.bounds_from_array(
                lon, internal_names.LON_STR, internal_names.LON_BOUNDS_STR)
        if lat_bounds is None:
            lat_bounds = cls.bounds_from_array(
                lat, internal_names.LAT_STR, internal_names.LAT_BOUNDS_STR)
        # Compute the surface area.
        dlon = cls.diff_bounds(utils.vertcoord.to_radians(lon_bounds,
                                                          is_delta=True), lon)
        sinlat_bounds = np.sin(utils.vertcoord.to_radians(lat_bounds,
                                                          is_delta=True))
        dsinlat = np.abs(cls.diff_bounds(sinlat_bounds, lat))
        sfc_area = dlon*dsinlat*(r_e**2)
        # Rename the coordinates such that they match the actual lat / lon.
        try:
            sfc_area = sfc_area.rename(
                {internal_names.LAT_BOUNDS_STR: internal_names.LAT_STR,
                 internal_names.LON_BOUNDS_STR: internal_names.LON_STR})
        except ValueError:
            pass
        # Clean up: correct names and dimension order.
        sfc_area = sfc_area.rename('sfc_area')
        sfc_area[internal_names.LAT_STR] = lat
        sfc_area[internal_names.LON_STR] = lon
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
