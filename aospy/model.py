"""Functionality for representing data on disk of individual models."""
import logging

import numpy as np
import xarray as xr

from ._constants import RADIUS_EARTH
from . import internal_names
from . import utils


def _get_grid_attr(grid_objs, attr_name):
    """Get attribute from the grid_objs file(s)."""
    for xds in grid_objs:
        try:
            return getattr(xds, attr_name)
        except AttributeError:
            pass


def _rename_coords(ds, attrs):
    """Rename coordinates to aospy's internal names."""
    for name_int, names_ext in attrs.items():
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


def _bounds_from_array(arr, dim_name, bounds_name):
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


def _diff_bounds(bounds, coord):
    """Get grid spacing by subtracting upper and lower bounds."""
    try:
        return bounds[:, 1] - bounds[:, 0]
    except IndexError:
        diff = np.diff(bounds, axis=0)
        return xr.DataArray(diff, dims=coord.dims, coords=coord.coords)


def _grid_sfc_area(lon, lat, lon_bounds=None, lat_bounds=None):
    """Calculate surface area of each grid cell in a lon-lat grid."""
    # Compute the bounds if not given.
    if lon_bounds is None:
        lon_bounds = _bounds_from_array(
            lon, internal_names.LON_STR, internal_names.LON_BOUNDS_STR)
    if lat_bounds is None:
        lat_bounds = _bounds_from_array(
            lat, internal_names.LAT_STR, internal_names.LAT_BOUNDS_STR)
    # Compute the surface area.
    dlon = _diff_bounds(utils.vertcoord.to_radians(lon_bounds, is_delta=True),
                        lon)
    sinlat_bounds = np.sin(utils.vertcoord.to_radians(lat_bounds,
                                                      is_delta=True))
    dsinlat = np.abs(_diff_bounds(sinlat_bounds, lat))
    sfc_area = dlon*dsinlat*(RADIUS_EARTH**2)
    # Rename the coordinates such that they match the actual lat / lon.
    try:
        sfc_area = sfc_area.rename(
            {internal_names.LAT_BOUNDS_STR: internal_names.LAT_STR,
             internal_names.LON_BOUNDS_STR: internal_names.LON_STR})
    except ValueError:
        pass
    # Clean up: correct names and dimension order.
    sfc_area = sfc_area.rename(internal_names.SFC_AREA_STR)
    sfc_area[internal_names.LAT_STR] = lat
    sfc_area[internal_names.LON_STR] = lon
    return sfc_area.transpose()


class Model(object):
    """An object that describes a single climate or weather model.

    Each `Model` object is associated with a parent `Proj` object and also with
    one or more child `Run` objects.

    If aospy is being used to work with non climate- or weather-model data, the
    `Model` object can be used e.g. to represent a gridded observational
    product, with its child `Run` objects representing different released
    versions of that dataset.

    Attributes
    ----------
    name : str
        The model's name
    description : str
        A description of the model
    proj : {None, aospy.Proj}
        The model's parent aospy.Proj object
    runs : list
        A list of this model's child Run objects
    default_runs : list
        The default subset of child run objects on which to perform
        calculations via `aospy.Calc` with this model if not otherwise
        specified
    grid_file_paths : list
        The paths to netCDF files stored on disk from which the model's
        coordinate data can be taken.
    default_start_date, default_end_date : datetime.datetime
        The default start and end dates of any calculations using this Model
    """

    def __init__(self, name=None, description=None, proj=None,
                 grid_file_paths=None, default_start_date=None,
                 default_end_date=None, runs=None, default_runs=None,
                 load_grid_data=False, grid_attrs=None):
        """
        Parameters
        ----------
        name : str
            The model's name.  This must be unique from that of any other
            `Model` objects being used by the parent `Proj`.
        description : str, optional
            A description of the model.  This is not used internally by
            aospy; it is solely for the user's information.
        proj : {None, aospy.Proj}, optional
            The parent Proj object.  When the parent `Proj` object is
            instantiated with this Model included in its `models` attribute,
            this will be over-written with that `Proj` object.
        grid_file_paths : {None, sequence of strings}, optional
            The paths to netCDF files stored on disk from which the model's
            coordinate data can be taken.
        default_start_date : {None, `datetime.datetime`}, optional
            Default start date of calculations to be performed using
            this Model.
        default_end_date : {None, `datetime.datetime`}, optional
            Default end date of calculations to be performed using
            this Model.
        runs : {None, sequence of aospy.Run objects}, optional
            The child run objects of this Model
        default_runs : {None, sequence of aospy.Run objects}, optional
            The subset of this Model's runs over which to perform calculations
            by default.
        load_grid_data : bool, optional (default False)
            Whether or not to load the grid data specified by 'grid_file_paths'
            upon initilization
        grid_attrs : dict, optional (default None)
            Dictionary mapping aospy internal names of grid attributes to their
            corresponding names used in a particular model.
            E.g. ``{TIME_STR: 'T'}``.  While aospy checks for a number of
            alternative names for grid attributes used by various models,
            it is not possible to anticipate all possible names.  This option
            allows the user to explicitly tell aospy which variables correspond
            to which internal names (internal names not provided in this
            dictionary will be attempted to be found in the usual way).  For a
            list of built-in alternative names see
            :ref:`the table here <built-in-alternative-names>`.

        See Also
        --------
        aospy.DataLoader, aospy.Proj, aospy.Run

        Notes
        -----
        A side-effect of instantiating a Model object is that the `parent`
        attribute of all of the model's `Run` objects is set to that model.

        """
        if isinstance(name, str) and name:
            self.name = name
        else:
            raise ValueError("Non-empty string value of `name` is required")
        self.description = '' if description is None else description
        self.proj = proj

        grid_file_paths = [] if grid_file_paths is None else grid_file_paths
        self.grid_file_paths = grid_file_paths

        self.default_start_date = default_start_date
        self.default_end_date = default_end_date

        self.runs = runs
        [setattr(run, 'parent', self) for run in self.runs]

        if default_runs is None:
            self.default_runs = []
        else:
            self.default_runs = default_runs

        self.grid_attrs = grid_attrs

        self._grid_data_is_set = False
        if load_grid_data:
            self.set_grid_data()
            self._grid_data_is_set = True

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

    def _set_mult_grid_attr(self):
        """
        Set multiple attrs from grid file given their names in the grid file.
        """
        grid_objs = self._get_grid_files()

        if self.grid_attrs is None:
            self.grid_attrs = {}

        # Override GRID_ATTRS with entries in grid_attrs
        attrs = internal_names.GRID_ATTRS.copy()
        for k, v in self.grid_attrs.items():
            if k not in attrs:
                raise ValueError(
                    'Unrecognized internal name, {!r}, specified for a '
                    'custom grid attribute name.  See the full list of '
                    'valid internal names below:\n\n{}'.format(
                        k, list(internal_names.GRID_ATTRS.keys())))
            attrs[k] = (v, )

        for name_int, names_ext in attrs.items():
            for name in names_ext:
                grid_attr = _get_grid_attr(grid_objs, name)
                if grid_attr is not None:
                    TIME_STR = internal_names.TIME_STR
                    renamed_attr = _rename_coords(grid_attr, attrs)
                    if ((TIME_STR not in renamed_attr.dims) and
                       (TIME_STR in renamed_attr.coords)):
                        renamed_attr = renamed_attr.drop(TIME_STR)
                    setattr(self, name_int, renamed_attr)
                    break

    def set_grid_data(self):
        """Populate the attrs that hold grid data."""
        if self._grid_data_is_set:
            return
        self._set_mult_grid_attr()
        if not np.any(getattr(self, 'sfc_area', None)):
            try:
                sfc_area = _grid_sfc_area(self.lon, self.lat, self.lon_bounds,
                                          self.lat_bounds)
            except AttributeError:
                sfc_area = _grid_sfc_area(self.lon, self.lat)
            self.sfc_area = sfc_area
        try:
            self.levs_thick = utils.vertcoord.level_thickness(self.level)
        except AttributeError:
            self.level = None
            self.levs_thick = None
        self._grid_data_is_set = True
