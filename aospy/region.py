"""Functionality pertaining to aggregating data over geographical regions."""
from collections import namedtuple
import logging

import numpy as np

from .internal_names import (
    LAND_MASK_STR,
    LAT_STR,
    LON_STR,
    SFC_AREA_STR,
    YEAR_STR
)
from .utils.longitude import _maybe_cast_to_lon


def _get_land_mask(data, do_land_mask, land_mask_str=LAND_MASK_STR):
    if not do_land_mask:
        return 1
    try:
        land_mask = data[land_mask_str].copy()
    except AttributeError:
        # TODO: Implement aospy built-in land mask to default to.
        msg = ("No land mask found.  Using empty mask, which amounts to "
               "no land or ocean mask being applied.  Regions that use a "
               "land or ocean mask will therefore NOT be accurately "
               "computed.")
        logging.warning(msg)
        return 1
    try:
        percent_bool = land_mask.units.lower() in ('%', 'percent')
    except AttributeError:
        percent_bool = np.any(land_mask > 1)
    if percent_bool:
        land_mask *= 0.01
        logging.debug("Converting land mask from 0-100 to 0.0-1.0")
    if do_land_mask is True:
        return land_mask
    if do_land_mask == 'ocean':
        return 1. - land_mask
    if do_land_mask in ('strict_land', 'strict_ocean'):
        raise NotImplementedError
    msg = ("'do_land_mask' value of '{0}' is not one of the valid "
           "choices: [True, False, 'ocean', 'strict_land', "
           "'strict_ocean']").format(do_land_mask)
    raise ValueError(msg)


class BoundsRect(namedtuple('BoundsRect', ['west', 'east', 'south', 'north'])):
    """Bounding longitudes and latitudes of a given lat-lon rectangle."""
    def __new__(cls, west, east, south, north):
        new_west = _maybe_cast_to_lon(west, strict=True)
        new_east = _maybe_cast_to_lon(east, strict=True)
        return super(BoundsRect, cls).__new__(cls, new_west, new_east,
                                              south, north)

    def __repr__(self):
        return ("BoundsRect(west={0}, east={1}, south={2}, "
                "north={3}".format(self.west, self.east, self.south,
                                   self.north))


class Region(object):
    """Geographical region over which to perform averages and other reductions.

    Each `Proj` object includes a list of `Region` objects, which is used by
    `Calc` to determine which regions over which to perform time reductions
    over region-average quantities.

    Region boundaries are specified as either a single "rectangle" in latitude
    and longitude or the union of multiple such rectangles.  In addition, a
    land or ocean mask can be applied.

    Attributes
    ----------
    name : str
        The region's name
    description : str
        A description of the region
    mask_bounds : tuple
        The coordinates definining the lat-lon rectangle(s) that define(s)
        the region's boundaries
    do_land_mask
        Whether values occurring over land, ocean, or neither are excluded from
        the region, and whether the included points must be strictly land or
        ocean or if fractional land/ocean values are included.

    See Also
    --------
    aospy.Calc.region_calcs

    """

    def __init__(self, name='', description='', west_bound=None,
                 east_bound=None, south_bound=None, north_bound=None,
                 mask_bounds=None, do_land_mask=False):
        """Instantiate a Region object.

        Note that longitudes spanning (-180, 180), (0, 360), or any other range
        are all supported: -180 to 0, 180 to 360, etc. are interpreted as the
        western hemisphere, and 0-180, 360-540, etc. are interpreted as the
        eastern hemisphere.  This is true both of the region definition and of
        any data upon which the region mask is applied.

        E.g. suppose some of your data is defined on a -180 to 180 longitude
        grid, some of it is defined on a 0 to 360 grid, and some of it is
        defined on a -70 to 290 grid.  A single Region object will work with
        all three of these.

        Conversely, latitudes are always treated as -90 as the South Pole, 0 as
        the Equator, and 90 as the North Pole.  Latitudes larger than 90 are
        not physically meaningful.

        Parameters
        ----------
        name : str
             The region's name.  This must be unique from that of any other
             `Region` objects being used by the overlying `Proj`.
        description : str, optional
             A description of the region.  This is not used internally by
             aospy; it is solely for the user's information.
        west_bound, east_bound : { scalar, aospy.Longitude }, optional
             The western and eastern boundaries of the region.  All input
             longitudes are casted to ``aospy.Longitude`` objects, which
             essentially maps them to a 180W to 180E grid.  The region's
             longitudes always start at ``west_bound`` and move toward the east
             until reaching ``east_bound``.  This means that there are two
             distinct cases:

             - If, after this translation, ``west_bound`` is less than
               ``east_bound``, the region includes the points east of
               ``west_bound`` and west of ``east_bound``.
             - If ``west_bound`` is greater than ``east_bound``, then the
               region is treated as wrapping around the dateline, i.e. it's
               western-most point is ``east_bound``, and it includes all points
               moving east from there until ``west_bound``.

             If the region boundaries are more complicated than a single
             lat-lon rectangle, use ``mask_bounds`` instead.

        south_bound, north_bound : scalar, optional
             The southern, and northern boundaries, respectively, of the
             region.  If the region boundaries are more complicated than a
             single lat-lon rectangle, use `mask_bounds` instead.
        mask_bounds : sequence, optional
             Each element is a length-4 sequence of the format `(west_bound,
             east_bound, south_bound, north_bound)`, where each of these
             `_bound` arguments is of the form described above.
        do_land_mask, bool or str, optional
             Determines what, if any, land mask is applied in addition to the
             mask defining the region's boundaries.  Default `False`.  Must be
             one of ``False``, ``True``, 'ocean', 'strict_land', or
             'strict_ocean':

             - ``True``: apply the data's full land mask
             - ``False``: apply no mask
             - 'ocean': mask out land rather than ocean
             - 'strict_land': mask out all points that are not 100% land
             - 'strict_ocean': mask out all points that are not 100% ocean

        Examples
        --------
        Define a region spanning the entire globe:

        >>> globe = Region(name='globe', west_bound=0, east_bound=360,
        ...                south_bound=-90, north_bound=90, do_land_mask=False)

        Longitudes are handled as cyclic, so this definition could have
        equivalently used `west_bound=-180, east_bound=180` or `west_bound=200,
        east_bound=560`, or anything else that spanned 360 degrees total.

        Define a region corresponding to land in the mid-latitudes, which we'll
        define as land points within 30-60 degrees latitude in both
        hemispheres.  Because this region is the union of multiple lat-lon
        rectangles, it has to be defined using `mask_bounds`:

        >>> land_mid_lats = Region(name='land_midlat', do_land_mask=True,
        ...                        mask_bounds=[(-180, 180, 30, 60),
        ...                                     (-180, 180, -30, -60)])

        Define a region spanning the southern Tropical Atlantic ocean, which
        we'll take to be all ocean points between 60W and 30E and between the
        Equator and 30S:

        >>> atl_south_trop = Region(name='atl_sh_trop', west_bound=-60,
        ...                         east_bound=30, south_bound=-30,
        ...                         north_bound=0, do_land_mask='ocean')

        Define the "opposite" region, i.e. all ocean points in the southern
        Tropics *outside* of the Atlantic.  We simply swap ``west_bound`` and
        ``east_bound`` of the previous example:

        >>> non_atl_south_trop = Region(name='non_atl_sh_trop', west_bound=30,
        ...                             east_bound=-60, south_bound=-30,
        ...                             north_bound=0, do_land_mask='ocean')

        """
        self.name = name
        self.description = description
        self.do_land_mask = do_land_mask

        if mask_bounds is None:
            self.mask_bounds = tuple([BoundsRect(west_bound, east_bound,
                                                 south_bound, north_bound)])
        else:
            bounds = []
            for rect_bounds in mask_bounds:
                if len(rect_bounds) != 4:
                    raise ValueError("Each element of `mask_bounds` must be a "
                                     "length-4 array with values (west, east, "
                                     "south, north).  Value given: "
                                     "{}".format(rect_bounds))
                else:
                    bounds.append(BoundsRect(*rect_bounds))
            self.mask_bounds = tuple(bounds)

    def __str__(self):
        return 'Geographical region "' + self.name + '"'

    def _make_mask(self, data, lon_str=LON_STR, lat_str=LAT_STR):
        """Construct the mask that defines a region on a given data's grid."""
        mask = False
        for west, east, south, north in self.mask_bounds:
            if west < east:
                mask_lon = (data[lon_str] > west) & (data[lon_str] < east)
            else:
                mask_lon = (data[lon_str] < west) | (data[lon_str] > east)
            mask_lat = (data[lat_str] > south) & (data[lat_str] < north)
            mask |= mask_lon & mask_lat
        return mask

    def mask_var(self, data, lon_cyclic=True, lon_str=LON_STR,
                 lat_str=LAT_STR):
        """Mask the given data outside this region.

        Parameters
        ----------
        data : xarray.DataArray
            The array to be regionally masked.
        lon_cyclic : bool, optional (default True)
            Whether or not the longitudes of ``data`` span the whole globe,
            meaning that they should be wrapped around as necessary to cover
            the Region's full width.
        lon_str, lat_str : str, optional
            The names of the longitude and latitude dimensions, respectively,
            in the data to be masked.  Defaults are
            ``aospy.internal_names.LON_STR`` and
            ``aospy.internal_names.LON_STR``, respectively.

        Returns
        -------
        xarray.DataArray
            The original array with points outside of the region masked.

        """
        # TODO: is this still necessary?
        if not lon_cyclic:
            if self.west_bound > self.east_bound:
                raise ValueError("Longitudes of data to be masked are "
                                 "specified as non-cyclic, but Region's "
                                 "definition requires wraparound longitudes.")
        masked = data.where(self._make_mask(data, lon_str=lon_str,
                                            lat_str=lat_str))
        return masked

    def ts(self, data, lon_cyclic=True, lon_str=LON_STR, lat_str=LAT_STR,
           land_mask_str=LAND_MASK_STR, sfc_area_str=SFC_AREA_STR):
        """Create yearly time-series of region-averaged data.

        Parameters
        ----------
        data : xarray.DataArray
            The array to create the regional timeseries of
        lon_cyclic : { None, True, False }, optional (default True)
            Whether or not the longitudes of ``data`` span the whole globe,
            meaning that they should be wrapped around as necessary to cover
            the Region's full width.
        lat_str, lon_str, land_mask_str, sfc_area_str : str, optional
            The name of the latitude, longitude, land mask, and surface area
            coordinates, respectively, in ``data``.  Defaults are the
            corresponding values in ``aospy.internal_names``.

        Returns
        -------
        xarray.DataArray
            The timeseries of values averaged within the region and within each
            year, one value per year.

        """
        data_masked = self.mask_var(data, lon_cyclic=lon_cyclic,
                                    lon_str=lon_str, lat_str=lat_str)
        sfc_area = data[sfc_area_str]
        sfc_area_masked = self.mask_var(sfc_area, lon_cyclic=lon_cyclic,
                                        lon_str=lon_str, lat_str=lat_str)
        land_mask = _get_land_mask(data, self.do_land_mask,
                                   land_mask_str=land_mask_str)
        weights = sfc_area_masked * land_mask
        # Mask weights where data values are initially invalid in addition
        # to applying the region mask.
        weights = weights.where(np.isfinite(data))
        weights_reg_sum = weights.sum(lon_str).sum(lat_str)
        data_reg_sum = (data_masked * sfc_area_masked *
                        land_mask).sum(lat_str).sum(lon_str)
        return data_reg_sum / weights_reg_sum

    def av(self, data, lon_str=LON_STR, lat_str=LAT_STR,
           land_mask_str=LAND_MASK_STR, sfc_area_str=SFC_AREA_STR):
        """Time-average of region-averaged data.

        Parameters
        ----------
        data : xarray.DataArray
            The array to compute the regional time-average of
        lat_str, lon_str, land_mask_str, sfc_area_str : str, optional
            The name of the latitude, longitude, land mask, and surface area
            coordinates, respectively, in ``data``.  Defaults are the
            corresponding values in ``aospy.internal_names``.

        Returns
        -------
        xarray.DataArray
            The region-averaged and time-averaged data.

        """
        ts = self.ts(data, lon_str=lon_str, lat_str=lat_str,
                     land_mask_str=land_mask_str, sfc_area_str=sfc_area_str)
        if YEAR_STR not in ts.coords:
            return ts
        else:
            return ts.mean(YEAR_STR)

    def std(self, data, lon_str=LON_STR, lat_str=LAT_STR,
            land_mask_str=LAND_MASK_STR, sfc_area_str=SFC_AREA_STR):
        """Temporal standard deviation of region-averaged data.

        Parameters
        ----------
        data : xarray.DataArray
            The array to compute the regional time-average of
        lat_str, lon_str, land_mask_str, sfc_area_str : str, optional
            The name of the latitude, longitude, land mask, and surface area
            coordinates, respectively, in ``data``.  Defaults are the
            corresponding values in ``aospy.internal_names``.

        Returns
        -------
        xarray.DataArray
            The temporal standard deviation of the region-averaged data

        """
        ts = self.ts(data, lon_str=lon_str, lat_str=lat_str,
                     land_mask_str=land_mask_str, sfc_area_str=sfc_area_str)
        if YEAR_STR not in ts.coords:
            return ts
        else:
            return ts.std(YEAR_STR)
