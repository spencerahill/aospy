"""Functionality pertaining to aggregating data over geographical regions."""
import logging

from . import internal_names


def _add_to_mask(data, lat_bounds, lon_bounds):
    """Add mask spanning given lat-lon rectangle."""
    mask_lat = ((data[internal_names.LAT_STR] > lat_bounds[0]) &
                (data[internal_names.LAT_STR] < lat_bounds[1]))
    return mask_lat & ((data[internal_names.LON_STR] > lon_bounds[0]) &
                       (data[internal_names.LON_STR] < lon_bounds[1]))


def _make_mask(data, mask_bounds):
    """Construct the mask that defines this region."""
    # For each set of bounds add to the conditional.
    mask = False
    for lat_bounds, lon_bounds in mask_bounds:
        mask |= _add_to_mask(data, lat_bounds, lon_bounds)
    return mask


def _get_land_mask(data, do_land_mask):
    if not do_land_mask:
        return 1
    try:
        land_mask = data.land_mask.copy()
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
        logging.debug("Converting land mask from 0-100 to 0.0-1.0")
    except AttributeError:
        # Wrong for the edge case where no grid cell is 100% land.
        percent_bool = land_mask.max() == 100
    if percent_bool:
        land_mask *= 0.01
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


def _sum_over_lat_lon(arr):
    """Sum an array over the latitude and longitude dimensions."""
    return arr.sum(internal_names.LAT_STR).sum(internal_names.LON_STR)


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

    def __init__(self, name='', description='', lon_bounds=[], lat_bounds=[],
                 mask_bounds=[], do_land_mask=False):
        """Instantiate a Region object.

        If a region spans across the endpoint of the data's longitude array
        (i.e.  it crosses the Prime Meridian for data with longitudes spanning
        0 to 360), it must be defined as the union of two sections extending to
        the east and to the west of the Prime Meridian.

        Parameters
        ----------
        name : str
             The region's name.  This must be unique from that of any other
             `Region` objects being used by the overlying `Proj`.
        description : str, optional
             A description of the region.  This is not used internally by
             aospy; it is solely for the user's information.
        lon_bounds, lat_bounds : length-2 sequence, optional
             The longitude and latitude bounds of the region.  If the region
             boundaries are more complicated than a single lat-lon rectangle,
             use `mask_bounds` instead.
        mask_bounds : sequence, optional
             Each element is a length-2 tuple of the format `(lat_bounds,
             lon_bounds)`, where each of `lat_bounds` and `lon_bounds` are
             of the form described above.
        do_land_mask : { False, True, 'ocean', 'strict_land', 'strict_ocean'},
                       optional
             Determines what, if any, land mask is applied in addition to the
             mask defining the region's boundaries.  Default `False`.

             - True: apply the data's full land mask
             - False: apply no mask
             - 'ocean': mask out land rather than ocean
             - 'strict_land': mask out all points that are not 100% land
             - 'strict_ocean': mask out all points that are not 100% ocean

        Examples
        --------
        Define a region spanning the entire globe

        >>> globe = Region(name='globe', lat_bounds=(-90, 90),
        ...                lon_bounds=(0, 360), do_land_mask=False)

        Define a region corresponding to the Sahel region of Africa, which
        we'll define as land points within 10N-20N latitude and 18W-40E
        longitude.  Because this region crosses the 0 degrees longitude point,
        it has to be defined using `mask_bounds` as the union of two lat-lon
        rectangles.

        >>> sahel = Region(name='sahel', do_land_mask=True,
        ...                mask_bounds=[((10, 20), (0, 40)),
        ...                             ((10, 20), (342, 360))])

        """
        self.name = name
        self.description = description
        if lon_bounds and lat_bounds and not mask_bounds:
            self.mask_bounds = [(lat_bounds, lon_bounds)]
        else:
            self.mask_bounds = mask_bounds
        self.do_land_mask = do_land_mask

    def __str__(self):
        return 'Geographical region "' + self.name + '"'

    __repr__ = __str__

    def mask_var(self, data):
        """Mask the data of the given variable outside the region."""
        return data.where(_make_mask(data, self.mask_bounds))

    def ts(self, data):
        """Create time-series of region-average data."""
        data_masked = self.mask_var(data)
        sfc_area = data.sfc_area
        land_mask = _get_land_mask(data, self.do_land_mask)
        weights = _sum_over_lat_lon((self.mask_var(sfc_area)*land_mask))
        return (_sum_over_lat_lon(data_masked*sfc_area*land_mask) /
                weights)

    def av(self, data):
        """Time average of region-average time-series."""
        ts_ = self.ts(data)
        if 'year' not in ts_.coords:
            return ts_
        return ts_.mean('year')

    def std(self, data):
        """Standard deviation of region-average time-series."""
        ts_ = self.ts(data)
        if 'year' not in ts_.coords:
            return ts_
        return ts_.std('year')
