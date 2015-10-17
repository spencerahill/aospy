"""region.py: Region class and region_inst()."""
import xray


class Region(object):
    """Geographical region."""
    def __init__(self, name='', description='', lon_bounds=[], lat_bounds=[],
                 mask_bounds=[], land_mask=False):
        """Instantiate a Region object."""
        self.name = name
        self.description = description
        if lon_bounds and lat_bounds and not mask_bounds:
            self.mask_bounds = [(lat_bounds, lon_bounds)]
        else:
            self.mask_bounds = mask_bounds
        self.land_mask = land_mask

    def __str__(self):
        return 'Geographical region "' + self.name + '"'

    __repr__ = __str__

    def _add_to_mask(self, data, lat_bounds, lon_bounds):
        """Add mask spanning given lat-lon rectangle."""
        mask_lat = ((data['lat'] > lat_bounds[0]) &
                    (data['lat'] < lat_bounds[1]))
        return mask_lat & ((data['lon'] > lon_bounds[0]) &
                           (data['lon'] < lon_bounds[1]))

    def make_mask(self, data):
        """Construct the mask that defines this region."""
        # For each set of bounds add to the conditional.
        mask = False
        for lat_bounds, lon_bounds in self.mask_bounds:
            mask |= self._add_to_mask(data, lat_bounds, lon_bounds)
        return mask

    def mask_var(self, data):
        """Mask the data of the given variable outside the region."""
        return data.where(self.make_mask(data))

    def _get_sfc_area(self, model, dims, coords):
        try:
            return xray.DataArray(model.sfc_area, dims=dims, coords=coords)
        except:
            raise

    def _get_land_mask(self, model, dims, coords):
        try:
            lmask = xray.DataArray(model.land_mask, dims=dims, coords=coords)
        except:
            # S. Hill 2015-10-14: Eventually aospy will have a built-in land
            # mask array that it can use in case the object doesn't have one
            # of its own.  For now the object /must/ have one.
            raise
        if self.land_mask in (True, 'land'):
            return lmask
        if self.land_mask == 'ocean':
            return 1. - lmask
        if self.land_mask in ('strict_land', 'strict_ocean'):
            raise NotImplementedError
        return 1

    def ts(self, data, model):
        """Create time-series of region-average data."""
        data_masked = self.mask_var(data)
        dims = ['lat', 'lon']
        coords = [data_masked.coords[c] for c in dims]
        sfc_area = self._get_sfc_area(model, dims, coords)
        land_mask = self._get_land_mask(model, dims, coords)
        weights = (self.mask_var(sfc_area)*land_mask).sum('lat').sum('lon')
        return(data_masked*sfc_area*land_mask).sum('lat').sum('lon') / weights

    def av(self, data, model):
        """Time average of region-average time-series."""
        ts_ = self.ts(data, model)
        if 'year' not in ts_.coords:
            return ts_
        return ts_.mean('year')

    def std(self, data, model):
        """Standard deviation of region-average time-series."""
        ts_ = self.ts(data, model)
        if 'year' not in ts_.coords:
            return ts_
        return ts_.std('year')
