"""region.py: Region class and region_inst()."""
import xray


class Region(object):
    """Geographical region."""
    def __init__(self):
        """Instantiate a Region object."""
        self.name = ''
        self.lon_bounds = []
        self.lat_bounds = []
        self.mask_bounds = []

    def __str__(self):
        return 'Geographical region "' + self.name + '"'

    __repr__ = __str__

    def _add_to_mask(self, data, latb, lonb):
        mask_lat = (data['lat'] > latb[0]) & (data['lat'] < latb[1])
        mask_latlon = mask_lat & ((data['lon'] > lonb[0]) &
                                  (data['lon'] < lonb[1]))
        return mask_latlon

    def make_mask(self, data):
        """Construct the mask that defines this region."""
        # For each set of bounds add to the conditional.
        mask = False
        try:
            for bounds in self.mask_bounds:
                mask |= self._add_to_mask(data, bounds[0], bounds[1])
        except:
            mask |= self._add_to_mask(data, self.lat_bnds, self.lon_bnds)

        # No landmask for now.
        return mask

    def mask_var(self, data):
        """Mask the data of the given variable outside the region."""
        return data.where(self.make_mask(data))

    def ts(self, data, model):
        """Create time-series of region average-data."""
        data_masked = self.mask_var(data)
        dims = ['lat', 'lon']
        coords = [data_masked.coords[c] for c in dims]
        sfc_area = xray.DataArray(model.sfc_area, dims=dims, coords=coords)
        weights = self.mask_var(sfc_area)
        # Take the area average
        return ((data_masked*model.sfc_area).sum('lat').sum('lon') /
                weights.sum('lat').sum('lon'))

    def av(self, data, model):
        """ Time average of region-average data."""
        ts_ = self.ts(data, model)
        if 'year' not in ts_.coords:
            return ts_
        else:
            return ts_.mean('year')

    def std(self, data, model):
        """Standard deviation of time-series data"""
        ts_ = self.ts(data, model)
        if 'year' not in ts_.coords:
            return ts_
        else:
            return ts_.std('year')
