"""region.py: Region class and region_inst()."""
import numpy as np

from .utils import get_parent_attr


class Region(object):
    """Geographical region."""
    def __init__(self):
        pass

    def __str__(self):
        return 'Geographical region "' + self.name + '"'

    def _add_to_mask(self, mask, latb, lonb, lat, lon):
        """Return specified lat-lon rectangle as 2D grid."""
        lons, lats = np.meshgrid((lon > lonb[0]) & (lon < lonb[1]),
                                 (lat > latb[0]) & (lat < latb[1]))
        reg = lons*lats
        return np.where(reg, reg, mask)

    def make_mask(self, model):
        """Create region mask for the given model."""
        # Start with empty mask.
        lat = get_parent_attr(model, 'lat')
        lon = get_parent_attr(model, 'lon')
        mask = np.zeros((lat.size, lon.size))
        # Use region bounds stored in self.mask_bounds if available.
        try:
            # Add to mask for each set of bounds specified.
            for bounds in self.mask_bounds:
                mask = self._add_to_mask(mask, bounds[0], bounds[1],
                                         model.lat, model.lon)
        # Otherwise use self.lat_bnds and self.lon_bnds attributes.
        except AttributeError:
            mask = self._add_to_mask(mask, self.lat_bnds, self.lon_bnds,
                                     model.lat, model.lon)

        # Apply land or ocean mask as needed.
        finally:
            if model.land_mask is None:
                return mask
            elif self.land_mask in [True, 'land']:
                return mask*model.land_mask
            elif self.land_mask == 'strict_land':
                return mask*np.where(model.land_mask == 1., 1., 0.)
            elif self.land_mask == 'ocean':
                return mask*(1. - model.land_mask)
            elif self.land_mask == 'strict_ocean':
                return mask*np.where(model.land_mask == 0., 1., 0.)
            else:
                return mask

    def mask_var(self, data, model):
        """Mask the data of the given variable outside the region."""
        # Interpolate region to model grid.
        reg_mask = self.make_mask(model)
        # Mask input values where region mask is zero. Assumes dimensions are
        # (time, level, lat, lon), i.e. data.ndim=4.
        return np.ma.masked_where(
            np.tile(reg_mask == 0., (data.shape[0], data.shape[1], 1, 1)), data
        )

    def ts(self, data, model):
        """Create a time-series of region-average data."""
        # Mask the data outside the region, and flatten lat/lon into 1 dim.
        if data.ndim == 3:
            data = data[:,np.newaxis,:,:]
        data = self.mask_var(data, model)
        data = data.reshape(data.shape[0], data.shape[1], -1)
        # Get the region mask for the given model's grid.
        reg_mask = self.make_mask(model)
        # At gridpoints where the region is not totally masked, weight by that
        # point's surface area and the mask value.
        weights = np.ma.masked_where(reg_mask == 0,
                                     model.sfc_area*reg_mask).ravel()
        # Average over the region at each timestep and at each level.
        weights = np.tile(weights, (data.shape[0], data.shape[1], 1))
        avg = np.squeeze(np.ma.average(data, weights=weights, axis=-1))
        # If time axis had length 1, squeeze removed it, so put back on.
        if data.shape[0] == 1:
            avg = avg[np.newaxis,:]
        # If result is singleton, then turn into float, not numpy array.
        if isinstance(avg, np.ma.core.MaskedArray) and avg.ndim == 0:
            avg = float(avg)
        return avg

    def av(self, data, model):
        """Time average of region-average data."""
        ts_ = self.ts(data, model)
        if isinstance(ts_, float):
            return ts_
        else:
            mean = np.ma.mean(ts_, axis=0)
            return np.squeeze(mean)

    def std(self, data, model):
        """Standard deviation of time-series data."""
        ts_ = self.ts(data, model)
        if isinstance(ts_, float):
            return ts_
        else:
            stdev = np.ma.std(ts_, axis=0)
            return np.squeeze(stdev)
