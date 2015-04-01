"""
Package for atmospheric & oceanic science data analysis, management, and
visualization.
"""
def _set_attr_from_init_kwargs(obj, kwargs):
    """
    Given a dict of keyword arguments and their values as input for an
    __init__() call, set each attribute of the object with the name and value
    of the inputted arguments.
    """
    for kwarg_key, kwarg_value in kwargs.iteritems():
        setattr(obj, kwarg_key, kwarg_value)

def _get_parent_attr(obj, attr):
    """
    Check if the object has the given attribute and it is non-empty.  If not,
    check each parent object for the attribute and use the first one found.
    """
    # Get the attribute from the object itself, if it exists.
    try:
        return getattr(obj, attr)
    # Otherwise, loop through the parent objects until finding the attr.
    except AttributeError:
        for parent_obj in ('run', 'model', 'proj'):
            try:
                return getattr(getattr(obj, parent_obj), attr)
            except AttributeError:
                pass

def _set_named_attr_dict(obj, attr, attr_list):
    """Set attribute lists that are dicts of {name: value} type."""
    setattr(obj, attr, {item.name: item for item in attr_list})

def _set_parent_object_of_child(parent_obj, parent_attr, child):
    setattr(child, parent_attr, parent_obj)

def _set_nc_grid_objs(obj):
    """Set the nc_grid_objs of an aospy object."""
    from netCDF4 import Dataset, MFDataset
    nc_grid_paths = _get_parent_attr(obj, 'nc_grid_paths')
    nc_objs = []
    for path in nc_grid_paths:
        try:
            nc_obj = MFDataset(path)
        except IOError:
            nc_obj = Dataset(path)
        finally:
            nc_objs.append(nc_obj)
    setattr(obj, 'nc_grid_objs', tuple(nc_objs))

def _set_attr_from_nc_grid(obj, attr, attr_name):
    """Set attribute that comes from an nc_grid file."""
    for nc in obj.nc_grid_objs:
        try:
            val = nc.variables[attr_name][:]
        except KeyError:
            pass
        else:
            setattr(obj, attr, val)
            break
    # If not in any nc_grid file, set to None
    if not hasattr(obj, attr):
        setattr(obj, attr, None)
def _set_mult_nc_grid_attr(obj):
    """
    Set multiple attrs from grid file given their names in the grid file.
    """
    grid_attrs = {
        'lat':         ('lat', 'latitude', 'LATITUDE', 'y'),
        'lat_bounds':  ('latb', 'lat_bnds', 'lat_bounds'),
        'lon':         ('lon', 'longitude', 'LONGITUDE', 'x'),
        'lon_bounds':  ('lonb', 'lon_bnds', 'lon_bounds'),
        'level':       ('level', 'lev', 'plev'),
        'time':        ('time', 'TIME'),
        'time_st':     ('average_T1',),
        'time_end':    ('average_T2',),
        'time_dur':    ('average_DT',),
        'time_bounds': ('time_bounds', 'time_bnds'),
        'sfc_area':    ('area', 'sfc_area'),
        'zsurf':       ('zsurf',),
        'land_mask':   ('land_mask',),
        'pk':          ('pk',),
        'bk':          ('bk',),
        'phalf':       ('phalf',)
    }
    for attr, attr_names in grid_attrs.iteritems():
        for name in attr_names:
            _set_attr_from_nc_grid(obj, attr, name)

def calc_grid_sfc_area(lonb, latb, gridcenter=False):
    """Calculate surface area of each grid cell in a lon-lat grid."""
    import numpy as np
    from aospy.constants import r_e

    def diff_latlon_bnds(array):
        import numpy as np
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

def calc_levs_thick(levs):
    """
    Calculates the thickness, in Pa, of each pressure level.

    Assumes that the pressure values given are at the center of that model
    level, except for the lowest value (typically 1000 hPa), which is the
    bottom boundary. The uppermost level extends to 0 hPa.

    """
    import numpy as np
    # Bottom level extends from levs[0] to halfway betwen levs[0]
    # and levs[1].
    lev_thick = [0.5*(levs[0] - levs[1])]
    # Middle levels extend from halfway between [k-1], [k] and [k], [k+1].
    for k in range(1, levs.size-1):
        lev_thick.append(0.5*(levs[k-1] - levs[k+1]))
    # Top level extends from halfway between top two levels to 0 hPa.
    lev_thick.append(0.5*(levs[-2] + levs[-1]))
    # Convert to numpy array and from hectopascals (hPa) to Pascals (Pa).
    return np.array(lev_thick)*100.


class Proj(object):
    """
    Project parameters: models, regions, directories, etc.

    """
    def __init__(self, **kwargs):
        # Set attributes from inputted keyword arguments.
        _set_attr_from_init_kwargs(self, kwargs)
        # Set the vars dict if a var list was given.
        if 'vars' in kwargs.keys():
            _set_named_attr_dict(self, 'vars', kwargs['vars'])
            for var in self.vars.values():
                _set_parent_object_of_child(self, 'proj', var)

    def __str__(self):
        return 'Project instance "' + self.name + '"'

class Model(object):
    """Model and parameters of local data stored and desired."""
    def __init__(self, **kwargs):
        _set_attr_from_init_kwargs(self, kwargs)
        # Use the inputted names and netCDF filepath to create grid data.
        _set_nc_grid_objs(self)
        _set_mult_nc_grid_attr(self)
        self._set_sfc_area()
        # Populate the runs dictionary with the specified list of runs.
        if 'runs' in kwargs.keys():
            _set_named_attr_dict(self, 'runs', kwargs['runs'])
            # Create dict of default runs to use.
            _set_named_attr_dict(
                self, 'default_runs',
                kwargs.get('default_runs', kwargs['runs'])
                )

    def __str__(self):
        return 'Model instance "' + self.name + '"'

    def _set_levs_thick(self):
        """
        Set the 1D array holding the pressure thickness of the model levels.
        """
        if self.level:
            self.levs_thick = calc_levs_thick(self.level)
        else:
            self.levs_thick = None

    def _set_sfc_area(self):
        """Set the 2D array holding the surface area of gridboxes."""
        if getattr(self, 'sfc_area', None) is not None:
            return
        else:
            try:
                sfc_area = calc_grid_sfc_area(
                    self.lon_bounds, self.lat_bounds, gridcenter=False
                )
            except AttributeError:
                sfc_area = calc_grid_sfc_area(
                    self.lon, self.lat, gridcenter=True
                )
            self.sfc_area = sfc_area

    def _get_dt(self):
        """Get the model's time spacing."""
        import numpy as np
        if self.time_dur is not None:
            return self.time_dur
        elif self.time_bounds is not None:
            return self.time_bounds[:,1] - self.time_bounds[:,0]
        else:
            return np.ones(np.size(self.time))

    def _set_runs(self, runs_list):
        """Create the Model's dictionary of Runs."""
        self.runs = {run.name: run for run in runs_list}
        for run in runs_list:
            run.model = self

class Run(object):
    """Model run parameters."""
    def __init__(self, **kwargs):
        _set_attr_from_init_kwargs(self, kwargs)
        self._set_direc()

    def __str__(self):
        return 'Run instance "%s"' % self.name

    def _set_direc(self):
        """Set the list of paths containing the Run's netCDF data."""
        if all([hasattr(self, em) for em in 
                ['ens_mem_prefix', 'ens_mem_ext', 'ens_mem_suffix']]):
            self.direc_nc = [self.ens_mem_prefix + ext + self.ens_mem_suffix
                             for ext in self.ens_mem_ext]

class Var(object):
    """Model-native variables."""
    def __init__(self, **kwargs):
        _set_attr_from_init_kwargs(self, kwargs)
        if 'vars_list' in kwargs.keys():
            _set_named_attr_dict(self, 'vars', kwargs['vars_list'])

    def __str__(self):
        return 'Var instance "' + self.name + '"'

    def convert_to_plot_units(self, data):
        """
        Multiply the given data by the plotting units conversion if it exists.
        """
        try:
            return data*self.plot_units_conv
        except:
            pass

    def mask_unphysical(self, data):
        """Mask data array where values are outside physically valid range."""
        from numpy import ma
        try:
            return ma.masked_outside(data, self.valid_range[0],
                                     self.valid_range[1])
        except AttributeError:
            return data

    def _get_nc_one_dir(self, name, intvl_type):
        """Get the names of netCDF files when all in same directory."""
        direc_nc = _get_parent_attr(self, 'direc_nc')
        nc_files = _get_parent_attr(self, 'nc_files')
        try:
            return [direc_nc + '/' + nc_files[name]]
        except TypeError:
            return [direc_nc + '/' + nc_file for nc_file in nc_files[name]]
        except KeyError:
            if type(nc_files[intvl_type]) in (list, tuple):
                files = [direc_nc + '/' + nc_file for nc_file in
                         nc_files[intvl_type]]
            elif type(nc_files[intvl_type]) is str:
                files = [direc_nc + '/' + nc_files[intvl_type]]
            return files

    def _get_nc_gfdl_dir_struct(self, name, ens_mem, intvl, start_yr, 
                                end_yr, intvl_type, dtype):
        """Get the names of netCDF files stored in GFDL standard directory
        structure and names.
        """
        from aospy.io import nc_name_gfdl
        direc_nc = _get_parent_attr(self, 'direc_nc')
        domain = _get_parent_attr(self, 'domain')
        nc_dur = _get_parent_attr(self, 'nc_dur')
        nc_start_yr = _get_parent_attr(self, 'nc_start_yr')
        separator = '_' if dtype == 'av' else '/'
        if type(ens_mem) is int:
            direc_nc = direc_nc[ens_mem]
        direc = (direc_nc + '/' + domain + '/' + dtype + '/'  + 
                 intvl_type + separator + str(nc_dur) + 'yr/')
        files = [direc + nc_name_gfdl(name, domain, dtype, intvl_type,
                                      yr, intvl, nc_start_yr, nc_dur)
                 for yr in range(start_yr, end_yr + 1)]
        return files

    def _get_nc(self, ens_mem, intvl, start_yr, end_yr, intvl_type, dtype):
        """
        Create an MFDataset for the variable spanning all desired timesteps.

        Files chosen depend on the specified variables and time intvl and the
        attributes of the netCDF files.
        """
        from netCDF4 import MFDataset
        # from aospy.io import dmget_nc
        from aospy.io import hsmget_nc
        # If all files in one directory, need to use user-specified file
        # names for each data type or each variable.
        nc_dir_structure = _get_parent_attr(self, 'nc_dir_structure')
        # Cycle through possible names.
        names = [self.name]
        if hasattr(self, 'alt_names'):
            names += list(self.alt_names)
        for name in names:
            try:
                if nc_dir_structure == 'one_dir':
                    files = self._get_nc_one_dir(name, intvl_type)
                elif nc_dir_structure.lower() == 'gfdl':
                    files = self._get_nc_gfdl_dir_struct(
                        name, ens_mem, intvl, start_yr, 
                        end_yr, intvl_type, dtype
                    )
                # Remove duplicate entries.
                files_list = list(set(files))
                files_list.sort()
                # Retrieve files from archive using desired system calls.
                # dmget_nc(files_list)
                hsmget_retcode = hsmget_nc(files_list)
                # Return netCDF4.MFDataset object containing the data.
                return MFDataset(files_list)
            except (RuntimeError, KeyError):
                pass
        else:
            raise IOError("Could not find files for variable '%s'." % self.name)

    def data_yr_range(self, yr_range_in, dtype):
        """Get year range and 'effective' # of years of some given data."""
        if yr_range_in in ('default', None):
            yr_range = _get_parent_attr(self, 'default_yr_range')
            if yr_range == 'all':
                start_yr = _get_parent_attr(self, 'nc_start_yr')
                end_yr = _get_parent_attr(self, 'nc_end_yr')
            else:
                start_yr, end_yr = yr_range
        else:
            start_yr, end_yr = yr_range_in
        # If input data is time-series, need to know how many
        # years. If input data is time-averaged, there's only one time
        # index.
        if dtype == 'ts':
            num_yr = end_yr - start_yr + 1
        else:
            num_yr = 1
        return start_yr, end_yr, num_yr

    def _get_data_vals(self, ens_mem, level, intvl, intvl_type, dtype,
                       start_yr, end_yr, indices):
        """Get the desired time indices of a run's netCDF data."""
        import numpy as np
        from numpy import ma
        # Get index of desired level, if applicable.
        if level and self.def_vert:
            lev_ind = np.where(self.model.level == level)
        # Average over ensemble members if desired.
            
        if ens_mem == 'avg' and _get_parent_attr(self, 'ens_mem_prefix'):
            ens_avg_dat = []
            for i in range(len(_get_parent_attr(self, 'ens_mem_ext'))):
                nc = self._get_nc(i, intvl, start_yr, end_yr,
                                 intvl_type, dtype)
                nc_vals = nc.variables[self.name][indices]
                # Remove unphysical values that can occur in 'av' data.
                # if dtype == 'av':
                    # nc_vals = self.mask_unphysical(nc_vals)
                # Get desired vertical level.
                if level and self.def_vert:
                    ens_avg_dat.append(nc_vals[:,lev_ind])
                else:
                    ens_avg_dat.append(nc_vals)
            vals = ma.mean(ens_avg_dat, axis=0)
        # Otherwise just get the values from the desired ensemble member.
        else:
            # If ensemble averaging was input but the run doesn't have
            # ensemble members, then let the _get_nc() function know this.
            if ens_mem == 'avg':
                ens_mem_nc = None
            else:
                ens_mem_nc = ens_mem
            nc = self._get_nc(ens_mem_nc, intvl, start_yr,
                             end_yr, intvl_type, dtype)
            # Variable names can differ.
            try:
                vals = nc.variables[self.name][:]
            except KeyError:
                for alt_name in self.alt_names:
                    try:
                        vals = nc.variables[alt_name][:]
                    except KeyError:
                        pass
                    else:
                        break
            finally:
                vals = vals[indices]        
            # Remove unphysical values that can occur in 'av' data.
            if dtype == 'av':
                vals = self.mask_unphysical(vals)
            # Get desired vertical level.
            if level and self.def_vert:
                vals = vals[:,lev_ind]
        # Get durations of the desired timesteps.
        for dt_name in ('average_DT',):
            try:
                dt = nc.variables[dt_name][:]
            except:
                pass
            else:
                break
        else:
            dt = self.model._get_dt()
        nc.close()
        return vals, dt[indices]

    def var_time_avg(self, ens_mem, level, intvl, intvl_type,
                     dtype, yr_range, **kwargs):
        """Calculate time averages of data and save externally."""
        import numpy as np
        from numpy import ma
        from aospy.io import _intvl_indices_and_label, regions_calcs, \
                             data_av_stat, get_timesteps, save_file
        # Get start and end years and 'effective' number of year timesteps.
        start_yr, end_yr, num_yr = self.data_yr_range(yr_range, dtype)
        # Get time interval.
        label, intvl = _intvl_indices_and_label(intvl, intvl_type)
        # Get time indices of the input data corresponding to the
        # desired timesteps.
        nc_start_yr = _get_parent_attr(self, 'nc_start_yr')
        nc_end_yr = _get_parent_attr(self, 'nc_end_yr')
        nc_dur = _get_parent_attr(self, 'nc_dur')
        nc_start_month = _get_parent_attr(self, 'nc_start_month')
        if nc_start_month is None:
            nc_start_month = 1
        nc_end_month = _get_parent_attr(self, 'nc_end_month')
        if nc_end_month is None:
            nc_end_month = 12
        indices = get_timesteps(intvl, start_yr, end_yr, nc_start_yr,
                                nc_end_yr, nc_dur, intvl_type, dtype,
                                nc_start_month=nc_start_month,
                                nc_end_month=nc_end_month)
        # Get the desired data values and time durations of the timesteps.
        vals, dt = self._get_data_vals(ens_mem, level, intvl, intvl_type,
                                       dtype, start_yr, end_yr, indices)
        # Add empty dimensions to data without time- or
        # level-dependence. Makes subsequent syntax simpler.
        if not self.def_vert:
            vals = vals[:,np.newaxis]
        # Split time dimension into two: one for each year, and one
        # for each timestep within a year.
        dt_ts = np.reshape(dt, (num_yr, intvl.size))
        dt_ts = dt_ts[:,:,np.newaxis,np.newaxis,np.newaxis]
        # Reshape the data array using these new time dimensions.
        intvl_ts = vals.reshape((num_yr, intvl.size, vals.shape[-3],
                                 vals.shape[-2], vals.shape[-1]))
        # Average over all timesteps within each year to get an annual
        # time-series.
        loc_ts = ma.multiply(intvl_ts, dt_ts).sum(axis=1) / dt_ts.sum(axis=1)
        # Compute and save all the desired computations on the data.
        calcs_out = kwargs.get('calcs_out', ['av'])
        files = data_av_stat(loc_ts, calcs_out)
        for out_dtype, data in files.iteritems():
            save_file(np.squeeze(data), self.proj, self.model, self.run,
                      ens_mem, self, level, intvl, out_dtype, yr_range)
        # Region averages and other computations.
        regions = kwargs.get('regions', ['globe'])
        regions_calcs(regions, calcs_out, loc_ts, self.proj, self.model,
                      self.run, ens_mem, self, level,
                      intvl, yr_range, label)

class Region(object):
    """Geographical region."""
    def __init__(self, **kwargs):
        _set_attr_from_init_kwargs(self, kwargs)

    def __str__(self):
        return 'Geographical region "' + self.name + '"'

    def _add_to_mask(self, mask, latb, lonb, lat, lon):
        """Return specified lat-lon rectangle as 2D grid."""
        import numpy as np
        lons, lats = np.meshgrid((lon > lonb[0]) & (lon < lonb[1]),
                                 (lat > latb[0]) & (lat < latb[1]))
        reg = lons*lats
        return np.where(reg, reg, mask)

    def make_mask(self, model):
        """Create region mask for the given model."""
        import numpy as np
        # Start with empty mask.
        lat = _get_parent_attr(model, 'lat')
        lon = _get_parent_attr(model, 'lon')
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
                return  mask*np.where(model.land_mask == 0., 1., 0.)
            else:
                return mask

    def mask_var(self, data, var, model):
        """Mask the data of the given variable outside the region."""
        from numpy import tile
        from numpy.ma import masked_where
        # Interpolate region to model grid.
        reg_mask = self.make_mask(model)
        # Mask input values where region mask is zero. Assumes dimensions are
        # (time, level, lat, lon), i.e. data.ndim=4.
        return masked_where(tile(reg_mask == 0.,
                                 (data.shape[0], data.shape[1], 1, 1)), data)

    def ts(self, data, var, model, level):
        """Create a time-series of region-average data."""
        from numpy import tile, squeeze
        from numpy.ma import masked_where, average
        # Mask the data outside the region, and flatten lat/lon into 1 dim.
        data = self.mask_var(data, var, model)
        data = data.reshape(data.shape[0], data.shape[1], -1)
        # Get the region mask for the given model's grid.
        reg_mask = self.make_mask(model)
        # At gridpoints where the region is not totally masked, weight by that
        # point's surface area and the mask value.
        weights = masked_where(reg_mask == 0, model.sfc_area*reg_mask).ravel()
        # Tile the weights to be the same shape as the data. Required by the
        # numpy.ma.average function.
        weights = tile(weights, (data.shape[0], data.shape[1], 1))
        # Average over the region at each timestep and at each level.
        out = squeeze(average(data, weights=weights, axis=-1))
        if out.shape == (1,):
            out = out[0]
        return out

    def av(self, data, var, model, level):
        """Time average of region-average data."""
        from numpy import squeeze
        out =  squeeze(self.ts(data, var, model, level).mean(axis=0))
        if out.shape == (1,):
            out = out[0]
        return out


    def std(self, data, var, model, level):
        """Standard deviation of time-series data."""
        from numpy import squeeze
        out = squeeze(self.ts(data, var, model, level).std(axis=0))
        if out.shape == (1,):
            out = out[0]
        return out

import av_stat, calcs, constants, io, plotting, regions, variables, projects
from print_table import print_table
