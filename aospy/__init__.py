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
        for parent_obj in ('var', 'run', 'model', 'proj'):
            try:
                return getattr(getattr(obj, parent_obj), attr)
            except AttributeError:
                pass

def _set_named_attr_dict(obj, attr, attr_list):
    """Set attribute lists that are dicts of {name: value} type."""
    setattr(obj, attr, {item.name: item for item in attr_list})

def _set_parent_object_of_child(parent_obj, parent_attr, child):
    setattr(child, parent_attr, parent_obj)

def _get_nc_grid_objs(obj):
    """Get the nc_grid_objs of an aospy object."""
    from netCDF4 import Dataset, MFDataset
    nc_grid_paths = _get_parent_attr(obj, 'nc_grid_paths')
    nc_objs = []
    for path in nc_grid_paths:
        try:
            nc_obj = Dataset(path)
        except TypeError:
            nc_obj = MFDataset(path)
        finally:
            nc_objs.append(nc_obj)
    return tuple(nc_objs)

def _set_attr_from_nc_grid(nc_grid_obj, obj, attr, attr_name):
    """Set attribute that comes from an nc_grid file."""
    for nc in nc_grid_obj:
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
    nc_grid_obj = _get_nc_grid_objs(obj)
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
        'phalf':       ('phalf',),
        'pfull':       ('pfull',)
    }
    for attr, attr_names in grid_attrs.iteritems():
        for name in attr_names:
            _set_attr_from_nc_grid(nc_grid_obj, obj, attr, name)
    # Close file objects.
    for nc in nc_grid_obj:
        nc.close()

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
    from aospy.calcs import to_pascal
    # Bottom level extends from levs[0] to halfway betwen levs[0]
    # and levs[1].
    levs = to_pascal(levs)
    lev_thick = [0.5*(levs[0] - levs[1])]
    # Middle levels extend from halfway between [k-1], [k] and [k], [k+1].
    for k in range(1, levs.size-1):
        lev_thick.append(0.5*(levs[k-1] - levs[k+1]))
    # Top level extends from halfway between top two levels to 0 hPa.
    lev_thick.append(0.5*(levs[-2] + levs[-1]))
    # Convert to numpy array and from hectopascals (hPa) to Pascals (Pa).
    return np.array(lev_thick)


class Proj(object):
    """
    Project parameters: models, regions, directories, etc.

    """
    def __init__(self, verbose=True, **kwargs):
        import time
        self.verbose = verbose
        if self.verbose:
            print ("Initializing Project instance: %s (%s)"
                   % (kwargs['name'], time.ctime()))
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
    def __init__(self, units=False, **kwargs):
        """Model-native variables."""
        _set_attr_from_init_kwargs(self, kwargs)
        if 'vars_list' in kwargs.keys():
            _set_named_attr_dict(self, 'vars', kwargs['vars_list'])
        # Copy units from 
        if type(units) is Units:
            self._Units = units
            for var_attr, units_attr in zip(
                    ('units', 'plot_units', 'plot_units_conv', 'vert_int_units',
                     'vert_int_plot_units', 'vert_int_plot_units_conv'),
                    ('units', 'plot', 'plot_conv', 'vert_int',
                     'vert_int_plot', 'vert_int_plot_conv')
            ):
                setattr(self, var_attr, getattr(units, units_attr))
        else:
            self.units = units

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

class Units(object):
    vint_str = r'kg m$^{-2}$'
    def __init__(self, units='', plot=False, plot_conv=1., vert_int=False,
                 vert_int_plot=False, vert_int_plot_conv=False):
        """String representation of physical units and conversion methods.""" 
        self.units = units
        self.plot = plot if plot else units
        self.plot_conv = plot_conv
        if vert_int:
            self.vert_int = vert_int
        else:
            self.vert_int = ' '.join(
                [Units.vint_str, units]).replace('  ', ' ')
        if vert_int_plot:
            self.vert_int_plot = vert_int_plot
        else:
            self.vert_int_plot = ' '.join(
                [Units.vint_str, self.plot]).replace('  ', ' ')
        if vert_int_plot_conv:
            self.vert_int_plot_conv = vert_int_plot_conv
        else:
            self.vert_int_plot_conv = plot_conv


class Calc(object):
    def __init__(self, proj=None, model=None, run=None, ens_mem=None, var=None,
                 yr_range=None, region=None, intvl_in=None, intvl_out=None,
                 dtype_in_time=None, dtype_in_vert=None, dtype_out_time=None,
                 dtype_out_vert=None, level=None, skip_time_inds=False, 
                 yr_chunk_len=False, verbose=True):
        """Class for executing, saving, and loading a single computation."""
        from aospy.io import (_proj_inst, _model_inst, _run_inst, _var_inst,
                              _month_indices)
        # Turn strings into tuples.
        if type(proj) in (str, Proj):
            proj = tuple([proj])
        if type(model) in (str, Model):
            model = tuple([model])
        if type(run) in (str, Run):
            run = tuple([run])
        # Make tuples the same length.
        if len(proj) == 1 and (len(model) > 1 or len(run) > 1):
            proj = tuple(list(proj)*len(run))
        if len(model) == 1 and len(run) > 1:
            model = tuple(list(model)*len(run))
        assert len(model) == len(run)
        assert len(proj) == len(model)

        # Convert string names to aospy objects.
        self.proj  = tuple([_proj_inst(pr) for pr in proj])
        self.model = tuple([_model_inst(mod, pr) for
                           (mod, pr) in zip(model, self.proj)])
        self.run   = tuple([_run_inst(rn, mod, pr) for (rn, mod, pr)
                            in zip(run, self.model, self.proj)])

        self.proj_str = '_'.join(set([p.name for p in self.proj]))
        self.model_str = '_'.join(set([m.name for m in self.model]))
        run_names = [r.name for r in self.run]
        self.run_str = '_'.join(set(run_names))
        self.run_str_full = '_'.join(run_names)

        self.var = _var_inst(var)
        self.name = self.var.name
        self.domain = self.var.domain

        self.verbose = verbose
        self._print_verbose("\nInitializing Calc instance: %s", self.__str__())

        self._set_nc_attrs()

        # if not self.nc_start_month:
            # self.nc_start_month = 1
        # if not self.nc_end_month:
            # self.nc_end_month = 12
        if type(ens_mem) is int:
            self.direc_nc = self.direc_nc[ens_mem]

        # Some vars are computed as functions of other vars.
        try:
            self.function = self.var.func
        except AttributeError:
            self.function = lambda x: x
        try:
            self.vars = self.var.vars
        except AttributeError:
            self.vars = (self.var,)

        self.ens_mem = ens_mem
        self.level = level
        self.intvl_in = intvl_in
        self.intvl_out = intvl_out
        self.dtype_in_time = dtype_in_time
        self.dtype_in_vert = dtype_in_vert        
        if type(dtype_out_time) in (list, tuple):
            self.dtype_out_time = tuple(dtype_out_time)
        else:
            self.dtype_out_time = tuple([dtype_out_time])
        self.dtype_out_vert = dtype_out_vert
        self.region = region

        self.yr_range = yr_range
        self.start_yr, self.end_yr = self._get_yr_range()
        self.num_yr = self._get_num_yr()
        self.months = _month_indices(intvl_out)
        self.yr_chunk_len = yr_chunk_len
        self.chunk_ranges = self._make_time_chunks()

        self.dir_scratch = self._dir_scratch()
        self.dir_archive = self._dir_archive()
        self.file_name = {d: self._file_name(d) for d in self.dtype_out_time}
        self.path_scratch = {d: '/'.join([self.dir_scratch,
                                          self.file_name[d]]).replace('//', '/')
                             for d in self.dtype_out_time}
        self.path_archive = {d: '/'.join([self.dir_archive,
                                          self.file_name[d]]).replace('//', '/')
                             for d in self.dtype_out_time}

        if not skip_time_inds:
            self._set_time_dt()

        self.data_out = {}

    def __str__(self):
        """String representation of the object."""
        return "Calc object: " + ', '.join(
            (self.name, self.proj_str, self.model_str, self.run_str_full)
        )

    def __call__(self, dtype_out):
        try:
            data = self.data_out[dtype_out]
        except AttributeError:
            raise AttributeError("'data_out' attribute has not been set for "
                                 "aospy.Calc object '%s'." % self)
        except KeyError:
            raise KeyError("No data exists for dtype '%s' for aospy.Calc "
                           "object '%s'." % (dtype_out, self))
        else:
            return data

    def __add__(self, other):
        """Add aospy.Calc.data_out, numpy.array, int, or float to the array."""
        import numpy as np
        try:
            data_self = getattr(self, 'data_out')
        except AttributeError:
            raise AttributeError("Object '%s' lacks a 'data_out' attr" % self)
        # aospy.Calc, numpy array, and int/flot objects each use a different
        # named attribute to store their values
        for attr in ('data_out', '__array__', 'real'):
            try:
                data_other = getattr(other, 'data_out')
            except AttributeError:
                pass
            else:
                break
        return np.add(data_self, data_other)
    
    def _print_verbose(self, *args):
        """Print diagnostic message."""
        import time
        if not self.verbose:
            pass
        else:
            try:
                print args[0] % args[1], '(%s)' % time.ctime()
            except IndexError:
                print args[0], '(%s)' % time.ctime()
        
    def _get_yr_range(self):
        """Set the object's span of years."""
        if self.yr_range == 'default':
            start_yr, end_yr = _get_parent_attr(self.run[0], 'default_yr_range')
        elif self.yr_range == 'all':
            start_yr = _get_parent_attr(self.run[0], 'nc_start_yr')
            end_yr = _get_parent_attr(self.run[0], 'nc_end_yr')
        else:
            start_yr, end_yr = self.yr_range
        return start_yr, end_yr

    def _get_num_yr(self):
        """Compute effective number of years in the input data."""
        if self.dtype_in_time in ('ts', 'inst'):
            num_yr = self.end_yr - self.start_yr + 1
        else:
            num_yr = 1
        return num_yr

    def _set_nc_attrs(self):
        for attr in ('nc_start_yr', 'nc_end_yr', 'nc_dur', 'direc_nc',
                     # 'nc_start_month', 'nc_end_month', 
                     # 'ens_mem_prefix', 'ens_mem_ext'
                     'nc_files', 'nc_dir_struc', 'default_yr_range'):
            attr_val = tuple([_get_parent_attr(rn, attr) for rn
                              in self.run])
            setattr(self, attr, attr_val)

    def _set_time_dt(self):
        """Get time and dt arrays at needed time indices."""
        from aospy.io import _get_time as io_get_time
        # Use the first var in the list that is an aospy.Var object.
        nc_var = self.vars[0]
        for var in self.vars:
            if type(var) is Var:
                nc_var = var
                break
        with self._get_nc(nc_var, self.start_yr, self.end_yr) as nc:
            time_obj = nc.variables['time']
            inds, time = io_get_time(
                time_obj[:], time_obj.units, time_obj.calendar,
                self.start_yr, self.end_yr, self.months, indices=True
            )
            self.time = time
            self.time_units = time_obj.units
            self.calendar = time_obj.calendar
            self.time_inds = inds
            self.dt = self._get_dt(nc, inds)

    def _get_dt(self, nc, indices):
        """Get durations of the desired timesteps."""
        import numpy as np
        if self.dtype_in_time == 'inst':
            return np.ones(np.shape(indices))
        else:
            for dt_name in ('average_DT',):
                try:
                    dt = nc.variables[dt_name]
                except KeyError:
                    pass
                else:
                    break
            else:
                dt = self.model[0]._get_dt()
            return dt[indices]

    def _reshape_time_indices(self, array, start_yr, end_yr):
        """Reshape time array to have year- and within-year axes.

        2015-04-23: This might not work for sub-monthly data spanning leap
        years or other calendar idiosyncracies.  For example, suppose using
        3 hourly data for DJF spanning a leap year and non leap-year.  The
        leap year February will have 4 more timesteps than the non-leap year.
        """
        import numpy as np
        reshaped = np.reshape(array, (end_yr - start_yr + 1, -1))
        return reshaped[:,:,np.newaxis,np.newaxis,np.newaxis]

    def _make_time_chunks(self):
        """Create tuple of (start, end) pairs based on specified year chunks."""
        import numpy as np
        if self.yr_chunk_len:
            dur = self.yr_chunk_len - 1
            st_yrs = range(self.start_yr, self.end_yr + 1, self.yr_chunk_len)
            end_yrs = range(self.start_yr + dur, self.end_yr + dur + 1,
                            self.yr_chunk_len)
            if len(end_yrs) == len(st_yrs) - 1:
                end_yrs.append(self.end_yr)
            dt_yrs = np.subtract(end_yrs, st_yrs)
        else:
            st_yrs, end_yrs = [self.start_yr], [self.end_yr]
        return zip(st_yrs, end_yrs)

    def _dir_scratch(self):
        """Create the string of the data directory on the scratch filesystem."""
        import os
        from aospy.io import _ens_label
        ens_label = _ens_label(self.ens_mem)
        dir_scratch = '/'.join(
            ['/work', os.getenv('USER'), self.proj_str, self.model_str,
             self.run_str, ens_label, self.var.name]
        ).replace('//', '/')
        return dir_scratch

    def _dir_archive(self):
        """Create the string of the data directory on the archive filesystem."""
        import os
        from aospy.io import _ens_label
        ens_label = _ens_label(self.ens_mem)
        dir_archive = '/'.join(
            ['/archive', os.getenv('USER'), self.proj_str, 'data',
             self.model_str, self.run_str, ens_label]
        ).replace('//', '/')
        return dir_archive

    def _file_name(self, dtype_out_time):
        """Create the name of the aospy file."""
        import aospy.io as io
        out_lbl = io._data_out_label(self.intvl_out, dtype_out_time, 
                                     dtype_vert=self.dtype_out_vert)
        in_lbl = io._data_in_label(self.intvl_in, self.dtype_in_time,
                                   self.dtype_in_vert)
        ens_lbl = io._ens_label(self.ens_mem)
        yr_lbl = io._yr_label((self.start_yr, self.end_yr))
        file_name = '.'.join(
            [self.name, out_lbl, in_lbl, self.model_str, self.run_str_full,
             ens_lbl, yr_lbl, 'p']
        ).replace('..','.')
        return file_name

    def _get_nc_one_dir(self, name, n=0):
        """Get the names of netCDF files when all in same directory."""
        try:
            files = [self.direc_nc[n] + '/' + self.nc_files[n][name]]
        except TypeError:
            files = [self.direc_nc[n] + '/' + nc_file for
                     nc_file in self.nc_files[n][name]]
        # except KeyError:
            # if type(self.nc_files[n][self.intvl_in]) in (list, tuple):
                # files = [self.direc_nc[n] + '/' + nc_file for nc_file in
                         # self.nc_files[n][self.intvl_in]]
            # elif type(self.nc_files[self.intvl_in]) is str:
                # files = [self.direc_nc[n] + '/' +
                         # self.nc_files[n][self.intvl_in]]
        # Remove duplicate entries.
        files = list(set(files))
        files.sort()
        return files

    def _get_nc_gfdl_dir_struct(self, name, start_yr, end_yr, n=0):
        """
        Get the names of netCDF files stored in GFDL standard directory
        structure and names.
        """
        from aospy.io import nc_name_gfdl
        domain = self.domain
        dtype_lbl = self.dtype_in_time
        if self.dtype_in_vert == 'sigma' and name != 'ps':
            domain += '_level'
        if self.dtype_in_time == 'inst':
            domain += '_inst'
            dtype_lbl = 'ts'
        # 2015-05-12 using 'ts' or 'inst' data only for now
        # separator = '_' if self.dtype_in_time == 'av' else '/'
        separator = '/' 
        if 'av_from_' in self.dtype_in_time:
            dtype = self.dtype_in_time.replace('av_from_', '')
            dtype_lbl = dtype
        else:
            dtype = self.dtype_in_time
        direc = (self.direc_nc[n] + '/' + domain + '/' + dtype_lbl + '/' +
                 self.intvl_in + separator + str(self.nc_dur[n]) + 'yr/')
        files = [direc + nc_name_gfdl(name, domain, dtype,
                                      self.intvl_in, yr, self.intvl_out,
                                      self.nc_start_yr[n], self.nc_dur[n])
                 for yr in range(start_yr, end_yr + 1)]
        # Remove duplicate entries.
        files = list(set(files))
        files.sort()
        return files

    def _get_nc(self, var, start_yr=False, end_yr=False, n=0):
        """
        Create an MFDataset for the variable spanning all desired timesteps.

        Files chosen depend on the specified variables and time intvl and the
        attributes of the netCDF files.
        """
        from netCDF4 import MFDataset
        from aospy.io import hsmget_nc, dmget_nc
        # Cycle through possible names.
        names = [var.name]
        if hasattr(var, 'alt_names'):
            names += list(var.alt_names)
        # Get directory and file info from the corresponding Run.
        # if n:
        # Get the data.
        for name in names:
            if self.nc_dir_struc[n] == 'one_dir':
                files = self._get_nc_one_dir(name, n=n)
            elif self.nc_dir_struc[0].lower() == 'gfdl':
                files = self._get_nc_gfdl_dir_struct(name, start_yr,
                                                     end_yr, n=n)
            try:
                # Retrieve files from archive using desired system calls.
                dmget_nc(files)
                # hsmget_retcode = hsmget_nc(files)
                # Return netCDF4.MFDataset object containing the data.
                return MFDataset(files)
            except (RuntimeError, KeyError):
                pass
        else:
            raise IOError("Could not find files for variable '%s'." % var)

    def _get_pressure_vals(self, var, start_yr, end_yr, n=0):
        """Get pressure array, whether sigma or standard levels."""
        from aospy.io import _var_inst
        from aospy.calcs import pfull_from_sigma, dp_from_sigma
        self._print_verbose("Getting pressure data: %s", var)
        if self.dtype_in_vert == 'pressure':
            if var == 'p':
                data = self.model[n].level
            elif var == 'dp':
                data = calc_levs_thick(self.model[n].level)
        if self.dtype_in_vert == 'sigma':
            bk = self.model[n].bk
            pk = self.model[n].pk
            ps_obj = _var_inst('ps')
            ps = self._get_var_data(ps_obj, start_yr, end_yr, eddy=False)
            if var == 'p':
                data = pfull_from_sigma(bk, pk, ps)
            elif var == 'dp':
                data = dp_from_sigma(bk, pk, ps)
        return data

    def _get_data_subset(self, data, region=False, eddy=False, time=False,
                         vert=False, lat=False, lon=False, n=0):
        """Subset the data array to the specified time/level/lat/lon, etc."""
        import numpy as np
        if region:
            if type(region) is str:
                data = data[region]
            elif type(region) is Region:
                data = data[region.name]
        if np.any(time):
            data = data[time]
            if 'av_from_' in self.dtype_in_time:
                data = np.mean(data, axis=0)[np.newaxis,:]
        if np.any(eddy):
            data -= np.ma.average(data, axis=0, weights=eddy)
        if np.any(vert):
            if self.dtype_in_vert != 'sigma':
                if np.max(self.model[n].level) > 1e4:
                    # Convert from Pa to hPa.
                    lev_hpa = self.model[n].level*1e-2
                else:
                    lev_hpa = self.model[n].level
                level_index = np.where(lev_hpa == self.level)
                if 'ts' in self.dtype_out_time:
                    data = np.squeeze(data[:,level_index])
                else:
                    data = np.squeeze(data[level_index])
        if np.any(lat):
            pass
        if np.any(lon):
            pass
        return data

    def _get_var_data(self, var, start_yr, end_yr, n=0, region=False,
                      eddy=False, vert=False, lat=False, lon=False):
        """Get the needed data from one aospy.Var."""
        import numpy as np
        from aospy.io import _get_time
        self._print_verbose("\tGetting data from netCDF files: %s", var)
        with self._get_nc(var, start_yr, end_yr, n=n) as nc:
            # Variable names can differ.
            try:
                data = nc.variables[var.name][:]
            except KeyError:
                for alt_name in var.alt_names:
                    try:
                        data = nc.variables[alt_name][:]
                    except KeyError:
                        pass
                    else:
                        break
            time = nc.variables['time']
            t_array, t_units, t_cal = time[:], time.units, time.calendar
            t_inds = _get_time(t_array, t_units, t_cal, start_yr, end_yr,
                             self.months, indices='only')
            if eddy:
                eddy = self._get_dt(nc, t_inds)
        data = self._get_data_subset(
            data, region=region, eddy=eddy, time=t_inds,
            vert=self.level, lat=lat, lon=lon
        )
        # Interpolate data at sigma half levels to full levels.
        if self.dtype_in_vert == 'sigma' and var.def_vert == 'phalf':
            data = 0.5*(data[:,:-1] + data[:,1:])
        # if self.dtype_in_time == 'av':
            # data = self.mask_unphysical()
        # To simplify broadcasting, add dim to non-vertical data.
        if not var.def_vert:
            data = data[:,np.newaxis]
        return data

    def _get_all_vars_data(self, start_yr, end_yr, eddy=False):
        """Get the needed data from all of the vars in the calculation."""
        import numpy as np
        all_vals = []
        for n, var in enumerate(self.vars):
            # If only 1 run, use it to load all data.
            # Otherwise assume that # runs == # vars to load.
            if len(self.run) == 1:
                n = 0
            # Pressure handled specially due to complications from sigma vs. p.
            if var in ('p', 'dp'):
                data = self._get_pressure_vals(var, start_yr, end_yr)
            # Pass numerical constants as is.
            elif type(var) in (float, int):
                data = var
            # aospy.Var objects remain.
            elif var.name in ('lat', 'lon', 'level', 'pk',
                              'bk' 'sfc_area', 'time'):
                # Get grid, time, etc. arrays directly from model object
                data = getattr(self.model[0], var.name)
            else:
                data = self._get_var_data(var, start_yr, end_yr, n=n, eddy=eddy)
            all_vals.append(data)
        return all_vals

    def _local_ts(self, dp, dt, num_yr, *data_in):
        """Compute the function on the given gridded input data."""
        import numpy as np
        from aospy.calcs import int_dp_g
        result = self.function(*data_in)
        # Apply spatial reductions methods.
        if self.dtype_out_vert == 'vert_int':
            result = int_dp_g(result, dp)[:,np.newaxis,:,:]
        # If already averaged, pass data on.  Otherwise do time averaging.
        if 'av' in self.dtype_in_time:
            return result
        try:
            result = result.reshape((num_yr, -1, result.shape[-3],
                                     result.shape[-2], result.shape[-1]))
        except ValueError:
            result = result.reshape((num_yr, -1, result.shape[-2], 
                                     result.shape[-1]))
        # Average within each year, yielding a yearly time-series.
        try:
            return np.ma.multiply(result, dt).sum(axis=1) / dt.sum(axis=1)
        except ValueError:
            dt = dt[:,:,:,:,0]
            return np.ma.multiply(result, dt).sum(axis=1) / dt.sum(axis=1)

    def _time_reduce(self, loc_ts):
        """Compute all desired calculations on a local time-series."""
        files = {}
        if 'ts' in self.dtype_out_time:
            files.update({'ts': loc_ts})
        if 'av' in self.dtype_out_time:
            files.update({'av': loc_ts.mean(axis=0)})
        if 'eddy.av' in self.dtype_out_time:
            files.update({'eddy.av': loc_ts.mean(axis=0)})
        if 'std' in self.dtype_out_time:
            files.update({'std': loc_ts.std(axis=0)})
        if 'eddy.std' in self.dtype_out_time:
            files.update({'eddy.std': loc_ts.std(axis=0)})
        # Zonal asymmetry.
        if any('zasym' in out_type for out_type in self.dtype_out_time):
            # '.T'=transpose; makes numpy broadcasting syntax work.
            znl_ts = loc_ts.mean(axis=-1)
            zasym_ts = (loc_ts.T - znl_ts.T).T
            if 'zasym.ts' in self.dtype_out_time:
                files.update({'zasym.ts': zasym_ts})
            if 'zasym.av' in self.dtype_out_time:
                files.update({'zasym.av': zasym_ts.mean(axis=0)})
            if 'zasym.std' in self.dtype_out_time:
                files.update({'zasym.std': zasym_ts.std(axis=0)})
        # Zonal mean.
        if any('znl' in out_type for out_type in self.dtype_out_time):
            if 'znl.ts' in self.dtype_out_time:
                files.update({'znl.ts': znl_ts})
            if 'znl.av' in self.dtype_out_time:
                files.update({'znl.av': znl_ts.mean(axis=0)})
            if 'znl.std' in self.dtype_out_time:
                files.update({'znl.std': znl_ts.std(axis=0)})
        return files
    
    def region_calcs(self, loc_ts, eddy=False, n=0):
        """Region-averaged computations.  Execute and save to external file."""
        from aospy.io import _region_inst
        calcs_reg = ('ts', 'av', 'std')
            #, 'spvar.ts', 'spvar.av', 'znl.spvar.ts',
            # 'znl.spvar.av', 'zasym.ts', 'zasym.av', 'zasym.std',
            # 'zasym.spvar.ts', 'zasym.spvar.av'
        regions = [_region_inst(reg) for reg in self.region]
        # Perform each calculation for each region.
        for calc in calcs_reg:
            calc_name = ('reg.' + calc)
            if eddy:
                calc_name = '.'.join(['eddy', calc_name])
            if calc_name in self.dtype_out_time:
                reg_dat = {}
                for region in regions:
                    # Just pass along the data if averaged already.
                    if 'av' in self.dtype_in_time:
                        data_out = region.ts(loc_ts, self.model[n])
                    # Otherwise perform the calculation.
                    else:
                        method = getattr(region, calc)
                        data_out = method(loc_ts, self.model[n])
                    reg_dat.update({region.name: data_out})
                self.save(reg_dat, calc_name)

    def _compute_chunk(self, start_yr, end_yr, eddy=False):
        """Perform the calculation on the given chunk of times."""
        import time
        from aospy.io import _get_time
        self._print_verbose("Computing desired timeseries from netCDF data for "
                            "years %d-%d.", (start_yr, end_yr))
        inds = _get_time(
            self.time, self.time_units, self.calendar,
            start_yr, end_yr, self.months, indices='only'
        )
        dt = self.dt[inds]
        dt = self._reshape_time_indices(dt, start_yr, end_yr)
        data_in = self._get_all_vars_data(start_yr, end_yr, eddy=eddy)
        if self.dtype_out_vert == 'vert_int':
            dp = self._get_pressure_vals('dp', start_yr, end_yr)
        else:
            dp = False
        return self._local_ts(dp, dt, end_yr - start_yr + 1, *data_in)

    def compute(self, eddy=False):
        """Perform all desired calculations on the data and save externally."""
        import time
        import numpy as np
        # Compute the local time series for each chunk and then combine chunks.
        if all(['eddy' in do for do in self.dtype_out_time]) and eddy is False:
            self._print_verbose("Computing and saving eddy outputs.")
            eddy = True
        full_ts = [self._compute_chunk(start_yr, end_yr, eddy=eddy)
                   for start_yr, end_yr in self.chunk_ranges]
        if len(full_ts) == 1:
            full_ts = full_ts[0]
        else:
            self._print_verbose("Combining output data from all time chunks.")
            full_ts = np.ma.concatenate(full_ts, axis=0)
        # Apply time reduction methods on gridded data and save.
        self._print_verbose("Applying desired time-reduction methods.")
        reduced = self._time_reduce(full_ts)
        self._print_verbose("Writing desired gridded outputs to disk.")
        for dtype_out_time, data in reduced.iteritems():
            self.save(np.squeeze(data), dtype_out_time, 
                      dtype_out_vert=self.dtype_out_vert)
        # Apply time reduction methods to regional averages and save.
        if any(['reg' in do for do in self.dtype_out_time]) and self.region:
            self._print_verbose("Computing and saving regional outputs.")
            self.region_calcs(full_ts, eddy=eddy)
        # Perfom eddy computations.
        if any(['eddy' in do for do in self.dtype_out_time]) and eddy is False:
            self._print_verbose("Computing and saving eddy outputs.")
            self.compute(eddy=True)
    
    def _save_to_scratch(self, data, dtype_out_time, dtype_out_vert=False):
        """Save the data to the scratch filesystem."""
        import os
        import cPickle
        path = self.path_scratch[dtype_out_time]
        if not os.path.isdir(self.dir_scratch):
            os.makedirs(self.dir_scratch)
        with open(path, 'a+') as file_scratch:
            # Update, rather than overwrite, existing regional data.
            if 'reg' in dtype_out_time:
                # Open the existing dictionary if it exists.
                try:
                    reg_data = cPickle.load(file_scratch)
                except EOFError:
                    reg_data = {}
                # Add the new data to the dictionary.
                reg_data.update(data)
                data_out = reg_data
            else:
                data_out = data
        with open(path, 'w') as file_scratch:
            cPickle.dump(data_out, file_scratch)

    def _save_to_archive(self, dtype_out_time, dtype_out_vert=False):
        """Add the data to the tar file in /archive."""
        import os
        import tarfile
        if not os.path.isdir(self.dir_archive):
            os.makedirs(self.dir_archive)
        path = self.path_archive[dtype_out_time]
        with tarfile.open(path + 'data.tar', 'a') as tar:
            tar.add(self.path_scratch[dtype_out_time], arcname=path)

    def _update_data_out(self, data, dtype):
        """Append the data of the given dtype_out to the data_out attr."""
        try:
            self.data_out.update({dtype: data})
        except AttributeError:
            self.data_out = {dtype: data}

    def save(self, data, dtype_out_time, dtype_out_vert=False, 
             scratch=True, archive=True):
        """Save aospy data to data_out attr and to an external file."""
        self._update_data_out(data, dtype_out_time)
        if scratch:
            self._save_to_scratch(data, dtype_out_time, 
                                  dtype_out_vert=dtype_out_vert)
        if archive:
            self._save_to_archive(dtype_out_time, 
                                  dtype_out_vert=dtype_out_vert)
        print '\t%s' % self.path_scratch[dtype_out_time]

    def _load_from_scratch(self, dtype_out_time, dtype_out_vert=False):
        """Load aospy data saved on scratch file system."""
        import cPickle
        with open(self.path_scratch[dtype_out_time], 'r') as data:
            data_vals = cPickle.load(data)
        return data_vals

    def _load_from_archive(self, dtype_out_time, dtype_out_vert=False):
        """Load data save in tarball on archive file system."""
        import subprocess
        import tarfile
        path = '/'.join([self.dir_archive, 'data.tar']).replace('//','/')
        subprocess.call(['dmget'] + [path])
        with tarfile.open(path, 'r') as data_tar:
            data_vals = load(
                data_tar.extractfile(self.file_name[dtype_out_time])
            )
        return data_vals
    
    def load(self, dtype_out_time, dtype_out_vert=False, region=False, 
             time=False, vert=False, lat=False, lon=False, plot_units=False):
        """Load the data from the object if possible or from disk."""
        # Grab from the object if its there.
        try:
            data = self.data_out[dtype_out_time]
        except (AttributeError, KeyError):
            # Otherwise get from disk.  Try scratch first, then archive.
            try:
                data = self._load_from_scratch(dtype_out_time, dtype_out_vert)
            except OSError:
                data = self._load_from_archive(dtype_out_time, dtype_out_vert)
            # Subset the array and convert units as desired.
            if any((region, time, vert, lat, lon)):
                data = self._get_data_subset(data, region=region, time=time,
                                             vert=vert, lat=lat, lon=lon)
        if plot_units:
            if dtype_out_vert == 'vert_int':
                conv = self.var.vert_int_plot_units_conv
            else:
                conv = self.var.plot_units_conv
            data_out = data*conv
        else:
            data_out = data
        # Copy the array to self.data_out for ease of future access.
        self._update_data_out(data_out, dtype_out_time)
        return data_out
            
    # def _load_ann_cycle(self, do_znl_mean=False, region=False):
    #     """Load monthly values Jan-Dec into a single numpy array."""
    #     import numpy as np
    #     from aospy.io import _var_inst
    #     var = _var_inst(var)
    #     # Make year label.
    #     if yr_range == 'default':
    #         try:
    #             yr_range = run.default_yr_range
    #         except AttributeError:
    #             yr_range = model.default_yr_range
    #     # Take zonal mean if quantity is defined zonally.
    #     if do_znl_mean and var.def_lon:
    #         return np.array([load_file(proj, model, run, ens_mem, var, level, i,
    #                                    data_type, yr_range,
    #                                    region=region).mean(axis=-1)
    #                          for i in range(1,13)])
    #     else:
    #         return np.array([load_file(proj, model, run, ens_mem, var, level, i,
    #                                    data_type, yr_range, region=region)
    #                          for i in range(1,13)])

    # def load_plot_data(self, do_znl_mean=False, region=False, lats=False,
    #                    lons=False):
    #     """Load and prep data for plotting."""
    #     import numpy as np
    #     from aospy.io import _time_label, _ens_label, _yr_label, _var_label,
    #     # Load data. Could be single run or list/tuple/etc. of runs.
    #     if type(run) not in [list, tuple]:
    #         if intvl == 'ann_cycle':
    #             data = load_ann_cycle(proj, model, run, ens_mem, var, level,
    #                                   data_type, yr_range, do_znl_mean=False,
    #                                   region=region)
    #         else:
    #             data = load_file(proj, model, run, ens_mem, var, level,
    #                              intvl, data_type, yr_range, region=region,
    #                              lats=lats, lons=lons)
    #     # If multiple runs, use list comprehension to load each one's data.
    #     else:
    #         if intvl == 'ann_cycle':
    #             data = [load_ann_cycle(proj, model, rn, ens_mem, var, level,
    #                                    data_type, yr_range, do_znl_mean=False,
    #                                    region=region)
    #                     for rn in run]
    #         else:
    #             data = [load_file(proj, model, rn, ens_mem, var, level,
    #                               intvl, data_type, yr_range, region=region,
    #                               lats=lats, lons=lons)
    #                     for rn in run]
    #         if len(data) == 1:
    #             data = data[0]
    #         else:
    #             # For data other than stdev, assume last run is the
    #             # control and all preceding are perturbations; take sum of
    #             # (pert minus control) over all perturbation runs.
    #             if 'std' not in data_type:
    #                 data = np.sum(data[:-1], axis=0) - data[-1]*(len(data[:-1]))
    #             # For stdev data, assume runs are independent and
    #             # therefore compute the standard deviation of their sum or
    #             # difference by taking the square root of their squared
    #             # stdev.
    #             else:
    #                 if yr_range == 'default':
    #                     try:
    #                         start_yr, end_yr = run.default_yr_range
    #                     except AttributeError:
    #                         start_yr, end_yr = model.default_yr_range
    #                 data = np.sqrt(np.sum(np.power(data, 2), axis=0) /
    #                                (end_yr - start_yr + 1))
    #     # Convert to plotting units.
    #     try:
    #         # 2015-02-17: Temporary hack to account for different units in this
    #         #             particular dataset.
    #         if var.name in ('precip', 'prec_conv', 'prec_ls', 'p-e', 'evap'):
    #             # CRU output is in mm, not mm/day, so divide by the avg
    #             # month length (365.25/12=30.4375).  This is imperfect:
    #             # need to divide by each month's length.
    #             if model.name == 'cru':
    #                 data /= 30.4375
    #             # Similarly, U. Delaware data is in cm.
    #             elif model.name == 'udel':
    #                 data /= 3.4375
    #             elif model.name not in ('landflux-eval', 'landflux-eval95',
    #                                     'cmap', 'prec_l'):
    #                 data *= var.plot_units_conv
    #         else:
    #             data *= var.plot_units_conv
    #     except AttributeError:
    #         pass
    #     return data


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

    def mask_var(self, data, model):
        """Mask the data of the given variable outside the region."""
        from numpy import tile
        from numpy.ma import masked_where
        # Interpolate region to model grid.
        reg_mask = self.make_mask(model)
        # Mask input values where region mask is zero. Assumes dimensions are
        # (time, level, lat, lon), i.e. data.ndim=4.
        return masked_where(tile(reg_mask == 0.,
                                 (data.shape[0], data.shape[1], 1, 1)), data)

    def ts(self, data, model):
        """Create a time-series of region-average data."""
        import numpy as np
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
        # Tile the weights to be the same shape as the data. Required by the
        # numpy.ma.average function.
        weights = np.tile(weights, (data.shape[0], data.shape[1], 1))
        # Average over the region at each timestep and at each level.
        out = np.squeeze(np.ma.average(data, weights=weights, axis=-1))
        if type(out) is np.ma.core.MaskedArray and out.ndim == 0:
            out = float(out)
        return out

    def av(self, data, model):
        """Time average of region-average data."""
        import numpy as np
        out = np.squeeze(self.ts(data, model).mean(axis=0))
        return out

    def std(self, data, model):
        """Standard deviation of time-series data."""
        import numpy as np
        out = np.squeeze(self.ts(data, model).std(axis=0))
        return out

import calcs, constants, io, plotting, regions, variables, projects
from print_table import print_table
from main import main

        # Average over ensemble members if desired.
        # if self.ens_mem == 'avg' and self.ens_mem_prefix:
        #     ens_avg_dat = []
        #     for i in range(len(self.ens_mem_ext)):
        #         nc = self._get_nc()
        #         nc_vals = nc.variables[self.name][indices]
        #         # Remove unphysical values that can occur in 'av' data.
        #         if self.dtype_in_time == 'av':
        #             nc_vals = self.var.mask_unphysical(nc_vals)
        #         # Get desired vertical level.
        #         if self.level and self.level != 'sigma' and self.def_vert:
        #             ens_avg_dat.append(nc_vals[:,self.lev_ind])
        #         else:
        #             self.ens_avg_dat.append(nc_vals)
        #     vals = np.ma.mean(ens_avg_dat, axis=0)
        # Otherwise just get the values from the desired ensemble member.
        # else:
        # If ensemble averaging was input but the run doesn't have
        # ensemble members, then let the _get_nc() function know this.
        # if self.ens_mem == 'avg':
        #     ens_mem_nc = None
        # else:
        #     ens_mem_nc = ens_mem
