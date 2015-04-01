"""
io submodule of aospy module.
"""
def _var_inst(var):
    """Convert string of an aospy.var name to an aospy.var instance."""
    import aospy.variables
    if type(var) is str:
        try:
            var_inst = getattr(aospy.variables, var)
        except AttributeError:
            raise AttributeError
    else:
        var_inst = var
    return var_inst

def _region_inst(region):
    """Convert string of an aospy.Region name to an aospy.Region instance."""
    import aospy.regions
    if type(region) is str:
        return getattr(aospy.regions, region)
    else:
        return region

def _proj_inst(proj):
    """Convert string of an aospy.Proj name to an aospy.Proj instance."""
    import aospy.projects
    if type(proj) is aospy.Proj:
        return proj
    elif type(proj) is str:
        try:
            return getattr(aospy.projects, proj)()
        except AttributeError:
            raise AttributeError('Not a recognized aospy.Proj instance: %s'
                                 % proj)

def _model_inst(model, parent_proj):
    """Convert string of an aospy.model name to an aospy.model instance."""
    from aospy.io import _proj_inst
    if type(model) is str:
        try:
            parent_proj = _proj_inst(parent_proj)
        except:
            pass
        finally:
            model_inst = parent_proj.models[model]
    else:
        model_inst =  model
    # Set the Model-Proj child-parent relation.
    model_inst.proj = parent_proj
    return model_inst

def _run_inst(run, parent_model, parent_proj):
    """Convert string matching an aospy.run name to an aospy.run instance."""
    from aospy.io import _model_inst
    if type(run) is str:
        if type(parent_model) is str:
            parent_model = _model_inst(parent_model, parent_proj)
        run_inst = parent_proj.models[parent_model.name].runs[run]
    elif type(run) in (list, tuple):
        parent_model = _model_inst(parent_model, parent_proj)
        run_inst = [parent_proj.models[parent_model.name].runs[rn]
                    for rn in run]
    else:
        run_inst = run
    return run_inst

def _aospy_inst(proj=False, model=False, run=False, var=False):
    """Convert string matching aospy object names to class instances."""
    from aospy.io import (_run_inst, _region_inst, _proj_inst, _model_inst,
                          _var_inst)
    def _to_list(obj):
        if type(obj) is str:
            list_obj = [obj]
        else:
            list_obj = obj
        return list_obj

    pr, md, rn, vr = [], [], [], []
    for p, m, r, v in zip(_to_list(proj), _to_list(model), _to_list(run),
                          _to_list(var)):
        if p:
            p = _proj_inst(p)
            pr.append(p)
        if m:
            m = _model_inst(m, p)
            md.append(m)
        if r:
            r = _run_inst(r, m, p)
            rn.append(r)
        if v:
            vr.append(_var_inst(v))

    if type(proj) is str or len(proj) == 1:
        pr = pr[0]
    if type(model) is str or len(model) == 1:
        md = md[0]
    if type(run) is str or len(run) == 1:
        rn = rn[0]
    if type(var) is str or len(var) == 1:
        vr = vr[0]

    return pr, md, rn, vr

def _ens_label(ens_mem):
    """Create label of the ensemble member for aospy data I/O."""
    if ens_mem is None:
        return ''
    elif ens_mem == 'avg':
        return 'ens_mean'
    else:
        return 'mem' + str(ens_mem + 1)

def _yr_label(yr_range):
    """Create label of start and end years for aospy data I/O."""
    assert yr_range is not None, "yr_range is None" 
    if yr_range[0] == yr_range[1]:
        return '{:04}'.format(yr_range[0])
    else:
        return '{:04}'.format(yr_range[0]) + '-' + '{:04}'.format(yr_range[1])

def _var_label(var, level):
    """
    Create label of variable name and potentially the desired level for aospy
    data I/O.
    """
    if type(var) == str:
        var_name = var
        defvert = False
    else:
        var_name = var.name
        defvert = var.def_vert
    if type(level) is str and defvert:
        return var_name + '.' + str(level)
    else:
        return var_name

def _znl_label(var):
    """Create label denoting zonal mean values for aospy data I/O."""
    try: 
        if var.def_lon:
            return 'znl'
        else:
            return ''
    except:
        return ''

def _time_label(intvl):
    """Create time interval label for aospy data I/O."""
    import numpy as np
    # Monthly labels are 2 digit integers: '01' for jan, '02' for feb, etc.
    if type(intvl) in [list, tuple, np.ndarray] and len(intvl) == 1:
        return '{:02}'.format(intvl[0]), np.array(intvl)
    elif type(intvl) == int and intvl in range(1,13):
        return '{:02}'.format(intvl), np.array([intvl])
    # Seasonal and annual time labels are short strings.
    else:
        labels = {'jfm': (1, 2, 3), 'fma': (2, 3, 4), 'mam': (3,  4,  5),
                  'amj': (4, 5, 6), 'mjj': (5, 6, 7), 'jja': (6,  7,  8),
                  'jas': (7, 8, 9), 'aso': (8, 9,10), 'son': (9, 10, 11),
                  'ond':(10,11,12), 'ndj': (11,12,1), 'djf': (1,  2, 12),
                  'jjas': (6,7,8,9), 'djfm': (12, 1, 2, 3), 
                  'ann': range(1,13)}
        for lbl, vals in labels.iteritems():
            if intvl == lbl or set(intvl) == set(vals):
                return lbl, np.array(vals)

def _intvl_indices_and_label(intvl, intvl_type):
    """Get the indices of the given time label and its string label."""
    import numpy as np
    from aospy.io import _time_label
    # Annual mean input data has only one index per year and the label 'ann'.
    # For all other input data, look up the correct label and indices for the
    # given interval.
    if intvl_type == 'annual':
        label, intvl_indices = 'ann', np.array([0])
    elif intvl_type == 'seasonal':
        label, intvl_indices = _time_label(intvl)
        intvl_indices = np.array([0])
    else:
        label, intvl_indices = _time_label(intvl)
    return label, intvl_indices

def save_file(data, proj, model, run, ens_mem, var, 
              level, intvl, out_dtype, yr_range):
    """Save aospy data to an external file."""
    from os.path import isdir
    from os import makedirs
    from cPickle import dump, load
    import tarfile
    # Create necessary string labels.
    ens_label = _ens_label(ens_mem)
    start_yr, end_yr, num_yr  = var.data_yr_range(yr_range, 'ts')
    yr_label = _yr_label((start_yr, end_yr))
    time_label, val = _time_label(intvl) 
    if time_label == '00':
        time_label = 'ann'
    path_out = (proj.direc_out + '/data/' + model.name + '/' +
                run.name + '/' + ens_label + '/').replace('//', '/')
    path_scratch = ('/work/s1h/' + proj.name + '/').replace('//', '/')
    file_name = (_var_label(var, level) + '.' + out_dtype + '.' + model.name +
                 '.' + run.name + '.' + ens_label + '.' + time_label +
                 '.' + yr_label + '.p').replace('..','.')
    # Save the data to the scratch filesystem.
    with open(path_scratch + file_name, 'a+') as file_work:
        # Update, rather than overwrite, existing regional data.
        if 'reg' in  out_dtype:
            # Open the existing dictionary if it exists.
            try:
                reg_data = load(file_work)
            except EOFError:
                reg_data = {}
            # Add the new data to the dictionary.
            reg_data.update(data)
            data_out = reg_data
        else:
            data_out = data
    # Export the data to the external file.
    with open(path_scratch + file_name, 'w') as file_work:
        dump(data_out, file_work)
    # Add the data to the tar file in /archive.
    if not isdir(path_out):
        makedirs(path_out)
    with tarfile.open(path_out+'data.tar', 'a') as tar:
        tar.add(path_scratch+file_name, arcname=path_out+file_name)
    print file_name

def data_av_stat(loc_ts, calcs_out):
    """Compute all desired calculations on a local time-series."""
    files = {}
    if 'ts' in calcs_out:
        files.update({'ts': loc_ts})
    if 'av' in calcs_out:
        files.update({'av': loc_ts.mean(axis=0)})
    if 'std' in calcs_out:
        files.update({'std': loc_ts.std(axis=0)})
    # Zonal mean.
    if any('znl' in out_type for out_type in calcs_out):
        znl_ts = loc_ts.mean(axis=-1)
        if 'znl.ts' in calcs_out:
            files.update({'znl.ts': znl_ts})
        if 'znl.av' in calcs_out:
            files.update({'znl.av': znl_ts.mean(axis=0)})
        if 'znl.std' in calcs_out:
            files.update({'znl.std': znl_ts.std(axis=0)})
    # Zonal asymmetry.
    if any('zasym' in out_type for out_type in calcs_out):
        # '.T'=transpose operator; just to make numpy broadcasting syntax work.
        zasym_ts = (loc_ts.T - znl_ts.T).T
        if 'zasym.ts' in calcs_out:
            files.update({'zasym.ts': zasym_ts})
        if 'zasym.av' in calcs_out:
            files.update({'zasym.av': zasym_ts.mean(axis=0)})
        if 'zasym.std' in calcs_out:
            files.update({'zasym.std': zasym_ts.std(axis=0)})
    return files

def regions_calcs(regions, calcs_out, loc_ts, proj, model, run, ens_mem, var, 
                  level, intvl, yr_range, label):
    """Region-averaged computations. Execute and save to external file."""
    from aospy.io import _region_inst
    calcs_reg = ['ts', 'av', 'std', 'spvar.ts', 'spvar.av', 'znl.spvar.ts',
                 'znl.spvar.av', 'zasym.ts', 'zasym.av', 'zasym.std',
                 'zasym.spvar.ts', 'zasym.spvar.av']
    regions = [_region_inst(region) for region in regions]
    # Perform each calculation for each region.
    for calc in calcs_reg:
        if 'reg.' + calc in calcs_out:
            reg_dat = {}
            for region in regions:
                method = region.__getattribute__(calc)
                reg_dat.update({region.name: method(loc_ts, var, model, label)})
            # Save the data to external file.
            save_file(reg_dat, proj, model, run, ens_mem, var, 
                      level, intvl, 'reg.' + calc, yr_range)

def prune_lat_lon(array, model, lats, lons):
    """Cut-off data outside desired lat-lon range."""
    # SAH 2015-03-10: have to pivot lons that span edge of array. See
    # "pivot" portion of aospy.plotting.plot_lon_vert.  Need to implement that
    # before this function will work.
    import numpy as np
    # Assume array dimensions are (time, lat, lon).
    if lats:
        latmin, latmax = np.min(lats), np.max(lats)
        lats_ind = np.where((model.lat > latmin) & (model.lat < latmax))
        array = array[:,lats_ind]
    if lons:
        lonmin, lonmax = np.min(lons), np.max(lons)
        lons_ind = np.where((model.lon > lonmin) & (model.lon < lonmax))
        array[:,:,lons_ind]
    return array

def load_file(proj, model, run, ens_mem, var, level,
              intvl, data_type, yr_range, region=False,
              lats=False, lons=False):
    """Load aospy data saved externally."""
    from cPickle import load
    import numpy as np
    from aospy.io import (_time_label, _ens_label, _yr_label, _var_label,
                          _var_inst, _proj_inst, _model_inst, _run_inst)
    # Make filename and filepath strings.
    time_label, index = _time_label(intvl)
    ens_label = _ens_label(ens_mem)
    # Get aospy objects if only the name was passed.
    proj = _proj_inst(proj)
    model = _model_inst(model, proj)
    run = _run_inst(run, model, proj)
    var = _var_inst(var)
    # Make year label.
    if yr_range == 'default':
        try:
            yr_range = run.default_yr_range
        except AttributeError:
            yr_range = model.default_yr_range    
    yr_label = _yr_label(yr_range)
    if var.def_vert:
        var_label = var.name
        # Get specified level, or all levels if none specified.
        if np.max(model.level) > 1e4:
            # Convert from Pa to hPa.
            lev_hpa = model.level*1e-2
        else:
            lev_hpa = model.level
        if level is None:
            level_index = np.where(lev_hpa > 0.)
        else:
            level_index = np.where(lev_hpa == level)
    else:
        var_label = _var_label(var, level)
    name = (var_label + '.' + data_type + 
            '.' + model.name + '.' + run.name + '.' + ens_label + '.' + 
            time_label + '.' + yr_label + '.p').replace('..', '.')
    # Get the data from /work if possible, otherwise from /archive.
    # try:
    with open('/work/s1h/' + proj.name + '/' +  name, 'r') as data:
        data_vals = load(data)
        if 'reg' in data_type and region:
            data_vals = data_vals[region]
        if level and var.def_vert:
            if 'ts' in data_type:
                data_vals = np.squeeze(data_vals[:,level_index])
            else:
                data_vals = np.squeeze(data_vals[level_index])                
        # if lats or lons:
            # data_vals = prune_lat_lon(data_vals, model, lats, lons)
        return data_vals
    # except:
    #     # Have to extract from '/archive/.../data.tar'.
    #     from subprocess import call
    #     import tarfile
    #     path_in = (proj.direc_out + '/data/' + model.name + '/' +
    #                run.name + '/' + ens_label + '/').replace('//', '/')
    #     call(['dmget'] + [path_in + 'data.tar'])
    #     with tarfile.open(path_in + 'data.tar', 'r') as data_tar:
    #         data_vals = load(data_tar.extractfile(name))
    #         if 'reg' in data_type and region:
    #             return data_vals[region]
    #         else:
    #             return data_vals

def load_ann_cycle(proj, model, run, ens_mem, var, level, 
                   data_type, yr_range, do_znl_mean=False, region=False):
    """Load monthly values Jan-Dec into a single numpy array."""
    import numpy as np
    from aospy.io import _var_inst
    var = _var_inst(var)
    # Make year label.
    if yr_range == 'default':
        try:
            yr_range = run.default_yr_range
        except AttributeError:
            yr_range = model.default_yr_range
    # Take zonal mean if quantity is defined zonally.
    if do_znl_mean and var.def_lon:
        return np.array([load_file(proj, model, run, ens_mem, var, level, i, 
                                   data_type, yr_range, 
                                   region=region).mean(axis=-1)
                         for i in range(1,13)])
    else:
        return np.array([load_file(proj, model, run, ens_mem, var, level, i, 
                                   data_type, yr_range, region=region) 
                         for i in range(1,13)])

def load_plot_data(proj, model, run, ens_mem, var, level, 
                   intvl, data_type, yr_range, 
                   do_znl_mean=False, region=False, lats=False, lons=False):
    """Load and prep data for plotting."""
    import numpy as np
    from aospy.io import (_time_label, _ens_label, _yr_label, _var_label,
                          _var_inst, _proj_inst, _model_inst, _run_inst)
    # Get aospy objects if only the name was passed.
    proj = _proj_inst(proj)
    model = _model_inst(model, proj)
    run = _run_inst(run, model, proj)
    var = _var_inst(var)
    # Load data. Could be single run or list/tuple/etc. of runs.
    if type(run) not in [list, tuple]:
        if intvl == 'ann_cycle':
            data = load_ann_cycle(proj, model, run, ens_mem, var, level, 
                                  data_type, yr_range, do_znl_mean=False, 
                                  region=region)
        else:
            data = load_file(proj, model, run, ens_mem, var, level, 
                             intvl, data_type, yr_range, region=region,
                             lats=lats, lons=lons)
    # If multiple runs, use list comprehension to load each one's data.
    else:
        if intvl == 'ann_cycle':
            data = [load_ann_cycle(proj, model, rn, ens_mem, var, level, 
                                   data_type, yr_range, do_znl_mean=False, 
                                   region=region) 
                    for rn in run]
        else:
            data = [load_file(proj, model, rn, ens_mem, var, level, 
                              intvl, data_type, yr_range, region=region,
                              lats=lats, lons=lons) 
                    for rn in run]
	if len(data) == 1:
            data = data[0]
	else:
            # For data other than stdev, assume last run is the
            # control and all preceding are perturbations; take sum of
            # (pert minus control) over all perturbation runs.
            if 'std' not in data_type:
                data = np.sum(data[:-1], axis=0) - data[-1]*(len(data[:-1]))
            # For stdev data, assume runs are independent and
            # therefore compute the standard deviation of their sum or
            # difference by taking the square root of their squared
            # stdev.
            else:
                if yr_range == 'default':
                    try:
                        start_yr, end_yr = run.default_yr_range
                    except AttributeError:
                        start_yr, end_yr = model.default_yr_range
                data = np.sqrt(np.sum(np.power(data, 2), axis=0) / 
                               (end_yr - start_yr + 1))
    # Convert to plotting units.
    try:
        # 2015-02-17: Temporary hack to account for different units in this
        #             particular dataset.
        if var.name in ('precip', 'prec_conv', 'prec_ls', 'p-e', 'evap'):
            # CRU output is in mm, not mm/day, so divide by the avg
            # month length (365.25/12=30.4375).  This is imperfect:
            # need to divide by each month's length.
            if model.name == 'cru':
                data /= 30.4375
            # Similarly, U. Delaware data is in cm.
            elif model.name == 'udel':
                data /= 3.4375
            elif model.name not in ('landflux-eval', 'landflux-eval95',
                                    'cmap', 'prec_l'): 
                data *= var.plot_units_conv
        else:
            data *= var.plot_units_conv
    except AttributeError:
        pass
    return data

def get_timesteps(intvl, start_yr, end_yr, nc_start_yr, nc_end_yr, nc_dur, 
                  intvl_type, data_type, nc_start_month=1, nc_end_month=12):
    """Determine desired time indices."""
    import numpy as np
    if intvl_type == 'annual':
        label, intvl = 'ann', [0]
    else:
        label, intvl = _time_label(intvl)
    # Given time range must be subset of available time range.
    assert start_yr >= nc_start_yr
    assert end_yr <= nc_end_yr
    if start_yr == nc_start_yr:
        assert np.min(intvl) >= nc_start_month
    if end_yr == nc_end_yr:
        assert np.max(intvl) <= nc_end_month
    # Number of netCDF files preceding the first needed netCDF file.
    num_prior_nc = (start_yr - nc_start_yr) / nc_dur
    # Number of years preceding the first needed netCDF file.
    yrs_prior = num_prior_nc * nc_dur
    # Start year of the first necessary netCDF file.
    first_nec_nc_st_yr = nc_start_yr + yrs_prior
    # Number of years into the first needed netCDF file that data starts.
    yrs_into_first_nc = start_yr - first_nec_nc_st_yr
    num_indices = (end_yr + 1 - start_yr)*12
    # ts data has a value for each month or year for each year.
    if data_type == 'ts':
        if intvl_type == 'monthly':
            indices = []
            for i in intvl:
                start_index = yrs_into_first_nc*12 + (i - nc_start_month)
                end_index = start_index + num_indices
                for t in range(start_index, end_index, 12):
                    indices.append(t)
            indices.sort()
            return indices
        # Only one timestep per year for annual ts data.
        elif intvl_type == 'annual':
            return np.arange(start_yr - nc_start_yr - nc_dur*prior_nc, 
                             end_yr - nc_start_yr - nc_dur*prior_nc + 1)
    # av data has one time averaged value for the designated time interval.
    elif data_type == 'av':
        return np.array([0])
    # av_ts data has one time averaged value for each month.
    elif data_type == 'av_ts':
        # January is month 1 but index 0, and so on, so subtract one.
        return np.array(intvl) - 1

def nc_name_gfdl(name, domain, data_type, intvl_type, data_yr,
                 intvl, nc_start_yr, nc_dur):
    """Determines the gfdl_file name of GFDL model data output."""
    from aospy.io import _time_label, _intvl_indices_and_label
    # Determine starting year of netCDF file to be accessed.
    extra_yrs = (data_yr - nc_start_yr) % nc_dur
    nc_yr = data_yr - extra_yrs
    # Determine file name. Two cases: time series (ts) or time-averaged (av).
    if data_type == 'ts':
        if intvl_type == 'annual':
            if nc_dur == 1:
                gfdl_file = (domain + '.{:04}'.format(nc_yr) +
                        '.' + name + '.nc')
            else:
                gfdl_file = (domain + '.{:04}'.format(nc_yr) +
                        '-{:04}'.format(nc_yr+nc_dur-1)
                        + '.' + name + '.nc')
        elif intvl_type == 'monthly':
            gfdl_file = (domain + '.{:04}'.format(nc_yr) + '01-' +
                    '{:04}'.format(int(nc_yr+nc_dur-1)) +
                    '12.' + name + '.nc')
    elif data_type == 'av':
        if intvl_type in ['annual', 'ann']:
            label = 'ann'
        elif intvl_type in ['seasonal', 'seas']:
            label, val = _intvl_indices_and_label(intvl, intvl_type)
            label = label.upper()
        elif intvl_type in ['monthly', 'mon']:
            label, val = _time_label(intvl)
        if nc_dur == 1:
            gfdl_file = domain + '.{:04}'.format(nc_yr) + '.' + label +'.nc'
        else:
            gfdl_file = (domain + '.{:04}'.format(nc_yr) + '-' +
                    '{:04}'.format(int(nc_yr+nc_dur-1)) +
                    '.' + label + '.nc')
    elif data_type == 'av_ts':
        gfdl_file = (domain + '.{:04}'.format(nc_yr) + '-' +
                '{:04}'.format(int(nc_yr+nc_dur-1)) + '.01-12.nc')
    return gfdl_file
 
def dmget_nc(files_list):
    """Call GFDL command 'dmget' to access archived files."""
    from subprocess import call
    call(['dmget'] + files_list)

def hsmget_nc(files_list):
    """Call GFDL command 'hsmget' to access archived files."""
    import os, subprocess, string
    tmpdir = os.environ['TMPDIR']
    workdir = '/work/' + os.environ['USER']
    # Assumes files are located somewhere on /archive.
    # Assumes files_list is list of absolute paths to the netCDF files.
    # Assumes that all files in files list are under same archive root.
    arch_loc = string.join(files_list[0].split('/')[:3],'/') + '/'
    files = [f.partition(arch_loc)[2] for f in files_list]
    # subprocess.call(['module', 'load', 'hsm'])
    retcode = subprocess.call(['hsmget', '-a', arch_loc, '-p', workdir,
                               '-w', tmpdir] + files + ['-q'])
    return retcode

def _get_time_avg_var(proj, model_name, run_name, var):
    """Get the desired Var object and designate its parent Run oject."""
    try:
        var = proj.vars[var.name]
    except KeyError:
        for alt_name in self.alt_names:
            try:
                var = proj.vars[alt_name]
            except KeyError:
                pass
            else:
                vals = vals[indices]
                break
    except AttributeError:
        var = proj.vars[var]
    var.proj = proj
    try:
        var.model = proj.models[model_name]
    except KeyError:
        var.model = model_name
    # try/except TypeError clause is hack solution to deal with plotting
    # functions inputting a list of run names rather than a single run name
    # string.
    try:
        var.run = var.model.runs[run_name]
    except KeyError:
        var.run = run_name
    except TypeError:
        var.run = var.model.runs[run_name[-1]]
    return var
    
def mult_var_calcs(proj, model, run, ens_mem, calc, level, intvl, 
                   intvl_type, dtype, yr_range, **kwargs):
    """Perform calculations involving one or more variables and save."""
    import numpy as np
    from aospy import _get_parent_attr
    # Get aospy objects if only the name was passed.
    if type(calc) == tuple:
        cal, params = calc
    else:
        cal = calc
        params = {}
    cal = _get_time_avg_var(proj, model, run, cal)
    # Copy name strings of aospy objects
    model_name, run_name = model, run
    # Get start and end years and 'effective' number of year timesteps.
    start_yr, end_yr, num_yr = cal.data_yr_range(yr_range, dtype)
    # Get time interval.
    label, intvl = _intvl_indices_and_label(intvl, intvl_type)
    # Load data.
    vars_calc = []
    for n, var in enumerate(cal.vars):
        var = _get_time_avg_var(proj, model, run, var)
        # If specified in params, use variable specific model/run/ens_mem.
        if params.get('models'):
            var_model = params['models'][n]
        else:
            var_model = proj.models[model_name]
        if params.get('runs'):
            var_run = params['runs'][n]
        else:
            var_run = var_model.runs[run_name]
        if params.get('ens_mems'):
            var_ens_mem = params['ens_mems'][n]
        else:
            var_ens_mem = ens_mem
        # For time-independent dimension variables (e.g. latitude), tile the
        # array so that it can be accessed in list comprehension below using
        # same syntax as time-dependent variables.
        if var.name in ['lat', 'lon', 'level', 'land_mask', 'sfc_area']:
            grid_val = var_model.__getattribute__(var.name)
            if dtype == 'ts':
                if var.name in ['lat', 'lon', 'level']:
                    vars_calc.append(np.tile(grid_val, (num_yr*len(intvl), 1)))
                elif var.name in ['land_mask', 'sfc_area']:
                    vars_calc.append(
                        np.tile(grid_val, (num_yr*len(intvl), 1, 1))
                    )
            elif dtype == 'av':
                vars_calc.append(grid_val)
        # For time-dependent vars, get values at the necessary timesteps.
        else:
            # Get time indices of the input data corresponding to the
            # desired timesteps.
            nc_start_yr = _get_parent_attr(var, 'nc_start_yr')
            nc_end_yr = _get_parent_attr(var, 'nc_end_yr')
            nc_dur = _get_parent_attr(var, 'nc_dur')
            nc_start_month = _get_parent_attr(var, 'nc_start_month')
            if nc_start_month is None:
                nc_start_month = 1
            nc_end_month = _get_parent_attr(var, 'nc_end_month')
            if nc_end_month is None:
                nc_end_month = 12
            indices = get_timesteps(
                intvl, start_yr, end_yr, nc_start_yr, nc_end_yr, nc_dur,
                intvl_type, dtype, nc_start_month=nc_start_month,
                nc_end_month=nc_end_month
            )
            # If ensemble mean, average over all members.
            if var_ens_mem == 'avg':
                ens_avg_data = []
                for i in range(len(var_run.ens_mem_ext)):
                    vals, dt = var._get_data_vals(
                        i, level, intvl, intvl_type,
                        dtype, start_yr, end_yr, indices
                    )
                    ens_avg_data.append(vals)
                vals = np.average(ens_avg_data, axis=0)
            # Otherwise just get the values from the desired ensemble member.
            else:
                vals, dt = var._get_data_vals(
                    var_ens_mem, level, intvl, intvl_type,
                    dtype, start_yr, end_yr, indices
                )
            vars_calc.append(vals)
    # Execute calculation at each time index.
    if dtype == 'ts':
        vals = np.ma.array([cal.func(*[var_calc[t] for var_calc in vars_calc])
                         for t in range(len(indices))])
        # For broadcasting purposes, add empty dimension if needed.
        if not cal.def_vert:
            vals = vals[:,np.newaxis]
    else:
        vals = cal.func(*vars_calc)[np.newaxis,:]
    # Split time dimension into two: one for each year, and one for each
    # timestep within a year.
    dt_ts = np.reshape(dt, (num_yr, intvl.size))
    dt_ts = dt_ts[:,:,np.newaxis,np.newaxis,np.newaxis]
    # # Reshape the data array using these new time dimensions.
    intvl_ts = vals.reshape((num_yr, intvl.size, vals.shape[-3],
                             vals.shape[-2], vals.shape[-1]))
    # Average over all timesteps within each year to get an annual time-series.
    loc_ts = np.ma.multiply(intvl_ts, dt_ts).sum(axis=1) / dt_ts.sum(axis=1)    
    # Compute and save all the desired computations on the data.
    calcs_out = kwargs.get('calcs_out', ['av'])
    files = data_av_stat(loc_ts, calcs_out)
    # Save.
    if params.get('save_model'):
        save_model = _model_inst(params['save_model'], proj)
    else:
        save_model = _model_inst(model, proj)
    if params.get('save_run'):
        save_run = _run_inst(params['save_run'], save_model, proj)
    else:
        save_run = _run_inst(run, save_model, proj)
    if params.get('save_ens_mem'):
        save_ens_mem = params['save_ens_mem']
    else:
        save_ens_mem = ens_mem
    for name, data in files.iteritems():
        save_file(data, proj, save_model, save_run, save_ens_mem, cal, 
                  level, intvl, name, yr_range)
    # Region averages and other computations.
    regions = kwargs.get('regions', ['globe'])
    regions_calcs(regions, calcs_out, loc_ts, proj, proj.models[model],
                  proj.models[model].runs[run], 
                  ens_mem, cal, level, intvl, yr_range, label)
