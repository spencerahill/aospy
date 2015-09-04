"""io.py: utility methods used internally by aospy for input/output, etc."""
import os
import string
import subprocess
import numpy as np
import netCDF4

def to_dup_list(x, n, single_to_list=True):
    """
    Convert singleton or iterable into length-n list.  If the input is
    a list, with length-1, its lone value gets duplicated n times.  If
    the input is a list with length-n, leave it the same.  If the
    input is any other data type, replicate it as a length-n list.
    """
    if isinstance(x, list):
        if len(x) == n:
            return x
        elif len(x) == 1:
            return x*n
        else:
            raise ValueError("Input %s must have length 1 or %d : len(%s) = %d"
                             % (x, n, x, len(x)))
    else:
        if n == 1:
            if single_to_list:
                return [x]
            else:
                return x
        else:
            return [x]*n


def _var_label(var, level):
    """
    Create label of variable name and potentially the desired level for aospy
    data I/O.
    """
    if isinstance(var, str):
        var_name = var
        defvert = False
    else:
        var_name = var.name
        defvert = var.def_vert
    if (defvert and level is not None) or level == 'sigma':
        return var_name + '.' + str(level)
    else:
        return var_name


def _data_in_label(intvl_in, dtype_in_time, dtype_in_vert=False):
    """Create string label specifying the input data of a calculation."""
    intvl_lbl = intvl_in
    time_lbl = dtype_in_time
    lbl = '_'.join(['from', intvl_lbl, time_lbl]).replace('__', '_')
    vert_lbl = dtype_in_vert if dtype_in_vert else False
    if vert_lbl:
        lbl = '_'.join([lbl, vert_lbl]).replace('__', '_')
    return lbl


def _data_out_label(time_intvl, dtype_time, dtype_vert=False):
    intvl_lbl = _time_label(time_intvl, return_val=False)
    time_lbl = dtype_time
    lbl = '.'.join([intvl_lbl, time_lbl]).replace('..', '.')
    vert_lbl = dtype_vert if dtype_vert else False
    if vert_lbl:
        lbl = '.'.join([lbl, vert_lbl]).replace('..', '.')
    return lbl


def _ens_label(ens_mem):
    """Create label of the ensemble member for aospy data I/O."""
    if ens_mem in (None, False):
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


def _znl_label(var):
    """Create label denoting zonal mean values for aospy data I/O."""
    try:
        if var.def_lon:
            return 'znl'
        else:
            return ''
    except:
        return ''


def _time_label(intvl, return_val=True):
    """Create time interval label for aospy data I/O."""
    # Monthly labels are 2 digit integers: '01' for jan, '02' for feb, etc.
    if type(intvl) in [list, tuple, np.ndarray] and len(intvl) == 1:
        label = '{:02}'.format(intvl[0])
        value = np.array(intvl)
    elif type(intvl) == int and intvl in range(1,13):
        label = '{:02}'.format(intvl)
        value = np.array([intvl])
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
                label = lbl
                value = np.array(vals)
                break
    if return_val:
        return label, value
    else:
        return label


def _month_indices(months, iterable=True):
    """Convert string labels for months to integer indices.

    :param months: String matching either 'ann' or some subset of
                   'jfmamjjasond'.  If 'ann', use all months.  Otherwise, use
                   the specified months.
    :type months: str or int
    """
    assert type(months) in (int, str)
    if type(months) is int:
        if iterable:
            return [months]
        else:
            return months
    elif months.lower() == 'ann':
        return range(1,13)
    else:
        first_letter = 'jfmamjjasond'*2
        # Python native indexing starts at 0, but within aospy months are
        # indexed starting at 1, so add 1.
        st_ind = first_letter.find(months.lower()) + 1
        return range(st_ind, st_ind + len(months))


def _get_time(time, units, calendar, start_yr, end_yr, months, indices=False):
    """Determine the indices of a time array falling in a specified interval.

    Given a start year, end year, and subset of each year, determine which of
    the input time array's values fall within that range.

    :param time: netCDF4 variable object specifying time
    :param start_yr, end_yr: Start and end years, inclusive, of desired time
                             range.
    :type start_yr, end_yr: int
    :param months: Subset of the annual cycle to sample.
    :type months: Iterable of ints in the range (1,13), inclusive.
    :param indices: Return an array of indices if True, otherwise return
                    the time array itself at those time indices.
    :type indices: bool
    """
    dates = netCDF4.num2date(time[:], units, calendar.lower())
    inds = [i for i, date in enumerate(dates) if (date.month in months) and
            (date.year in range(start_yr, end_yr+1))]
    if indices == 'only':
        return inds
    elif indices:
        return inds, time[inds]
    else:
        return time[inds]


def prune_lat_lon(array, model, lats, lons):
    """Cut-off data outside desired lat-lon range."""
    # SAH 2015-03-10: have to pivot lons that span edge of array. See
    # "pivot" portion of aospy.plotting.plot_lon_vert.  Need to implement that
    # before this function will work.
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


def nc_name_gfdl(name, domain, data_type, intvl_type, data_yr,
                 intvl, nc_start_yr, nc_dur):
    """Determines the gfdl_file name of GFDL model data output."""
    # Determine starting year of netCDF file to be accessed.
    extra_yrs = (data_yr - nc_start_yr) % nc_dur
    nc_yr = data_yr - extra_yrs
    # Determine file name. Two cases: time series (ts) or time-averaged (av).
    if data_type in ('ts', 'inst'):
        if intvl_type == 'annual':
            if nc_dur == 1:
                gfdl_file = '.'.join([domain, '{:04}'.format(nc_yr),
                                      name, 'nc'])
            else:
                gfdl_file = (domain + '.{:04}'.format(nc_yr) +
                             '-{:04}'.format(nc_yr+nc_dur-1)
                             + '.' + name + '.nc')
        elif intvl_type == 'monthly':
            gfdl_file = (domain + '.{:04}'.format(nc_yr) + '01-' +
                         '{:04}'.format(int(nc_yr+nc_dur-1)) +
                         '12.' + name + '.nc')
        elif intvl_type == 'daily':
            gfdl_file = (domain + '.{:04}'.format(nc_yr) + '0101-' +
                         '{:04}'.format(int(nc_yr+nc_dur-1)) +
                         '1231.' + name + '.nc')
        elif 'hr' in intvl_type:
            gfdl_file = '.'.join(
                [domain, '{:04}'.format(nc_yr) + '010100-' +
                 '{:04}'.format(nc_yr+nc_dur-1) + '123123', name, 'nc']
            )
    elif data_type == 'av':
        if intvl_type in ['annual', 'ann']:
            label = 'ann'
        elif intvl_type in ['seasonal', 'seas']:
            # 2015-04-29: This function is obsolete.  Should fix this
            # eventually.  But I almost never use seasonal data, so not a
            # priority.
            # label, val = _intvl_indices_and_label(intvl, intvl_type)
            label = label.upper()
        elif intvl_type in ['monthly', 'mon']:
            label, val = _time_label(intvl)
        if nc_dur == 1:
            gfdl_file = domain + '.{:04}'.format(nc_yr) + '.' + label + '.nc'
        else:
            gfdl_file = (domain + '.{:04}'.format(nc_yr) + '-' +
                         '{:04}'.format(int(nc_yr+nc_dur-1)) +
                         '.' + label + '.nc')
    elif data_type == 'av_ts':
        gfdl_file = (domain + '.{:04}'.format(nc_yr) + '-' +
                     '{:04}'.format(int(nc_yr+nc_dur-1)) + '.01-12.nc')
    return gfdl_file


def dmget(files_list):
    """Call GFDL command 'dmget' to access archived files."""
    subprocess.call(['dmget'] + files_list)


def hsmget_nc(files_list):
    """Call GFDL command 'hsmget' to access archived files."""
    # tmpdir = os.environ['TMPDIR']
    workdir = '/work/' + os.environ['USER']
    ptmpdir = '/ptmp/' + os.environ['USER']
    # Assumes files are located somewhere on /archive.
    # Assumes files_list is list of absolute paths to the netCDF files.
    # Assumes that all files in files list are under same archive root.
    arch_loc = string.join(files_list[0].split('/')[:3], '/') + '/'
    files = [f.partition(arch_loc)[2] for f in files_list]
    # subprocess.call(['module', 'load', 'hsm'])
    retcode = subprocess.call(['hsmget', '-a', arch_loc, '-p', workdir,
                               '-w', ptmpdir] + files + ['-q'])
    return retcode


def _get_time_avg_var(proj, model_name, run_name, var):
    """Get the desired Var object and designate its parent Run oject."""
    try:
        var = proj.vars[var.name]
    except KeyError:
        for alt_name in alt_names:
            try:
                var = proj.vars[alt_name]
            except KeyError:
                pass
            else:

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
