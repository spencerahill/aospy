"""io.py: utility methods used internally by aospy for input/output, etc."""
import os
import string
import subprocess

import numpy as np


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
        return '{:04d}'.format(yr_range[0])
    else:
        return '{:04d}-{:04d}'.format(*yr_range)


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
    elif type(intvl) == int and intvl in range(1, 13):
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
        for lbl, vals in labels.items():
            if intvl == lbl or set(intvl) == set(vals):
                label = lbl
                value = np.array(vals)
                break
    if return_val:
        return label, value
    else:
        return label


def data_in_name_gfdl(name, domain, data_type, intvl_type, data_yr,
                      intvl, data_in_start_yr, data_in_dur):
    """Determines the filename of GFDL model data output."""
    # Determine starting year of netCDF file to be accessed.
    extra_yrs = (data_yr - data_in_start_yr) % data_in_dur
    data_in_yr = data_yr - extra_yrs
    # Determine file name. Two cases: time series (ts) or time-averaged (av).
    if data_type in ('ts', 'inst'):
        if intvl_type == 'annual':
            if data_in_dur == 1:
                filename = '.'.join([domain, '{:04d}'.format(data_in_yr),
                                      name, 'nc'])
            else:
                filename = '.'.join([domain, '{:04d}-{:04d}'.format(
                    data_in_yr, data_in_yr + data_in_dur - 1
                ), name, 'nc'])
        elif intvl_type == 'monthly':
            filename = (domain + '.{:04d}'.format(data_in_yr) + '01-' +
                         '{:04d}'.format(int(data_in_yr+data_in_dur-1)) +
                         '12.' + name + '.nc')
        elif intvl_type == 'daily':
            filename = (domain + '.{:04d}'.format(data_in_yr) + '0101-' +
                         '{:04d}'.format(int(data_in_yr+data_in_dur-1)) +
                         '1231.' + name + '.nc')
        elif 'hr' in intvl_type:
            filename = '.'.join(
                [domain, '{:04d}010100-{:04d}123123'.format(
                    data_in_yr, data_in_yr + data_in_dur - 1), name, 'nc']
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
        if data_in_dur == 1:
            filename = domain + '.{:04d}'.format(data_in_yr) + '.' + label + '.nc'
        else:
            filename = (domain + '.{:04d}'.format(data_in_yr) + '-' +
                         '{:04d}'.format(int(data_in_yr+data_in_dur-1)) +
                         '.' + label + '.nc')
    elif data_type == 'av_ts':
        filename = (domain + '.{:04d}'.format(data_in_yr) + '-' +
                     '{:04d}'.format(int(data_in_yr+data_in_dur-1)) + '.01-12.nc')
    return filename


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


def get_data_in_direc_repo(data_in_direc, var_name, version=-1):
    """Determine the directory containing the needed netCDF files."""
    dir_prefix = os.path.realpath(data_in_direc)
    dirs_after = [os.path.join(dir_prefix, name)
                  for name in os.listdir(dir_prefix) if
                  (os.path.isdir(os.path.join(dir_prefix, name)) and
                   name[0] == 'v')]
    dirs_after.sort()
    # We assume that the directories are all of the format "v#[###...]",
    # the numbers typically specifiying a date YYYYMMDD, but occasionally
    # just a single number, e.g. 'v1' instead.  We assume we want the
    # latest/most recent version, which after sorting will be the last
    # entry in the list.  Or the version can be specified using the
    # 'version' kwarg.
    # Then append the run's direc_nc specification.
    return os.path.join(dir_prefix, dirs_after[version], var_name)
