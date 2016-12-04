"""io.py: utility methods used internally by aospy for input/output, etc."""
import imp
import os
import logging
import string
import subprocess

import numpy as np

from ..__config__ import user_path

logging.basicConfig(level=logging.INFO)


def load_user_data(name):
    """Load user data from aospy_path for given module name.

    File must be located in the `aospy_path` directory and be the same name
    as the desired aospy module subpackage, namely one of `regions`, `calcs`,
    `variables`, and `projects`.
    """
    return imp.load_source(
        name, '/'.join([user_path, name, '__init__.py']).replace('//', '/')
    )


def robust_bool(obj):
    try:
        return bool(obj)
    except ValueError:
        return obj.any()


def get_parent_attr(obj, attr, strict=False):
    """
    Check if the object has the given attribute and it is non-empty.  If not,
    check each parent object for the attribute and use the first one found.
    """
    attr_val = getattr(obj, attr, False)
    if robust_bool(attr_val):
        return attr_val

    else:
        for parent in ('parent', 'var', 'run', 'model', 'proj'):
            parent_obj = getattr(obj, parent, False)
            if parent_obj:
                return get_parent_attr(parent_obj, attr, strict=strict)

        if strict:
            raise AttributeError('Attribute %s not found in parent of %s'
                                 % (attr, obj))
        else:
            return None


def dict_name_keys(objs):
    """Create dict whose keys are the 'name' attr of the objects."""
    assert isinstance(objs, (tuple, list, dict, set))
    if isinstance(objs, (tuple, list, set)):
        try:
            return {obj.name: obj for obj in objs}
        except AttributeError as e:
            raise AttributeError(e)
    return objs


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
        msg = "Input {0} must have length 1 or {1}: len({0})= {2}"
        raise ValueError(msg.format(x, n, len(x)))
    if n == 1:
        if single_to_list:
            return [x]
        return x
    return [x]*n


def var_label(var, level):
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


def data_in_label(intvl_in, dtype_in_time, dtype_in_vert=False):
    """Create string label specifying the input data of a calculation."""
    intvl_lbl = intvl_in
    time_lbl = dtype_in_time
    lbl = '_'.join(['from', intvl_lbl, time_lbl]).replace('__', '_')
    vert_lbl = dtype_in_vert if dtype_in_vert else False
    if vert_lbl:
        lbl = '_'.join([lbl, vert_lbl]).replace('__', '_')
    return lbl


def data_out_label(time_intvl, dtype_time, dtype_vert=False):
    intvl_lbl = time_label(time_intvl, return_val=False)
    time_lbl = dtype_time
    lbl = '.'.join([intvl_lbl, time_lbl]).replace('..', '.')
    vert_lbl = dtype_vert if dtype_vert else False
    if vert_lbl:
        lbl = '.'.join([lbl, vert_lbl]).replace('..', '.')
    return lbl


def ens_label(ens_mem):
    """Create label of the ensemble member for aospy data I/O."""
    if ens_mem in (None, False):
        return ''
    elif ens_mem == 'avg':
        return 'ens_mean'
    else:
        return 'mem' + str(ens_mem + 1)


def yr_label(yr_range):
    """Create label of start and end years for aospy data I/O."""
    assert yr_range is not None, "yr_range is None"
    if yr_range[0] == yr_range[1]:
        return '{:04d}'.format(yr_range[0])
    else:
        return '{:04d}-{:04d}'.format(*yr_range)


def znl_label(var):
    """Create label denoting zonal mean values for aospy data I/O."""
    try:
        if var.def_lon:
            return 'znl'
        else:
            return ''
    except:
        return ''


def time_label(intvl, return_val=True):
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
        labels = {'jfm': (1, 2, 3),
                  'fma': (2, 3, 4),
                  'mam': (3, 4, 5),
                  'amj': (4, 5, 6),
                  'mjj': (5, 6, 7),
                  'jja': (6,  7,  8),
                  'jas': (7, 8, 9),
                  'aso': (8, 9, 10),
                  'son': (9, 10, 11),
                  'ond': (10, 11, 12),
                  'ndj': (11, 12, 1),
                  'djf': (1, 2, 12),
                  'jjas': (6, 7, 8, 9),
                  'djfm': (12, 1, 2, 3),
                  'ann': range(1, 13)}
        for lbl, vals in labels.items():
            if intvl == lbl or set(intvl) == set(vals):
                label = lbl
                value = np.array(vals)
                break
    if return_val:
        return label, value
    else:
        return label


def data_name_gfdl(name, domain, data_type, intvl_type, data_yr,
                      intvl, data_in_start_yr, data_in_dur):
    """Determine the filename of GFDL model data output."""
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
            label, val = time_label(intvl)
        if data_in_dur == 1:
            filename = (domain + '.{:04d}'.format(data_in_yr) +
                        '.' + label + '.nc')
        else:
            filename = (domain + '.{:04d}'.format(data_in_yr) + '-' +
                        '{:04d}'.format(int(data_in_yr+data_in_dur-1)) +
                        '.' + label + '.nc')
    elif data_type == 'av_ts':
        filename = (domain + '.{:04d}'.format(data_in_yr) + '-' +
                    '{:04d}'.format(int(data_in_yr+data_in_dur-1)) +
                    '.01-12.nc')
    return filename


def dmget(files_list):
    """Call GFDL command 'dmget' to access archived files."""
    try:
        subprocess.call(['dmget'] + files_list)
    except OSError:
        logging.warning('dmget command not found in this machine')


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


def get_data_direc_repo(data_direc, var_name, version=-1):
    """Determine the directory containing the needed netCDF files."""
    # Catch if version was specified as False or None.
    if version is None:
        version = -1
    dir_prefix = os.path.realpath(data_direc)
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
