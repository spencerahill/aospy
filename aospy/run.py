"""Functionality for storing attributes of a model run or obs product."""
from .utils.times import datetime_or_default


def _set_direc(data_direc, ens_mem_prefix, ens_mem_ext,
               ens_mem_suffix):
    """Set the list of paths containing the Run's netCDF data."""
    if all((ens_mem_prefix, ens_mem_ext, ens_mem_suffix)):
        return [ens_mem_prefix + ext + ens_mem_suffix
                for ext in ens_mem_ext]
    return data_direc


class Run(object):
    """Model run parameters."""
    def __init__(self, name=None, description=None, proj=None,
                 data_direc=None, data_dur=None,
                 data_start_date=None, data_end_date=None,
                 data_dir_struc='gfdl', data_suffix=None,
                 data_files=None, default_start_date=None,
                 default_end_date=None, ens_mem_prefix=None,
                 ens_mem_ext=None, ens_mem_suffix=None, tags=None,
                 idealized=False):
        """Instantiate a `Run` object."""
        self.name = '' if name is None else name
        self.description = '' if description is None else description
        self.proj = proj

        self.data_dur = data_dur
        self.data_dir_struc = data_dir_struc
        self.data_suffix = data_suffix
        self.data_files = {} if data_files is None else data_files

        self.data_start_date = datetime_or_default(data_start_date, None)
        self.data_end_date = datetime_or_default(data_end_date, None)
        self.default_start_date = datetime_or_default(default_start_date,
                                                      self.data_start_date)
        self.default_end_date = datetime_or_default(default_end_date,
                                                    self.data_end_date)

        self.tags = [] if tags is None else tags
        self.idealized = idealized
        self.ens_mem_prefix = ens_mem_prefix
        self.ens_mem_ext = ens_mem_ext
        self.ens_mem_suffix = ens_mem_suffix
        self.data_direc = _set_direc(data_direc, ens_mem_prefix, ens_mem_ext,
                                     ens_mem_suffix)

    def __str__(self):
        return 'Run instance "%s"' % self.name

    __repr__ = __str__
