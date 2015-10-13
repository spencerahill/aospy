"""`Run` class; for storing attributes of a model run or obs product."""
import numpy as np

from .timedate import TimeManager


class Run(object):
    """Model run parameters."""
    def _set_direc(self, direc_nc, ens_mem_prefix, ens_mem_ext,
                   ens_mem_suffix):
        """Set the list of paths containing the Run's netCDF data."""
        if all((ens_mem_prefix, ens_mem_ext, ens_mem_suffix)):
            return [ens_mem_prefix + ext + ens_mem_suffix
                    for ext in ens_mem_ext]
        return direc_nc

    def __init__(self, name='', description='', proj=False, direc_nc=False,
                 nc_dur=False, nc_start_date=False, nc_end_date=False,
                 nc_dir_struc='gfdl', nc_suffix=False, nc_files={},
                 default_date_range=False, ens_mem_prefix=False,
                 ens_mem_ext=False, ens_mem_suffix=False, tags=()):
        """Instantiate a `Run` object."""
        self.name = name
        self.description = description
        self.proj = proj

        self.nc_dur = nc_dur
        self.nc_dir_struc = nc_dir_struc
        self.nc_suffix = nc_suffix
        self.nc_files = nc_files
        self.nc_start_date = TimeManager.to_datetime(nc_start_date)
        self.nc_end_date = TimeManager.to_datetime(nc_end_date)
        try:
            self.default_date_range = tuple([TimeManager.to_datetime(d)
                                             for d in default_date_range])
        except:
            self.default_date_range = (self.nc_start_date, self.nc_end_date)

        self.tags = tags

        self.ens_mem_prefix = ens_mem_prefix
        self.ens_mem_ext = ens_mem_ext
        self.ens_mem_suffix = ens_mem_suffix
        self.direc_nc = self._set_direc(direc_nc, ens_mem_prefix,
                                        ens_mem_ext, ens_mem_suffix)

    def __str__(self):
        return 'Run instance "%s"' % self.name

    __repr__ = __str__
