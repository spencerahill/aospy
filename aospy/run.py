"""run.py: Run class of aospy for storing attributes of a GCM run."""

class Run(object):
    """Model run parameters."""
    def __init__(self, name='', description='', proj=False, direc_nc=False,
                 nc_dur=False, nc_start_yr=False, nc_end_yr=False,
                 nc_start_month=False, nc_end_month=False, nc_dir_struc='gfdl',
                 nc_suffix=False, nc_files={}, default_yr_range=False,
                 ens_mem_prefix=False, ens_mem_ext=False, ens_mem_suffix=False,
                 tags=(), read_mode='netcdf4', nc_start_day=None,
                 nc_end_day=None, default_time_range=None):
        self.name = name
        self.description = description
        self.proj = proj
        self.nc_dur = nc_dur
        self.nc_start_yr = nc_start_yr
        self.nc_end_yr = nc_end_yr
        self.nc_start_month = nc_start_month
        self.nc_end_month = nc_end_month
        self.nc_dir_struc = nc_dir_struc
        self.nc_suffix = nc_suffix
        self.nc_files = nc_files
        self.default_yr_range = default_yr_range
        self.tags = tags
        
        self.read_mode = read_mode
        self.nc_start_day = nc_start_day
        self.nc_end_day = nc_end_day
        self.default_time_range = default_time_range 

        self.ens_mem_prefix = ens_mem_prefix
        self.ens_mem_ext = ens_mem_ext
        self.ens_mem_suffix = ens_mem_suffix
        self._set_direc(direc_nc, ens_mem_prefix, ens_mem_ext, ens_mem_suffix)

    def __str__(self):
        return 'Run instance "%s"' % self.name

    __repr__ = __str__

    def _set_direc(self, direc_nc, ens_mem_prefix, ens_mem_ext, ens_mem_suffix):
        """Set the list of paths containing the Run's netCDF data."""
        if all((ens_mem_prefix, ens_mem_ext, ens_mem_suffix)):
            self.direc_nc = [ens_mem_prefix + ext + ens_mem_suffix
                             for ext in ens_mem_ext]
        else:
            self.direc_nc = direc_nc
