from .proj import Proj, proj_inst
from .model import Model, model_inst

class Run(object):
    """Model run parameters."""
    def __init__(self, name='', description='', proj=False, direc_nc=False,
                 nc_dur=False, nc_start_yr=False, nc_end_yr=False,
                 nc_dir_struc='gfdl', default_yr_range=False,
                 ens_mem_prefix=False, ens_mem_ext=False, ens_mem_suffix=False,
                 tags=()):
        self.name = name
        self.description = description
        self.proj = proj
        self.nc_dur = nc_dur
        self.nc_start_yr = nc_start_yr
        self.nc_end_yr = nc_end_yr
        self.nc_dir_struc = nc_dir_struc
        self.default_yr_range = default_yr_range
        self.tags = tags

        self.ens_mem_prefix = ens_mem_prefix
        self.ens_mem_ext = ens_mem_ext
        self.ens_mem_suffix = ens_mem_suffix
        self._set_direc(direc_nc, ens_mem_prefix, ens_mem_ext, ens_mem_suffix)

    def __str__(self):
        return 'Run instance "%s"' % self.name

    def _set_direc(self, direc_nc, ens_mem_prefix, ens_mem_ext, ens_mem_suffix):
        """Set the list of paths containing the Run's netCDF data."""
        if all((ens_mem_prefix, ens_mem_ext, ens_mem_suffix)):
            self.direc_nc = [ens_mem_prefix + ext + ens_mem_suffix
                             for ext in ens_mem_ext]
        else:
            self.direc_nc = direc_nc

def run_inst(run, parent_model=False, parent_proj=False):
    """Convert string matching an aospy.run name to an aospy.run instance."""
    orig_type = type(run)
    if parent_proj and type(parent_proj) is not Proj:
        parent_proj = proj_inst(parent_proj)
    if parent_model and type(parent_model) is not Model:
        parent_model = model_inst(parent_model, parent_proj)
    if type(run) is Run:
        pass
    elif type(run) is str:
        run = parent_model.runs[run]
    elif type(run) in (list, tuple):
        run = [run_inst(rn, mod, parent_proj) for (rn, mod)
               in zip(run, parent_model)]
        if orig_type is tuple:
            run = tuple(run)
    else:
        raise TypeError
    if parent_model:
        try:
            run.model = parent_model
        except AttributeError:
            pass
    return run

    
