"""proj.py: aospy.Proj class for organizing work in single project."""
import time

from .utils import dict_name_keys


class Proj(object):
    """Project parameters: models, regions, directories, etc."""
    def __init__(self, name, vars={}, models={}, default_models={}, regions={},
                 direc_out='', nc_dir_struc=False, verbose=True):
        self.verbose = verbose
        if self.verbose:
            print ("Initializing Project instance: %s (%s)"
                   % (name, time.ctime()))
        self.name = name
        self.direc_out = direc_out
        self.nc_dir_struc = nc_dir_struc

        self.vars = dict_name_keys(vars)
        self.models = dict_name_keys(models)
        self.default_models = dict_name_keys(default_models)
        self.regions = dict_name_keys(regions)

        for obj_dict in (self.vars, self.models, self.regions):
            for obj in obj_dict.values():
                setattr(obj, 'proj', self)

    def __str__(self):
        return 'Project instance "' + self.name + '"'
