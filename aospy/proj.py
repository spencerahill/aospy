"""proj.py: aospy.Proj class for organizing work in single project."""
import time

from . import utils


class Proj(object):
    """Project parameters: models, regions, directories, etc."""
    def __init__(self, name, vars={}, models={}, default_models={}, regions={},
                 direc_out='', tar_direc_out='', nc_dir_struc=False,
                 verbose=True):
        self.verbose = verbose
        if self.verbose:
            print ("Initializing Project instance: %s (%s)"
                   % (name, time.ctime()))
        self.name = name
        self.direc_out = direc_out
        self.tar_direc_out = tar_direc_out
        self.nc_dir_struc = nc_dir_struc

        if models:
            self.models = utils.io.dict_name_keys(models)
        else:
            self.models = {}
        if default_models == 'all':
            self.default_models = self.models
        elif default_models:
            self.default_models = utils.io.dict_name_keys(default_models)
        else:
            self.default_models = {}
        if regions:
            self.regions = utils.io.dict_name_keys(regions)
        else:
            self.regions = {}

        for obj_dict in (self.models, self.regions):
            for obj in obj_dict.values():
                setattr(obj, 'proj', self)

    def __str__(self):
        return 'Project instance "' + self.name + '"'

    __repr__ = __str__
