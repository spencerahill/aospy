import imp
import time

class Proj(object):
    """
    Project parameters: models, regions, directories, etc.

    """
    def __init__(self, name, vars={}, direc_out='',
                 nc_dir_struc='gfdl', verbose=True):
        self.verbose = verbose
        if self.verbose:
            print ("Initializing Project instance: %s (%s)"
                   % (name, time.ctime()))
        self.name = name
        self.direc_out = direc_out

        # vars becomes dict of form {var.name: var}
        assert type(vars) in (dict, list, tuple)
        if type(vars) in (list, tuple):
            self.vars = {v.name: v for v in vars}
        else:
            self.vars = vars
        # Set each Var's parent attribute to this project.
        for var in self.vars.values():
            setattr(var, 'proj', self)

    def __str__(self):
        return 'Project instance "' + self.name + '"'
