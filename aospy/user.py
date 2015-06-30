"""user.py: Load user's own library of aospy objects into aospy."""
import imp

from . import user_path, Proj, Model, Run, Var, Region

def proj_inst(proj):
    """Convert string of an aospy.Proj name to an aospy.Proj instance."""
    orig_type = type(proj)
    if type(proj) is Proj:
        return proj
    elif type(proj) is str:
        try:
            proj_module = imp.load_source(proj, (user_path + '/' +
                                          proj + '.py').replace('//','/'))
            proj_func = getattr(proj_module, proj)
            return proj_func()
        except AttributeError:
            raise AttributeError('Not a recognized aospy.Proj name: %s'
                                 % proj)
    elif type(proj) in (list, tuple):
        proj = [proj_inst(pr) for pr in proj]
        if orig_type is list:
            return proj
        else:
            return tuple(proj)
    else:
        raise TypeError

def model_inst(model, parent_proj=False):
    """Convert string of a Model name to a Model instance."""
    orig_type = type(model)
    if type(parent_proj) is str:
        parent_proj_out = proj_inst(parent_proj)
    if type(model) is Model:
        pass
    elif type(model) is str:
        model_out = parent_proj.models[model]
    elif type(model) in (list, tuple):
        model_out = [model_inst(mod, proj) for (mod, proj)
                     in zip(model, parent_proj_out)]
        # Retain original container type, i.e. list or tuple.
        if orig_type is tuple:
            model_out = tuple(model_out)
    if parent_proj:
        parent_proj_out = proj_inst(parent_proj)
        try:
            model_out.proj = parent_proj
        except AttributeError:
            pass
    return model_out


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


def var_inst(var):
    """Convert string of an aospy.var name to an aospy.var instance."""
    if type(var) is Var:
        var_out = var
    elif type(var) is str:
        try:
            var_out = getattr(variables, var)
        except AttributeError:
            raise AttributeError('Not a recognized aospy.Var name: %s'
                                 % var)
    elif type(var) in (list, tuple):
        var_out = [var_inst(v) for v in var]
        if type(var) is tuple:
            var_out = tuple(var_out)
    return var_out

def region_inst(region):
    """Convert string of an aospy.Region name to an aospy.Region instance."""
    if type(region) is str:
        return getattr(regions, region)
    else:
        return region

regions = imp.load_source(
    'regions', (user_path + '/regions/__init__.py').replace('//','/')
)
units = imp.load_source(
    'units', (user_path + '/units/__init__.py').replace('//','/')
)
variables = imp.load_source(
    'variables', (user_path + '/variables/__init__.py').replace('//','/')
)
runs = imp.load_source(
    'runs', (user_path + '/runs/__init__.py').replace('//','/')
)
models = imp.load_source(
    'models', (user_path + '/models/__init__.py').replace('//','/')
)
