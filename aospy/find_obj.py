"""Get aospy objects from an object containing them given their name string."""
from collections import Sequence

from . import Proj, Model, Run, Var, Region, Operator
from . import utils


def to_iterable(obj):
    """Return the object if already iterable, otherwise return it as a list."""
    try:
        zip([], obj)
    except TypeError:
        return [obj]
    else:
        return obj


def to_proj(proj, projs_module):
    """Convert string of an aospy.Proj name to an aospy.Proj instance."""
    orig_type = type(proj)
    if isinstance(proj, Proj):
        return proj

    elif isinstance(proj, str):
        try:
            return getattr(projs_module, proj)
        except AttributeError:
            raise AttributeError('Not a recognized Proj name: %s' % proj)

    elif isinstance(proj, Sequence):
        proj = [to_proj(pr, projs_module) for pr in proj]
        if orig_type is list:
            return proj
        else:
            return tuple(proj)

    else:
        raise TypeError


def to_model(model, proj, projs_module):
    """Convert string of a Model name to a Model instance."""
    orig_type = type(model)
    proj = to_proj(proj, projs_module)

    if isinstance(model, Model):
        return model

    elif isinstance(model, str):
        if model == 'all':
            model = proj.models.values()
        elif model == 'default':
            model = proj.default_models.values()
        else:
            model = proj.models[model]
        return model

    elif isinstance(model, Sequence):
        model = [to_model(mod, pr, projs_module) for (mod, pr)
                 in zip(model, utils.io.to_dup_list(proj, len(model)))]
        if orig_type is tuple:
            model = tuple(model)
        return model

    else:
        raise TypeError


def to_run(run, model, proj, projs_module):
    """Convert string matching an aospy.run name to an aospy.run instance."""
    orig_type = type(run)
    proj = to_proj(proj, projs_module)
    model = to_model(model, proj, projs_module)

    if isinstance(run, Run):
        return run

    elif isinstance(run, str):
        if run == 'default':
            return model.default_runs.values()
        elif run == 'all':
            return model.runs
        else:
            try:
                return model.runs[run]
            except KeyError:
                raise AttributeError("Model '{}' has no run '{}'.".format(
                    model, run))

    elif isinstance(run, Operator):
        operator = run.operator
        runs = to_run(run.objects, model, proj, projs_module)
        return Operator(operator, runs)

    elif isinstance(run, dict):
        # Assume run(s) is/are key(s), not value(s).
        vals = run.values()
        keys = to_run(run.keys(), model, proj, projs_module)
        return dict(zip(keys, vals))

    elif isinstance(run, Sequence):
        run = [to_run(rn, mod, pr, projs_module) for (rn, mod, pr)
               in zip(run, utils.io.to_dup_list(model, len(run)),
                      utils.io.to_dup_list(proj, len(run)))]
        if orig_type is tuple:
            run = tuple(run)
        return run

    else:
        raise TypeError


def to_var(var, vars_module):
    """Convert string of an aospy.var name to an aospy.var instance."""
    if isinstance(var, Var):
        var_out = var

    elif isinstance(var, str):
        try:
            var_out = getattr(vars_module, var)
        except AttributeError:
            raise AttributeError('Not a recognized Var name: %s' % var)

    elif isinstance(var, Sequence):
        var_out = [to_var(v, vars_module) for v in var]
        if isinstance(var, tuple):
            var_out = tuple(var_out)

    else:
        raise TypeError

    return var_out


def to_region(region, regions_module, proj=False):
    """Convert string of an aospy.Region name to an aospy.Region instance."""
    if isinstance(region, Region):
        return region

    elif isinstance(region, bool) or region is None:
        return region

    elif isinstance(region, str):
        if proj and region == 'all':
            return proj.regions
        else:
            return getattr(regions_module, region)

    elif isinstance(region, Sequence):
        region_out = [to_region(r, regions_module, proj=proj) for r in region]
        if isinstance(region, tuple):
            region_out = tuple(region_out)
        return region_out

    else:
        raise TypeError
