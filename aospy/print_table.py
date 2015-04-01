import cPickle
from aospy import projects
from aospy.io import load_plot_data

def print_data(proj, model, exp, ens_mem, var, level, intvl, dtype, yr_range,
               region):
    try:
        d = load_plot_data(proj, model, exp, ens_mem, var, level, intvl, dtype,
                           yr_range, region=region)
        if var in ('precip', 'prec_conv', 'prec_ls', 'p-e', 'evap'):
            # CRU output is in mm, not mm/day, so divide by the avg
            # month length (365.25/12=30.4375).  This is imperfect:
            # need to divide by each month's length.
            if model == 'cru':
                d/=30.4375
            # Similarly, U. Delaware data is in cm.
            elif model == 'udel':
                d/=3.4375
        try:
            print '%s\t%.2f' % (exp, round(d, 2))
        except TypeError:
            print exp, [round(val, 2) for val in d]
        del d
    except (KeyError, IOError):
        pass

def print_table(proj, models, exps, vars_, region, intvl,
                dtype='reg.av', yr_range='default', ens_mem=None, level=None):
    proj = getattr(projects, proj)()
    if models in ('default', ['default']):
        models = [getattr(model, 'name') for model in
                  proj.default_models.itervalues()]
    elif models in ('all', ['all']):
        models = [getattr(model, 'name') for model in
                  proj.models.itervalues()]
    for model in models:
        model = proj.models[model]
        if exps in ('default', ['default']):
            exps_conv = [getattr(run, 'name') for run in
                         model.default_runs.itervalues()]
        elif exps in ('all', ['all']):
            exps_conv = [getattr(run, 'name') for run in
                    model.runs.itervalues()]
        else:
            exps_conv = exps
        for var in vars_:
            print '\n%s\t%s\t%s' % (model, var, region)
            for exp in exps_conv:
                print_data(proj, model, exp, ens_mem, var, level,
                           intvl, dtype, yr_range, region)
