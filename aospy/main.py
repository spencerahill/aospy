"""Main script for 'aospy' module."""
import itertools
import aospy
from aospy.io import mult_var_calcs, _proj_inst, _get_time_avg_var, _model_inst

do_time_av = False
proj_av = 'aero_3agcm'
models_av = ['default']
runs_av = ['default']
ens_mem_av = [None]
vars_av = ['slp']
levs_av = [None]
regions_av = ['all']
intvls_av = ['jas']
intvl_type_av = ['monthly']
yr_ranges_av = ['default']
# yr_ranges_av = [(1979,2011), (1979, 2013)]
data_type_av = ['ts']
out_av = ['av', 'std', 'reg.av', 'reg.ts', 'reg.std']

do_calcs = True
proj_calc = 'aero_3agcm'
models_calc = ['am2', 'am3']
runs_calc = ['default']
ens_mem_calc = [None]
calcs = {'gross_moist_stab':{}, 'gross_dry_stab':{}, 'moisture_strat':{}}
levs_calc = [False]
regions_calc = ['all']
intvls_calc = ['jas']
intvl_type_calc = ['monthly']
yr_ranges_calc = ['default']
data_type_calc = ['ts']
out_calc = ['av', 'std', 'reg.av', 'reg.ts', 'reg.std']

# Print diagnostics.
if do_time_av:
    print '\n***Time averaging***'
    print '\nProject:', proj_av
    print 'Models:', models_av
    print 'Runs:', runs_av
    print 'Ensemble members:', ens_mem_av
    print 'Variables:', vars_av
    print 'Pressure levels:', levs_av
    print 'Geographical regions:', regions_av
    print 'Time intervals for averaging:', intvls_av
    print 'Time interval of input data:', intvl_type_av
    print 'Year ranges:', yr_ranges_av
    print 'Input data type:', data_type_av
    print 'Output data types:', out_av
if do_calcs:
    print '\n***Calculations***'
    print '\nProject:', proj_calc
    print 'Models:', models_calc
    print 'Runs:', runs_calc
    print 'Ensemble members:', ens_mem_calc
    print 'Calculations:', calcs
    print 'Pressure levels:', levs_calc
    print 'Geographical regions:', regions_calc
    print 'Time intervals for averaging:', intvls_calc
    print 'Time interval of input data:', intvl_type_calc
    print 'Year ranges:', yr_ranges_calc
    print 'Input data type:', data_type_calc
    print 'Output data types:', out_calc

if raw_input("\nProceed using these parameters? ").lower() == 'y':
    if do_time_av:
        proj_av = _proj_inst(proj_av)
        if models_av in ('all', ['all']):
            models_av = proj_av.models.values()
        elif models_av in ('default', ['default']):
            models_av = proj_av.default_models
        if regions_av in ('all', ['all']):
            regions_av = proj_av.regions
        kwargs_av = {'calcs_out': out_av, 'regions': regions_av}
        print '\n\tVariable time averages and statistics:'
        for model in models_av:
            if runs_av in ('default', ['default']):
                mod = _model_inst(model, proj_av)
                runs_a = mod.default_runs.keys()
            else:
                runs_a = runs_av

            for params in itertools.product(
                    [proj_av], [model], runs_a, vars_av, ens_mem_av, levs_av,
                    intvls_av, intvl_type_av, data_type_av, yr_ranges_av
            ):
                Var_av = _get_time_avg_var(proj_av, params[1], params[2],
                                          params[3])
                Var_av.var_time_avg(*params[4:], **kwargs_av)

    if do_calcs:
        proj_calc = _proj_inst(proj_calc)
        if models_calc in ('all', ['all']):
            models_calc = proj_calc.models.values()
        elif models_calc in ('default', ['default']):
            models_calc = proj_calc.default_models
        if regions_calc in ('all', ['all']):
            regions_calc = proj_calc.regions
        kwargs_calc = {'calcs_out': out_calc, 'regions': regions_calc}
        print '\n\tCalculations:'
        for model in models_calc:
            if runs_calc in ('default', ['default']):
                mod = _model_inst(model, proj_calc)
                runs_c = mod.default_runs.keys()
            else:
                runs_c = runs_calc
            for params in itertools.product(
                    [proj_calc], [model], runs_c, ens_mem_calc,
                    calcs.iteritems(), levs_calc, intvls_calc, intvl_type_calc,
                    data_type_calc, yr_ranges_calc
            ):
                mult_var_calcs(*params, **kwargs_calc)
else:
    print '\nExecution cancelled by user. '
