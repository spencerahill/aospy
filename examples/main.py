#! /usr/bin/env python
"""Template main script for executing calculations using aospy.

Before using this script
------------------------

It is best to copy this template into a separate directory before
populating it with the names from your own object library.  You can also always
get a fresh copy from https://github.com/spencerahill/aospy/

How to use this script
----------------------

As a user, there are only two places you need to modify the code:

1. (Only done once) Replace 'INSERT_YOUR_OBJ_LIB' with the name of your object
   library.  Consult the documentation for instructions on how to make your
   object library findable by Python if this statement generates errors.

   This is the only edit that should be made above the `if __name__ ==
   '__main__'` statement.  The other edits (described next) are all below that
   statement (at the end of the file).

2. (Done as often as desired) Replace the dummy project, model, run, var, and
   region names with the names of your objects that you want to perform
   calculations with.

   Also alter the other parameters -- date_range, etc. -- to your liking.

Running the script
------------------

Once the parameters are all set as desired, execute the script from the command
line ::

  ./main.py  # after `cd` to the directory where you've made your copy

"""
from __future__ import print_function
import datetime
import itertools
import logging

try:
    import colorama
except ImportError:
    pass
try:
    import multiprocess
except ImportError:
    pass

import aospy

# Modify this per instructions above.
from INSERT_YOUR_OBJ_LIB import projs, variables
# No other code should be modified between here and the 'BEGIN POTENTIAL
# MODIFICATIONS' comment below.


class ObjectsForCalc(tuple):
    """Container of aospy objects to be used for a single Calc."""
    def __new__(cls, *objects):
        return super(ObjectsForCalc, cls).__new__(cls, objects)

    def __init__(self, *objects):
        self.objects = objects

    def __repr__(self):
        return 'ObjectsForCalc' + repr(tuple(self.objects))


class MainParams(object):
    """Container for parameters specified in main routine."""
    pass


class MainParamsParser(object):
    """Interface between specified parameters and resulting CalcSuite."""
    def str_to_aospy_obj(self, proj, model, var, region):
        proj_out = aospy.to_proj(proj, self.projs)
        model_out = aospy.to_model(model, proj_out, self.projs)
        var_out = aospy.to_var(var, variables)
        region_out = aospy.to_region(region, self.projs, proj=proj_out)
        return proj_out, model_out, var_out, region_out

    def aospy_obj_to_iterable(self, proj, model, var, region):
        return [aospy.to_iterable(obj) for obj in (proj, model, var, region)]

    def str_to_aospy_iterable(self, proj, model, var, region):
        p, m, v, r = self.str_to_aospy_obj(proj, model, var, region)
        return self.aospy_obj_to_iterable(p, m, v, r)

    def create_child_run_obj(self, models, runs, proj):
        """Create child Run object(s) for each Model object."""
        run_objs = []
        for run in runs:
            for model in models:
                try:
                    run_obj = aospy.to_run(run, model, proj, self.projs)
                    if isinstance(run, ObjectsForCalc):
                        run_objs.append(ObjectsForCalc(run_obj))
                    else:
                        run_objs.append(run_obj)
                except AttributeError as ae:
                    logging.info(str(ae))
            # Retain the original type (e.g. list v. ObjectsForCalc).
            run_objs = type(runs)(run_objs)
        if 'cmip5' in [p.name for p in proj]:
            if isinstance(run_objs[0], list):
                return run_objs[0]
            return [run_objs[0]]
        # If flat list, return the list.  If nested, then flatten it.
        if all([isinstance(r, aospy.Run) for r in run_objs]):
            return run_objs
        return list(itertools.chain.from_iterable(run_objs))

    def __init__(self, main_params, projs):
        """Turn all inputs into aospy-ready objects."""
        self.__dict__ = vars(main_params)
        self.projs = projs
        self.proj, self.model, self.var, self.region = (
            self.str_to_aospy_iterable(main_params.proj, main_params.model,
                                       main_params.var, main_params.region)
            )
        self.run = self.create_child_run_obj(self.model, self.run, self.proj)
        self.region = [aospy.utils.io.dict_name_keys(self.region)]


class CalcSuite(object):
    """Creates suite of Calc objects based on inputted specifications. """
    def __init__(self, calc_suite_interface):
        self.__dict__ = vars(calc_suite_interface)

    def print_params(self):
        pairs = (
            ('Project', self.proj),
            ('Models', self.model),
            ('Runs', self.run),
            ('Ensemble members', self.ens_mem),
            ('Variables', self.var),
            ('Year ranges', self.date_range),
            ('Geographical regions', [r.values() for r in self.region]),
            ('Time interval of input data', self.intvl_in),
            ('Time interval for averaging', self.intvl_out),
            ('Input data time type', self.dtype_in_time),
            ('Input data vertical type', self.dtype_in_vert),
            ('Output data time type', self.dtype_out_time),
            ('Output data vertical type', self.dtype_out_vert),
            ('Vertical levels', self.level),
            ('Time offset', self.time_offset),
            ('Compute this data', self.compute),
        )
        print('')
        try:
            colorama.init()
            color_left = colorama.Fore.BLUE
            color_right = colorama.Fore.RESET
            reset = colorama.Style.RESET_ALL
        except NameError:
            color_left = ''
            color_right = ''
            reset = ''
        for left, right in pairs:
            print(color_left, left, ':', color_right, right)
        print(reset)

    def prompt_user_verify(self):
        try:
            input = raw_input
        except NameError:
            import builtins
            input = builtins.input
        if not input("Perform these computations? ").lower() in ('y', 'yes'):
            raise IOError('\n', 'Execution cancelled by user.')

    def create_params_all_calcs(self):
        attr_names = ('proj',
                      'model',
                      'run',
                      'ens_mem',
                      'var',
                      'date_range',
                      'level',
                      'region',
                      'intvl_in',
                      'intvl_out',
                      'dtype_in_time',
                      'dtype_out_time',
                      'dtype_in_vert',
                      'dtype_out_vert',
                      'time_offset')
        attrs = tuple([getattr(self, name) for name in attr_names])
        # Each permutation becomes a dictionary, with the keys being the attr
        # names and the values being the corresponding value for that
        # permutation.  These dicts can then be directly passed to the
        # CalcInterface class to make the Calc objects.
        permuter = itertools.product(*attrs)
        param_combos = []
        for permutation in permuter:
            param_combos.append(dict(zip(attr_names, permutation)))
        return param_combos

    def create_calcs(self, param_combos, exec_calcs=False):
        """Iterate through given parameter combos, creating needed Calcs."""
        calcs = []
        for params in param_combos:
            try:
                ci = aospy.CalcInterface(**params)
            except:
                raise
            calc = aospy.Calc(ci)
            if exec_calcs:
                try:
                    calc.compute()
                except RuntimeError as e:
                    logging.warn(repr(e))
                except IOError as e:
                    logging.warn(repr(e))
                except:
                    raise
            calcs.append(calc)
        return calcs

    def exec_calcs(self, calcs):
        out = []
        for calc in calcs:
            try:
                o = calc.compute()
            except RuntimeError as e:
                logging.warn(repr(e))
            else:
                out.append(o)
        return out

    def print_results(self, calcs):
        for calc in calcs:
            for region in calc.region.values():
                print([calc.load(self.dtype_out_time[0],
                                 # dtype_out_vert=params[-2],
                                 region=region, plot_units=True)])


def main(main_params, exec_calcs=True, prompt_verify=True):
    """Main script for interfacing with aospy."""
    # Instantiate objects and load default/all models, runs, and regions.
    cs = CalcSuite(MainParamsParser(main_params, projs))
    cs.print_params()
    if prompt_verify:
        try:
            cs.prompt_user_verify()
        except IOError as e:
            logging.warn(repr(e))
            return
    param_combos = cs.create_params_all_calcs()
    if main_params.parallelize and exec_calcs:
        calcs = cs.create_calcs(param_combos, exec_calcs=False)
        p = multiprocess.Pool()
        return p.map(lambda calc: calc.compute(), calcs)
    else:
        calcs = cs.create_calcs(param_combos, exec_calcs=exec_calcs)
    return calcs


if __name__ == '__main__':

    mp = MainParams()

    # BEGIN POTENTIAL MODIFICATIONS
    # Modifications are made below this point, per instructions above.
    # Consult `CalcInterface` API reference for further explanation
    # of each option and accepted values.

    # Name of the project.
    mp.proj = 'myproj'
    # List of the model names.
    mp.model = ['mymodel1']
    # List of the run names, or 'default', or 'all'.
    mp.run = ['myrun1', 'myrun2']
    # List of the variable names.
    mp.var = ['precip', 'pot_temp']
    # List of the region names, or 'all'.
    mp.region = 'all'

    # Start and end dates (inclusive).  List of tuples and/or 'default'.
    # If list of tuples, tuples are of form (start_date, end_date), where
    # start_date and end_date are datetime.datetime objects.
    mp.date_range = ['default']

    # Sub-annual time-interval to average over.  List of 'ann', seasonal string
    # (e.g. 'djf'), or month integer (1 for Jan, 2 for Feb, etc).
    mp.intvl_out = ['ann']
    # List of tuples, each listing the desired spatiotemporal reductions.
    mp.dtype_out_time = [('av', 'reg.av')]
    # List of desired vertical reductions to perform.
    mp.dtype_out_vert = [None]

    # List of time spacing of input data.
    mp.intvl_in = ['monthly']
    # List of time type of input data.
    mp.dtype_in_time = ['ts']
    # List of vertical data type of input data.
    mp.dtype_in_vert = [False]

    # List the time offset dictionaries (if desired) to apply to the input
    # data (e.g. [{'days': -15}]).
    mp.time_offset = [None]

    # Submit all calculations in parallel.  Requires 'multiprocess' package
    # (which can be obtained e.g. via `pip install multiprocess`).
    mp.parallelize = True

    # Don't modify this statement.
    calcs = main(mp)
