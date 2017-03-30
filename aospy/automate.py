"""Functionality for specifying and cycling through multiple calculations."""
from __future__ import print_function

import itertools
import logging
import pprint
import traceback

from .calc import Calc, CalcInterface
from .region import Region
from .var import Var

try:
    import multiprocess
except ImportError:
    pass


_OBJ_LIB_STR = 'library'
_PROJECTS_STR = 'projects'
_MODELS_STR = 'models'
_RUNS_STR = 'runs'
_REGIONS_STR = 'regions'
_VARIABLES_STR = 'variables'
_TAG_ATTR_MODIFIERS = dict(all='', default='default_')


def _get_attr_by_tag(obj, tag, attr_name):
    """Get attribute from an object via a string tag.

    Parameters
    ----------
    obj : object from which to get the attribute
    attr_name : str
        Unmodified name of the attribute to be found.  The actual attribute
        that is returned may be modified be 'tag'.
    tag : str
        Tag specifying how to modify 'attr_name' by pre-pending it with 'tag'.
        Must be a key of the _TAG_ATTR_MODIFIERS dict.

    Returns
    -------
    the specified attribute of obj
    """
    attr_name = _TAG_ATTR_MODIFIERS[tag] + attr_name
    return getattr(obj, attr_name)


def _permuted_dicts_of_specs(specs):
    """Create {name: value} dict, one each for every permutation.

    Each permutation becomes a dictionary, with the keys being the attr names
    and the values being the corresponding value for that permutation.  These
    dicts can then be directly passed to the CalcInterface class to make the
    Calc objects.

    """
    permuter = itertools.product(*specs.values())
    return [dict(zip(specs.keys(), perm)) for perm in permuter]


def _merge_dicts(*dict_args):
    """Merge the given dictionaries into single dict.

    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    From http://stackoverflow.com/a/26853961/1706640
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def _input_func_py2_py3():
    """Find function for reading user input that works on Python 2 and 3.

    See e.g. http://stackoverflow.com/questions/21731043
    """
    try:
        input = raw_input
    except NameError:
        import builtins
        input = builtins.input
    return input


class AospyException(Exception):
    """Base exception class for the aospy package."""
    pass


def _user_verify(input_func=_input_func_py2_py3(),
                 prompt='Perform these computations? [y/n] '):
    """Prompt the user for verification."""
    if not input_func(prompt).lower()[0] == 'y':
        raise AospyException('Execution cancelled by user.')


def _get_all_objs_of_type(type_, parent):
    """Get all attributes of the given type from the given object.

    Parameters
    ----------
    type_ : The desired type
    parent : The object from which to get the attributes with type matching
        'type_'

    Returns
    -------
    A list (possibly empty) of attributes from 'parent'
    """
    return set([obj for obj in parent.__dict__.values()
                if isinstance(obj, type_)])


class CalcSuite(object):
    """Suite of Calc objects generated from provided specifications."""

    _CORE_SPEC_NAMES = {_OBJ_LIB_STR, _PROJECTS_STR, _MODELS_STR, _RUNS_STR}
    _AUX_SPEC_NAMES = {_VARIABLES_STR,
                       _REGIONS_STR,
                       'date_ranges',
                       'input_time_intervals',
                       'input_time_datatypes',
                       'input_time_offsets',
                       'input_vertical_datatypes',
                       'output_time_intervals',
                       'output_time_regional_reductions',
                       'output_vertical_reductions'}
    _NAMES_SUITE_TO_CALC = {
        _PROJECTS_STR: 'proj',
        _MODELS_STR: 'model',
        _RUNS_STR: 'run',
        _VARIABLES_STR: 'var',
        _REGIONS_STR: 'region',
        'date_ranges': 'date_range',
        'input_time_intervals': 'intvl_in',
        'input_time_datatypes': 'dtype_in_time',
        'input_time_offsets': 'time_offset',
        'input_vertical_datatypes': 'dtype_in_vert',
        'output_time_intervals': 'intvl_out',
        'output_time_regional_reductions': 'dtype_out_time',
        'output_vertical_reductions': 'dtype_out_vert',
    }

    def __init__(self, calc_suite_specs):
        self._specs_in = calc_suite_specs
        self._obj_lib = self._specs_in[_OBJ_LIB_STR]

    def _get_requested_spec(self, obj, spec_name):
        """Helper to translate user specifications to needed objects."""
        requested = self._specs_in[spec_name]
        if isinstance(requested, str):
            return _get_attr_by_tag(obj, requested, spec_name)
        else:
            return requested

    def _permute_core_specs(self):
        """Generate all requested combinations of the core objects."""
        obj_trees = []
        projects = self._get_requested_spec(self._obj_lib, _PROJECTS_STR)
        for project in projects:
            models = self._get_requested_spec(project, _MODELS_STR)
            for model in models:
                runs = self._get_requested_spec(model, _RUNS_STR)
                for run in runs:
                    obj_trees.append({
                        self._NAMES_SUITE_TO_CALC[_PROJECTS_STR]: project,
                        self._NAMES_SUITE_TO_CALC[_MODELS_STR]: model,
                        self._NAMES_SUITE_TO_CALC[_RUNS_STR]: run,
                    })
        return obj_trees

    def _get_regions(self):
        """Get the requested regions."""
        if self._specs_in[_REGIONS_STR] == 'all':
            return [_get_all_objs_of_type(
                Region, getattr(self._obj_lib, 'regions', self._obj_lib)
            )]
        else:
            return [set(self._specs_in[_REGIONS_STR])]

    def _get_variables(self):
        """Get the requested variables."""
        if self._specs_in[_VARIABLES_STR] == 'all':
            return _get_all_objs_of_type(
                Var, getattr(self._obj_lib, 'variables', self._obj_lib)
            )
        else:
            return set(self._specs_in[_VARIABLES_STR])

    def _get_date_ranges(self):
        """Parse the input to get the desired date ranges."""
        if self._specs_in['date_ranges'] == 'default':
            return ['default']
        else:
            return self._specs_in['date_ranges']

    def _get_time_reg_reducts(self):
        """Parse the input to get the desired spatiotemporal reductions."""
        return [self._specs_in['output_time_regional_reductions']]

    def _get_aux_specs(self):
        """Get and pre-process all of the non-core specifications."""
        # Drop the "core" specifications, which are handled separately.
        specs = self._specs_in.copy()
        [specs.pop(core) for core in self._CORE_SPEC_NAMES]

        specs[_REGIONS_STR] = self._get_regions()
        specs[_VARIABLES_STR] = self._get_variables()
        specs['date_ranges'] = self._get_date_ranges()
        specs['output_time_regional_reductions'] = self._get_time_reg_reducts()

        return specs

    def _permute_aux_specs(self):
        """Generate all permutations of the non-core specifications."""
        # Convert to attr names that Calc is expecting.
        calc_aux_mapping = self._NAMES_SUITE_TO_CALC.copy()
        # Special case: manually add 'library' to mapping
        calc_aux_mapping[_OBJ_LIB_STR] = None
        [calc_aux_mapping.pop(core) for core in self._CORE_SPEC_NAMES]

        specs = self._get_aux_specs()
        for suite_name, calc_name in calc_aux_mapping.items():
            specs[calc_name] = specs.pop(suite_name)
        return _permuted_dicts_of_specs(specs)

    def _combine_core_aux_specs(self):
        """Combine permutations over core and auxilliary Calc specs."""
        all_specs = []
        for core_dict in self._permute_core_specs():
            for aux_dict in self._permute_aux_specs():
                all_specs.append(_merge_dicts(core_dict, aux_dict))
        return all_specs

    def create_calcs(self):
        """Generate a Calc object for each requested parameter combination."""
        return [Calc(CalcInterface(**sp)) for sp in
                self._combine_core_aux_specs()]


def _compute_or_skip_on_error(calc, compute_kwargs):
    """Execute the Calc, catching and logging exceptions, but don't re-raise.

    Prevents one failed calculation from stopping a larger requested set
    of calculations.
    """
    try:
        return calc.compute(**compute_kwargs)
    except Exception as e:
        msg = ("Skipping aospy calculation `{0}` due to error with the "
               "following traceback: \n{1}")
        logging.warn(msg.format(calc, traceback.format_exc()))
        return None


def _exec_calcs(calcs, parallelize=False, **compute_kwargs):
    """Execute the given calculations.

    Parameters
    ----------
    calcs : Sequence of ``aospy.Calc`` objects
    parallelize : bool, default False
        Whether to submit the calculations in parallel or not
    compute_kwargs : dict of keyword arguments passed to ``Calc.compute``

    Returns
    -------
    A list of the values returned by each Calc object that was executed.
    """
    if parallelize:
        pool = multiprocess.Pool()
        return pool.map(lambda calc:
                        _compute_or_skip_on_error(calc, compute_kwargs),
                        calcs)
    else:
        return [_compute_or_skip_on_error(calc, compute_kwargs)
                for calc in calcs]


def _print_suite_summary(calc_suite_specs):
    """Print summary of requested calculations."""
    return ('\nRequested aospy calculations:\n' +
            pprint.pformat(calc_suite_specs) + '\n')


def submit_mult_calcs(calc_suite_specs, exec_options=None):
    """Generate and execute all specified computations.

    Once the calculations are prepped and submitted for execution, any
    calculation that triggers any exception or error is skipped, and the rest
    of the calculations proceed unaffected.  This prevents an error in a single
    calculation from crashing a large suite of calculations.

    Parameters
    ----------
    calc_suite_specs : dict
        The specifications describing the full set of calculations to be
        generated and potentially executed.  Accepted keys and their values:

        library : module or package comprising an aospy object library
            The aospy object library for these calculations.
        projects : list of aospy.Proj objects
            The projects to permute over.
        models : 'all', 'default', or list of aospy.Model objects
            The models to permute over.  If 'all', use all models in the
            ``models`` attribute of each ``Proj``.  If 'default', use all
            models in the ``default_models`` attribute of each ``Proj``.
        runs : 'all', 'default', or list of aospy.Run objects
            The runs to permute over.  If 'all', use all runs in the
            ``runs`` attribute of each ``Model``.  If 'default', use all
            runs in the ``default_runs`` attribute of each ``Model``.
        variables : list of aospy.Var objects
            The variables to be calculated.
        regions : 'all' or list of aospy.Region objects
            The region(s) over which any regional reductions will be performed.
            If 'all', use all regions in the ``regions`` attribute of each
            ``Proj``.
        date_ranges : 'default' or tuple of datetime.datetime objects
            The range of dates (inclusive) over which to perform calculations.
            If 'default', use the ``default_start_date`` and
            ``default_end_date`` attribute of each ``Run``.
        output_time_intervals : {'ann', season-string, month-integer}
            The sub-annual time interval over which to aggregate.

            - 'ann' : Annual mean
            - season-string : E.g. 'JJA' for June-July-August
            - month-integer : 1 for January, 2 for February, etc.  Each one is
                  a separate reduction, e.g. [1, 2] would produce averages (or
                  other specified time reduction) over all Januaries, and
                  separately over all Februaries.

        output_time_regional_reductions : list of reduction string identifiers
            Unlike most other keys, these are not permuted over when creating
            the :py:class:`aospy.Calc` objects that execute the calculations;
            each :py:class:`aospy.Calc` performs all of the specified
            reductions.  Accepted string identifiers are:

            - Gridpoint-by-gridpoint output:

              - 'av' : Gridpoint-by-gridpoint time-average
              - 'std' : Gridpoint-by-gridpoint temporal standard deviation
              - 'ts' : Gridpoint-by-gridpoint time-series

            - Averages over each region specified via `region`:

              - 'reg.av', 'reg.std', 'reg.ts' : analogous to 'av', 'std', 'ts'

        output_vertical_reductions : {None, 'vert_av', 'vert_int'}, optional
            How to reduce the data vertically:

            - None : no vertical reduction
            - 'vert_av' : mass-weighted vertical average
            - 'vert_int' : mass-weighted vertical integral
        input_time_intervals : {'annual', 'monthly', 'daily', '#hr'}
            A string specifying the time resolution of the input data.  In
            '#hr' above, the '#' stands for a number, e.g. 3hr or 6hr, for
            sub-daily output.  These are the suggested specifiers, but others
            may be used if they are also used by the DataLoaders for the given
            Runs.
        input_time_datatypes : {'inst', 'ts', 'av'}
            What the time axis of the input data represents:

            - 'inst' : Timeseries of instantaneous values
            - 'ts' : Timeseries of averages over the period of each time-index
            - 'av' : A single value averaged over a date range

        input_vertical_datatypes : {False, 'pressure', 'sigma'}, optional
            The vertical coordinate system used by the input data:

            - False : not defined vertically
            - 'pressure' : pressure coordinates
            - 'sigma' : hybrid sigma-pressure coordinates

        input_time_offsets : {None, dict}, optional
            How to offset input data in time to correct for metadata errors

            - None : no time offset applied
            - dict : e.g. ``{'hours': -3}`` to offset times by -3 hours
              See :py:meth:`aospy.utils.times.apply_time_offset`.

    exec_options : dict or None (default None)
        Options regarding how the calculations are reported, submitted, and
        saved.  If None, default settings are used for all options.  Currently
        supported options (each should be either `True` or `False`):

        - prompt_verify : (default False) If True, print summary of
              calculations to be performed and prompt user to confirm before
              submitting for execution.
        - parallelize : (default False) If True, submit calculations in
              parallel.  This requires the `multiprocess` library, which can be
              installed via `pip install multiprocess`.
        - write_to_tar : (default True) If True, write results of calculations
              to .tar files, one for each :py:class:`aospy.Run` object.  These tar files have an
              identical directory structures the standard output relative to
              their root directory, which is specified via the `tar_direc_out`
              argument of each Proj object's instantiation.

    Returns
    -------
    A list of the return values from each :py:meth:`aospy.Calc.compute` call
        If a calculation ran without error, this value is the
        :py:class:`aospy.Calc` object itself, with the results of its
        calculations saved in its ``data_out`` attribute.  ``data_out`` is a
        dictionary, with the keys being the temporal-regional reduction
        identifiers (e.g. 'reg.av'), and the values being the corresponding
        result.

    If any error occurred during a calculation, the return value is None.

    Raises
    ------
    AospyException
        If the ``prompt_verify`` option is set to True and the user does not
        respond affirmatively to the prompt.

    """
    if exec_options is None:
        exec_options = dict()
    if exec_options.pop('prompt_verify', False):
        print(_print_suite_summary(calc_suite_specs))
        _user_verify()
    calc_suite = CalcSuite(calc_suite_specs)
    calcs = calc_suite.create_calcs()
    return _exec_calcs(calcs, **exec_options)
