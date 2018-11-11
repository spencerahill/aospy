from multiprocessing import cpu_count
from os.path import isfile
import shutil
import sys
import itertools

import distributed
import pytest

from aospy import Var, Proj
from aospy.automate import (_get_attr_by_tag, _permuted_dicts_of_specs,
                            _get_all_objs_of_type, _merge_dicts,
                            _input_func_py2_py3, AospyException,
                            _user_verify, CalcSuite, _MODELS_STR, _RUNS_STR,
                            _VARIABLES_STR, _REGIONS_STR,
                            _compute_or_skip_on_error, submit_mult_calcs,
                            _n_workers_for_local_cluster,
                            _prune_invalid_time_reductions)
from .data.objects import examples as lib
from .data.objects.examples import (
    example_proj, example_model, example_run, var_not_time_defined,
    condensation_rain, convection_rain, precip, ps, sphum, globe, sahel, bk,
    p, dp
)


@pytest.fixture
def obj_lib():
    return lib


@pytest.fixture
def all_vars():
    return [condensation_rain, convection_rain, precip, ps, sphum]


@pytest.fixture
def all_projects():
    return [example_proj]


@pytest.fixture
def all_models():
    return [example_model]


@pytest.fixture
def all_runs():
    return [example_run]


@pytest.fixture
def all_regions():
    return [globe, sahel]


@pytest.mark.parametrize(
    ('obj', 'tag', 'attr_name', 'expected'),
    [(example_proj, 'all', _MODELS_STR, [example_model]),
     (example_proj, 'default', _MODELS_STR, []),
     (example_model, 'all', _RUNS_STR, [example_run]),
     (example_model, 'default', _RUNS_STR, [])])
def test_get_attr_by_tag(obj, tag, attr_name, expected):
    actual = _get_attr_by_tag(obj, tag, attr_name)
    assert actual == expected


def test_get_attr_by_tag_invalid():
    with pytest.raises(KeyError):
        _get_attr_by_tag(example_proj, 'alll', _MODELS_STR)


@pytest.fixture
def calcsuite_specs():
    """Aux specs after being processed by CalcSuite."""
    return {
        'time_offset': [None],
        'date_range': ['default'],
        'intvl_in': ['monthly'],
        'region': [{globe, sahel}],
        'dtype_out_time': [['av', 'reg.av']],
        'dtype_in_vert': [False],
        'dtype_in_time': ['ts'],
        'var': [condensation_rain, convection_rain],
        'intvl_out': ['ann'],
        'dtype_out_vert': [None]
    }


def test_permuted_dict_of_specs(calcsuite_specs):
    actual = _permuted_dicts_of_specs(calcsuite_specs)
    expected = [
        {'time_offset': None,
         'date_range': 'default',
         'intvl_in': 'monthly',
         'region': {globe, sahel},
         'dtype_out_time': ['av', 'reg.av'],
         'dtype_in_vert': False,
         'dtype_in_time': 'ts',
         'var': condensation_rain,
         'intvl_out': 'ann',
         'dtype_out_vert': None},
        {'time_offset': None,
         'date_range': 'default',
         'intvl_in': 'monthly',
         'region': {globe, sahel},
         'dtype_out_time': ['av', 'reg.av'],
         'dtype_in_vert': False,
         'dtype_in_time': 'ts',
         'var': convection_rain,
         'intvl_out': 'ann',
         'dtype_out_vert': None}
    ]
    assert actual == expected


def test_merge_dicts():
    # no conflicts
    dict1 = dict(a=1)
    dict2 = {'b': 3, 43: False}
    dict3 = dict(c=['abc'])
    expected = {'a': 1, 'b': 3, 'c': ['abc'], 43: False}
    assert expected == _merge_dicts(dict1, dict2, dict3)

    # conflicts
    dict4 = dict(c=None)
    expected = {'a': 1, 'b': 3, 'c': None, 43: False}
    assert expected == _merge_dicts(dict1, dict2, dict3, dict4)


def test_input_func_py2_py3():
    result = _input_func_py2_py3()
    if sys.version.startswith('3'):
        import builtins
        assert result is builtins.input
    elif sys.version.startswith('2'):
        assert result is raw_input  # noqa: F821


def test_user_verify():
    with pytest.raises(AospyException):
        _user_verify(lambda x: 'no')
    _user_verify(lambda x: 'YES')


@pytest.mark.parametrize(
    ('type_', 'expected'),
    [(Var, [var_not_time_defined, condensation_rain, convection_rain,
            precip, ps, sphum, bk, p, dp]),
     (Proj, [example_proj])])
def test_get_all_objs_of_type(obj_lib, type_, expected):
    actual = _get_all_objs_of_type(type_, obj_lib)
    assert set(expected) == set(actual)


@pytest.fixture
def calcsuite_init_specs():
    return dict(
        library=lib,
        projects=[example_proj],
        models=[example_model],
        runs=[example_run],
        variables=[condensation_rain, convection_rain],
        regions='all',
        date_ranges='default',
        output_time_intervals=['ann'],
        output_time_regional_reductions=['av', 'reg.av'],
        output_vertical_reductions=[None],
        input_time_intervals=['monthly'],
        input_time_datatypes=['ts'],
        input_time_offsets=[None],
        input_vertical_datatypes=[False],
    )


@pytest.fixture
def calcsuite_init_specs_single_calc(calcsuite_init_specs):
    specs = calcsuite_init_specs.copy()
    specs['variables'] = [condensation_rain]
    specs['regions'] = [None]
    specs['output_time_regional_reductions'] = ['av']
    yield specs
    # Teardown procedure
    for direc in [example_proj.direc_out, example_proj.tar_direc_out]:
        shutil.rmtree(direc, ignore_errors=True)


@pytest.fixture
def calcsuite_init_specs_two_calcs(calcsuite_init_specs):
    specs = calcsuite_init_specs.copy()
    specs['variables'] = [condensation_rain, convection_rain]
    specs['regions'] = [None]
    specs['output_time_regional_reductions'] = ['av']
    yield specs
    # Teardown procedure
    for direc in [example_proj.direc_out, example_proj.tar_direc_out]:
        shutil.rmtree(direc, ignore_errors=True)


@pytest.fixture
def calc(calcsuite_init_specs_single_calc):
    return CalcSuite(calcsuite_init_specs_single_calc).create_calcs()[0]


def test_compute_or_skip_on_error(calc, caplog):
    result = _compute_or_skip_on_error(calc, dict(write_to_tar=False))
    assert result is calc

    calc.start_date = 'dummy'
    result = _compute_or_skip_on_error(calc, dict(write_to_tar=False))
    log_record = caplog.record_tuples[-1][-1]
    assert log_record.startswith("Skipping aospy calculation")
    assert result is None


@pytest.fixture
def external_client():
    # Explicitly specify we want only 4 workers so that when running on
    # Travis we don't request too many.
    cluster = distributed.LocalCluster(n_workers=4)
    client = distributed.Client(cluster)
    yield client
    client.close()
    cluster.close()


def assert_calc_files_exist(calcs, write_to_tar, dtypes_out_time):
    """Check that expected calcs were written to files"""
    for calc in calcs:
        for dtype_out_time in dtypes_out_time:
            assert isfile(calc.path_out[dtype_out_time])
            if write_to_tar:
                assert isfile(calc.path_tar_out)
            else:
                assert not isfile(calc.path_tar_out)


@pytest.mark.skipif(sys.version.startswith('2'),
                    reason='https://github.com/spencerahill/aospy/issues/259')
@pytest.mark.parametrize(
    ('exec_options'),
    [dict(parallelize=True, write_to_tar=False),
     dict(parallelize=True, write_to_tar=True)])
def test_submit_mult_calcs_external_client(calcsuite_init_specs_single_calc,
                                           external_client, exec_options):
    exec_options.update(client=external_client)
    calcs = submit_mult_calcs(calcsuite_init_specs_single_calc, exec_options)
    write_to_tar = exec_options.pop('write_to_tar', True)
    assert_calc_files_exist(
        calcs, write_to_tar,
        calcsuite_init_specs_single_calc['output_time_regional_reductions'])


@pytest.mark.skipif(sys.version.startswith('2'),
                    reason='https://github.com/spencerahill/aospy/issues/259')
@pytest.mark.parametrize(
    ('exec_options'),
    [dict(parallelize=False, write_to_tar=False),
     dict(parallelize=True, write_to_tar=False),
     dict(parallelize=False, write_to_tar=True),
     dict(parallelize=True, write_to_tar=True),
     None])
def test_submit_mult_calcs(calcsuite_init_specs_single_calc, exec_options):
    calcs = submit_mult_calcs(calcsuite_init_specs_single_calc, exec_options)
    if exec_options is None:
        write_to_tar = True
    else:
        write_to_tar = exec_options.pop('write_to_tar', True)
    assert_calc_files_exist(
        calcs, write_to_tar,
        calcsuite_init_specs_single_calc['output_time_regional_reductions'])


def test_submit_mult_calcs_no_calcs(calcsuite_init_specs):
    specs = calcsuite_init_specs.copy()
    specs['input_vertical_datatypes'] = []
    with pytest.raises(AospyException):
        submit_mult_calcs(specs)


@pytest.mark.skipif(sys.version.startswith('2'),
                    reason='https://github.com/spencerahill/aospy/issues/259')
@pytest.mark.parametrize(
    ('exec_options'),
    [dict(parallelize=True, write_to_tar=False),
     dict(parallelize=True, write_to_tar=True)])
def test_submit_two_calcs_external_client(calcsuite_init_specs_two_calcs,
                                          external_client, exec_options):
    exec_options.update(client=external_client)
    calcs = submit_mult_calcs(calcsuite_init_specs_two_calcs, exec_options)
    write_to_tar = exec_options.pop('write_to_tar', True)
    assert_calc_files_exist(
        calcs, write_to_tar,
        calcsuite_init_specs_two_calcs['output_time_regional_reductions'])


@pytest.mark.skipif(sys.version.startswith('2'),
                    reason='https://github.com/spencerahill/aospy/issues/259')
@pytest.mark.parametrize(
    ('exec_options'),
    [dict(parallelize=False, write_to_tar=False),
     dict(parallelize=True, write_to_tar=False),
     dict(parallelize=False, write_to_tar=True),
     dict(parallelize=True, write_to_tar=True),
     None])
def test_submit_two_calcs(calcsuite_init_specs_two_calcs, exec_options):
    calcs = submit_mult_calcs(calcsuite_init_specs_two_calcs, exec_options)
    if exec_options is None:
        write_to_tar = True
    else:
        write_to_tar = exec_options.pop('write_to_tar', True)
    assert_calc_files_exist(
        calcs, write_to_tar,
        calcsuite_init_specs_two_calcs['output_time_regional_reductions'])


def test_n_workers_for_local_cluster(calcsuite_init_specs_two_calcs):
    calcs = CalcSuite(calcsuite_init_specs_two_calcs).create_calcs()
    expected = min(cpu_count(), len(calcs))
    result = _n_workers_for_local_cluster(calcs)
    assert result == expected


@pytest.fixture
def calc_suite(calcsuite_init_specs):
    return CalcSuite(calcsuite_init_specs)


class TestCalcSuite(object):
    def test_init(self, calc_suite, calcsuite_init_specs, obj_lib):
        assert calc_suite._specs_in == calcsuite_init_specs
        assert calc_suite._obj_lib == obj_lib

    def test_permute_core_specs(self, calc_suite):
        expected = [dict(proj=example_proj, model=example_model,
                         run=example_run)]
        actual = calc_suite._permute_core_specs()
        assert expected == actual
        # TODO: cases w/ multiple projs and/or models and/or runs, with
        #       different default children for each

    def test_get_regions(self, calc_suite, all_regions):
        assert calc_suite._get_regions()[0] == set(all_regions)
        # TODO: case w/ not all regions
        # TODO: case w/ Region objects in 'regions' sub-module

    def test_get_variables(self, calc_suite, all_vars):
        assert not hasattr(calc_suite, 'variables')
        assert calc_suite._get_variables() == {condensation_rain,
                                               convection_rain}
        # TODO: case w/ Var objects in 'variables' sub-module
        # TODO: case w/ 'all'

    def test_get_aux_specs(self, calc_suite, all_regions):
        spec_names = [name for name in calc_suite._AUX_SPEC_NAMES
                      if name not in [_VARIABLES_STR, _REGIONS_STR]]
        expected = {name: calc_suite._specs_in[name] for name in spec_names}
        expected[_VARIABLES_STR] = {condensation_rain, convection_rain}
        expected[_REGIONS_STR] = [{globe, sahel}]
        expected['date_ranges'] = ['default']
        expected['output_time_regional_reductions'] = [['av', 'reg.av']]
        actual = calc_suite._get_aux_specs()
        assert actual == expected

    def test_permute_aux_specs(self, calc_suite, calcsuite_specs):
        expected = _permuted_dicts_of_specs(calcsuite_specs)
        actual = calc_suite._permute_aux_specs()
        assert len(actual) == len(expected)
        for act in actual:
            assert act in expected


@pytest.mark.parametrize('var', [var_not_time_defined, condensation_rain])
def test_prune_invalid_time_reductions(var):
    time_options = ['av', 'std', 'ts', 'reg.av', 'reg.std', 'reg.ts']
    spec = {
        'var': var,
        'dtype_out_time': None
    }
    assert _prune_invalid_time_reductions(spec) is None
    for i in range(1, len(time_options) + 1):
        for time_option in list(itertools.permutations(time_options, i)):
            spec['dtype_out_time'] = time_option
            if spec['var'].def_time:
                assert _prune_invalid_time_reductions(spec) == time_option
            else:
                assert _prune_invalid_time_reductions(spec) == []
