#!/usr/bin/env python
import os
import unittest
import numpy as np
import aospy

class AospyProjTestCase(unittest.TestCase):
    def delete_files_in_dir(self,direc, string):
        """In a directory, delete files whose names contain a given string."""
        for fl in os.listdir(direc):
            if string in fl:
                file_str = (direc + '/' + fl).replace('//', '/')
                os.remove(file_str)

    def setUp(self):
        self.Proj = gen_apc524()
        self.model = 'am2'
        self.Model = self.Proj.models[self.model]
        self.exp = 'control'
        self.Exp = self.Model.exps[self.exp]
        self.var = 't_surf'
        self.Var = self.Exp.vars[self.var]
        # self.delete_files_in_dir(aospy.path_in, 'dummy')

    def tearDown(self):
        """Several tests create files in the aospy.path_in directory with the
        name 'dummy' in them.  So delete any that remain.

        """
        self.delete_files_in_dir(aospy.path_in, 'dummy')


class TestProj(AospyProjTestCase):
    def test_init(self):
        self.assertEqual(self.Proj.name, 'apc524')
        self.assertEqual(self.Proj.type, 'Proj')
        self.assertIsNone(self.Proj.parent)
        self.assertEqual(self.Proj.path_in, aospy.path_in)
        self.assertEqual(self.Proj.path_out, aospy.path_out)

    def test_get_child(self):
        self.assertIs(self.Proj._get_child(model=self.model, exp=self.exp,
                                           var=self.var), self.Var)
        self.assertIs(self.Proj._get_child(obj_type='exp', model=self.model,
                                           exp=self.exp), self.Exp)
        self.assertIs(self.Proj._get_child(obj_type='model',
                                           model=self.model), self.Model)
        self.assertRaises(Exception, self.Proj._get_child, model='dummy')
        self.assertRaises(Exception, self.Proj._get_child, model=self.model,
                          exp='dummy')
        self.assertRaises(Exception, self.Proj._get_child, exp=self.exp,
                          var='dummy')

    def test_get_nc_bad_var_in(self):
        self.assertRaises(Exception, self.Proj._get_nc, 'dummy', 'dummy',
                          'dummy')

    def test_get_nc_bad_path_in(self):
        self.Proj.path_in = 'dummy'
        self.assertRaises(Exception, self.Proj._get_nc, self.model, self.exp,
                          self.var)

    def test_get_nc_empty_files_list(self):
        """Create dummy variable in the input netCDF directory to test that
        Proj._get_nc catches that the input list of netCDF data it
        generates is empty.
        """
        dummyvar = aospy.Var('dummy', self.Exp, None)
        self.assertRaises(Exception, self.Proj._get_nc, self.model,
                          self.exp, 'dummy')

    def test_get_nc_bad_files_list(self):
        """Create dummy variable and corresponding empty file in the input
        netCDF directory to test that Proj._get_nc catches that the
        generated list of netCDF data lacks any valid netCDF files.

        """
        dummyvar = aospy.Var('dummy', self.Exp, None)
        path = aospy.path_in + '/atmos.198901-198912.dummy.nc'
        with open(path, 'w') as filein:
            self.assertRaises(Exception, self.Proj._get_nc, self.model,
                              self.exp, 'dummy')
        os.remove(path)

    def test_get_nc_bad_file_in_list(self):
        """Create dummy blank file that includes the name of a valid var in
        order to test that Proj._get_nc catches the invalid file and
        ignores it.

        """
        nc_orig = self.Proj._get_nc(self.model, self.exp, self.var)
        path = aospy.path_in + '/dummy_' + self.var
        with open(path, 'w') as filein:
            nc_new = self.Proj._get_nc(self.model, self.exp, self.var)
            self.assertEqual(nc_orig._files, nc_new._files)
        os.remove(path)
        

class AospyModelTestCase(unittest.TestCase):
    def setUp(self):
        self.Model = aospy.Model('testmodel', None)

class TestModel(AospyModelTestCase):
    def test_init(self):
        self.assertEqual(self.Model.name, 'testmodel')
        self.assertEqual(self.Model.type, 'Model')
        self.assertEqual(self.Model.exps, {})
        self.assertEqual(self.Model.vars, {})        
        self.assertIsNone(self.Model.parent)


class AospyExpTestCase(unittest.TestCase):
    def setUp(self):
        self.Exp = aospy.Exp('testexp', None)

class TestExp(AospyExpTestCase):
    def test_init(self):
        self.assertEqual(self.Exp.name, 'testexp')
        self.assertEqual(self.Exp.type, 'Exp')
        self.assertEqual(self.Exp.vars, {})
        self.assertIsNone(self.Exp.parent)


class AospyVarTestCase(unittest.TestCase):
    def setUp(self):
        self.Var = aospy.Var('testvar', None, None)

class TestVar(AospyVarTestCase):
    def test_init(self):
        self.assertEqual(self.Var.name, 'testvar')
        self.assertEqual(self.Var.type, 'Var')
        self.assertIsNone(self.Var.parent)
    
    def test_set_func_args_no_func(self):
        self.assertIsNone(self.Var.func_args)

    def test_set_func_args_w_func(self):
        def dummy_func(arg):
            return arg
        dummy_var = aospy.Var('dummyvar', None, dummy_func, func_args='dummy')
        self.assertEqual(dummy_var.func_args, 'dummy')


class AospyCalcTestCase(unittest.TestCase):
    def setUp(self):
        self.Proj = gen_apc524()
        self.Calc = aospy.Calc(self.Proj, 'am2', 'control', 't_surf')        

class TestCalc(AospyCalcTestCase):
    def test_specs_accept(self):
        """Test if specs are passed in properly."""
        start = [2002,8,16,12]
        end = [2012,9,16,00]
        self.Calc.specs.update({'date_range': [start, end]})
        desired_range = [start, end]
        actual_range = self.Calc.specs['date_range']
        self.assertEqual(desired_range, actual_range)

    def test_specs_reject(self):
        """Test if catches error, returns all data when invalid date range is
        given.

        """
        start = [2002,7,16,14] # invalid dates
        end = [2008,0,0,0]
        self.Calc.specs.update({'date_range': [start, end]})
        all_data = self.Proj._get_nc(
            'am2','control','t_surf').variables['t_surf']
        actual_data = self.Calc._data
        self.assertEqual(len(all_data), len(actual_data))

    def test_metadata(self): 
        """Test if metadata are generated correctly."""
        model = 'am2'
        experiment = 'control'
        variable = 't_surf'
        desired_metadata = {'model': model, 'experiment': experiment,
                            'variable': variable, 'specifications': {}}
        self.assertEqual(desired_metadata, self.Calc._metadata)
    
    def test_exec_func_is_none(self):
        """Test Calc._exec_func when the function is None."""
        result = self.Calc._exec_func()
        np.testing.assert_equal(self.Calc._data, result)

    def test_exec_func_arg_is_self(self):
        """Test Calc._exec_func when func_args includes the Var itself."""
        def dummy_func(arg):
            return 3.*arg
        self.Calc.Var.function = dummy_func
        self.Calc.Var.func_args = ['t_surf'] # same as Calc.Var.name
        self.Calc._data = self.Calc._get_data()
        result = self.Calc._exec_func()
        np.testing.assert_equal(self.Calc._data*3., result)

    def test_exec_func_with_args(self):
        """Test Var that is a function of other Vars is computed correctly."""
        Calc2 = aospy.Calc(self.Proj, 'am2', 'control', 'ts_pr_wtd')
        result = Calc2._exec_func()
        self.assertEqual(np.shape(result), (24, 90, 144))

    def test_find_reduct_method(self):
        """Test if the correct reduction method is found."""
        self.Calc.time_reduct = 'stdev'
        result = self.Calc._apply_time_reduct(self.Calc._data)
        np.testing.assert_equal(result, aospy.reductions.stdev(self.Calc._data))

    def test_default_reduct(self):
        """Test if 'avg' time reduction method is used if none specified."""
        result = self.Calc._apply_time_reduct(self.Calc._data)
        np.testing.assert_equal(result, aospy.reductions.avg(self.Calc._data))

    def test_invalid_reduct(self):
        """Test if 'avg' method is used when invalid method is specified."""
        self.Calc.time_reduct = 'dummy'
        result = self.Calc._apply_time_reduct(self.Calc._data)
        np.testing.assert_equal(result, aospy.reductions.avg(self.Calc._data))

    def test_result_array(self):
        self.Calc.set_result_array()
        result = self.Calc._data
        time_red_result = aospy.reductions.avg(result)
        np.testing.assert_equal(time_red_result, self.Calc._result_array)

class AospyTestCase(unittest.TestCase):
    def setUp(self):
        self.Proj = gen_apc524()
        path = aospy.path_out + '/testDatabase'
        self.projDB = aospy.AospyDB(path, 'AospyDB')

    def tearDown(self):
        self.projDB.delete_all()
        self.projDB.close()

class TestAospy(AospyTestCase):
    def test_get_attr_up(self):
        self.assertRaises(Exception, aospy._get_attr_up, None, 'dummy')
        self.assertRaises(Exception, aospy._get_attr_up, self.Proj, 'dummy')

        model = self.Proj.models['am2']
        exp = model.exps['control']
        var = exp.vars['t_surf']

        self.assertEqual(aospy._get_attr_up(var, 'name'), 't_surf')
        self.assertEqual(
            aospy._get_attr_up(var, 'path_in'), aospy.path_in)
        self.assertEqual(
            aospy._get_attr_up(var, 'path_in'), aospy.path_in)
    
    def test_update_child_dict(self):
        name = 'dummy'
        dummy_model = aospy.Model(name, None)
        aospy._update_child_dict(self.Proj, dummy_model, 'models')
        self.assertIs(self.Proj.models[name], dummy_model)

        dummy_exp = aospy.Exp(name, None)
        aospy._update_child_dict(dummy_model, dummy_exp, 'exps')
        self.assertIs(dummy_model.exps[name], dummy_exp)

        dummy_var = aospy.Var(name, None, None)
        aospy._update_child_dict(dummy_exp, dummy_var, 'vars')
        self.assertIs(dummy_exp.vars[name], dummy_var)

    def test_set_parent_no_update(self):
        aospy._set_parent(self.Proj, 'dummy', update_child_dict=False)
        self.assertEqual(self.Proj.parent, 'dummy')

    def test_set_parent_do_update(self):
        dummyProj = aospy.Proj('dummyProj', None)
        dummyModel = aospy.Model('dummyModel', None)
        aospy._set_parent(dummyModel, dummyProj, update_child_dict=True)
        self.assertIs(dummyProj.models['dummyModel'], dummyModel)


if __name__ == '__main__':
    unittest.main()
