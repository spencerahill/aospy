#!/usr/bin/env python
"""Main test class for aospy."""

import sys
import unittest

import numpy as np

from aospy import Constant


class AospyProjTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestProj(AospyProjTestCase):
    pass


class AospyModelTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestModel(AospyModelTestCase):
    pass


class AospyVarTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestVar(AospyVarTestCase):
    pass


class AospyCalcTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestCalc(AospyCalcTestCase):
    pass


class AospyConstantTestCase(unittest.TestCase):
    def setUp(self):
        self._val1 = 1
        self._const1_units = 'dummy_units'
        self._const1_description = 'dummy_description'
        self._val2 = 2.
        self._const2_units = 'dummy_units2'

        self.const1 = Constant(self._val1, self._const1_units,
                               description=self._const1_description)
        self.const2 = Constant(self._val2, self._const2_units)

    def tearDown(self):
        pass


class TestConstant(AospyConstantTestCase):
    def test_init(self):
        self.assertEqual(self.const1.value, self._val1)
        self.assertEqual(self.const1.units, self._const1_units)
        self.assertEqual(self.const1.description, self._const1_description)
        self.assertEqual(self.const2.description, '')

    def test_add_two_consts(self):
        self.assertEqual(self._val1 + self._val1, self.const1 + self.const1)

    def test_add_const_scalar(self):
        self.assertEqual(self._val1 + self._val1, self.const1 + self._val1)
        self.assertEqual(self._val1 + self._val1, self._val1 + self.const1)

    def test_add_const_numpy_array(self):
        self.assertEqual(np.array(self._val1) + self.const1,
                         self._val1 + self.const1.value)

    def test_add_two_consts_units_mismatch(self):
        self.assertRaises(TypeError, self.const1.__add__, self.const2)
        self.assertRaises(TypeError, self.const1.__radd__, self.const2)

    def test_subtract_two_consts(self):
        self.assertEqual(self._val1 - self._val1, self.const1 - self.const1)

    def test_subtract_const_scalar(self):
        self.assertEqual(self._val1 - self._val1, self.const1 - self._val1)
        self.assertEqual(self._val1 - self._val1, self._val1 - self.const1)

    # TODO: Fix test and/or the code it is testing.
    @unittest.expectedFailure
    def test_subtract_two_consts_units_mismatch(self):
        self.assertRaises(TypeError, self.const1.__sub__, self.const2)
        self.assertRaises(TypeError, self.const1.__rsub__, self.const2)

    def test_multiply_two_consts(self):
        self.assertEqual(self._val1 * self._val2, self.const1 * self.const2)

    def test_multiply_const_scalar(self):
        self.assertEqual(self._val1 * self._val2, self.const1 * self._val2)
        self.assertEqual(self._val1 * self._val2, self._val1 * self.const2)

    @unittest.skip('TODO: Fix test and/or the code it is testing')
    def test_divide_two_consts(self):
        self.assertEqual(self._val1 / self._val2, self.const1 / self.const2)

    @unittest.skip('TODO: Fix test and/or the code it is testing')
    def test_divide_const_scalar(self):
        self.assertEqual(self._val1 / self._val2, self.const1 / self._val2)
        self.assertEqual(self._val1 / self._val2, self._val1 / self.const2)

    def test_power_two_consts(self):
        self.assertEqual(self._val1**self._val2, self.const1**self.const2)

    def test_power_const_scalar(self):
        self.assertEqual(self._val1**self._val2, self.const1**self._val2)
        self.assertEqual(self._val1**self._val2, self._val1**self.const2)


if __name__ == '__main__':
    sys.exit(unittest.main())
