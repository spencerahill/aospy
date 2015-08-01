#!/usr/bin/env python
"""Main test class for aospy."""

import sys
import unittest

from aospy import Constant

CONST_VALUE = 1
CONST_UNITS = 'dummy_units'
CONST_DESCRIPTION = 'dummy_description'
CONST2_VALUE = 2.
CONST2_UNITS = 'dummy_units2'


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


class AospyRunTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


class TestExp(AospyRunTestCase):
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
        self.const = Constant(CONST_VALUE, CONST_UNITS,
                              description=CONST_DESCRIPTION)
        self.const2 = Constant(CONST2_VALUE, CONST2_UNITS)

    def tearDown(self):
        pass


class TestConstant(AospyConstantTestCase):
    def test_init(self):
        self.assertEqual(self.const.value, CONST_VALUE)
        self.assertEqual(self.const.units, CONST_UNITS)
        self.assertEqual(self.const.description, CONST_DESCRIPTION)
        self.assertEqual(self.const2.description, '')

    def test_add_two_consts(self):
        self.assertEqual(CONST_VALUE + CONST_VALUE, self.const + self.const)

    def test_add_const_scalar(self):
        self.assertEqual(CONST_VALUE + CONST_VALUE, self.const + CONST_VALUE)
        self.assertEqual(CONST_VALUE + CONST_VALUE, CONST_VALUE + self.const)

    def test_add_two_consts_units_mismatch(self):
        self.assertRaises(TypeError, self.const.__add__, self.const2)
        self.assertRaises(TypeError, self.const.__radd__, self.const2)

    def test_subtract_two_consts(self):
        self.assertEqual(CONST_VALUE - CONST_VALUE, self.const - self.const)

    def test_subtract_const_scalar(self):
        self.assertEqual(CONST_VALUE - CONST_VALUE, self.const - CONST_VALUE)
        self.assertEqual(CONST_VALUE - CONST_VALUE, CONST_VALUE - self.const)

    def test_subtract_two_consts_units_mismatch(self):
        self.assertRaises(TypeError, self.const.__sub__, self.const2)
        self.assertRaises(TypeError, self.const.__rsub__, self.const2)

    def test_multiply_two_consts(self):
        self.assertEqual(CONST_VALUE * CONST_VALUE, self.const * self.const)

    def test_multiply_const_scalar(self):
        self.assertEqual(CONST_VALUE * CONST_VALUE, self.const * CONST_VALUE)
        self.assertEqual(CONST_VALUE * CONST_VALUE, CONST_VALUE * self.const)

    def test_divide_two_consts(self):
        self.assertEqual(CONST_VALUE / CONST_VALUE, self.const / self.const)

    def test_divide_const_scalar(self):
        self.assertEqual(CONST_VALUE / CONST_VALUE, self.const / CONST_VALUE)
        self.assertEqual(CONST_VALUE / CONST_VALUE, CONST_VALUE / self.const)

if __name__ == '__main__':
    sys.exit(unittest.main())
