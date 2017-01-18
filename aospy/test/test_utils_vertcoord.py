#!/usr/bin/env python
"""Test suite for aospy.utils module."""
import sys
import unittest

import numpy as np

import aospy.utils.vertcoord as vertcoord


class AospyUtilsTestCase(unittest.TestCase):
    def setUp(self):
        self.p_in_hpa = np.array([1000, 925, 850, 775, 700, 600, 500, 400, 300,
                                  200, 150, 100, 70, 50, 30, 20, 10],
                                 dtype=np.float64)
        self.p_in_pa = self.p_in_hpa*1e2
        self.p_top = 0
        self.p_bot = 1.1e5
        self.p_edges = 0.5*(self.p_in_pa[1:] + 0.5*self.p_in_pa[:-1])
        self.phalf = np.concatenate(([self.p_bot], self.p_edges, [self.p_top]))

    def tearDown(self):
        pass


class TestUtils(AospyUtilsTestCase):
    def test_to_pascal_scalar_positive(self):
        self.assertEqual(vertcoord.to_pascal(1e2), 1e4)
        self.assertEqual(vertcoord.to_pascal(1e5), 1e5)

    def test_to_pascal_scalar_negative(self):
        self.assertEqual(vertcoord.to_pascal(-1e2), -1e4)
        self.assertEqual(vertcoord.to_pascal(-1e5), -1e5)

    def test_to_pascal_array(self):
        np.testing.assert_array_equal(vertcoord.to_pascal(self.p_in_hpa),
                                      self.p_in_pa)
        np.testing.assert_array_equal(vertcoord.to_pascal(self.p_in_pa),
                                      self.p_in_pa)


if __name__ == '__main__':
    sys.exit(unittest.main())
