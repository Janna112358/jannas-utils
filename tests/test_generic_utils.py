#!/usr/bin/env python3

import unittest
import jannasutils as ju
import numpy as np

class Test_isIterable(unittest.TestCase):
    
    def test_scalar(self):
        self.assertFalse(ju.isIterable(1))
        self.assertFalse(ju.isIterable(1.))
    
    def test_array_like(self):
        self.assertTrue(ju.isIterable([1, 2]))
        self.assertTrue(ju.isIterable(np.array([1, 2])))
