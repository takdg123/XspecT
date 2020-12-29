import unittest, os, shutil, sys
import numpy as np
from gbm.background.binned import polynomial
py_version = sys.version_info[0]


counts = np.array([78.0, 58.0, 40.0, 26.0, 14.0, 6.0, 2.0, 0.0, 2.0, 6.0])
exposure = np.array([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
edges = np.array([-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
bounds = (edges[0:-1], edges[1:])

class TestPolynomialBackground(unittest.TestCase):
    
    def test_fit(self): 
        bkgd = polynomial(2)
        rates, errs = bkgd.fit(counts, exposure, bounds)
        test_arr = np.array([39.0, 29.0, 20.0, 13.0, 7.0, 3.0, 1.0, 0.0, 1.0, 3.0])
        for i in range(10):
            self.assertAlmostEqual(rates[i], test_arr[i], places=0)
        
        self.assertRaises(AssertionError, polynomial, -1)
        self.assertRaises(AssertionError, bkgd.fit, counts, exposure[:-1], edges)
        self.assertRaises(AssertionError, bkgd.fit, counts, exposure, edges[:-1])


    def test_interpolate(self):
        bkgd = polynomial(2)
        _, _ = bkgd.fit(counts, exposure, bounds)
        rates, errs = bkgd.interpolate(edges[:-1], edges[1:])    
        test_arr = np.array([39.0, 29.0, 20.0, 13.0, 7.0, 3.0, 1.0, 0.0, 1.0, 3.0])
        for i in range(10):
            self.assertAlmostEqual(rates[i], test_arr[i], places=0)
