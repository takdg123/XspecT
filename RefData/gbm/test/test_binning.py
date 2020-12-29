import unittest, os, shutil, sys
import numpy as np
from gbm.binning.binned import combine_by_factor, rebin_by_time, combine_into_one
py_version = sys.version_info[0]

counts = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
exposure = np.array([1.024, 1.01, 1.00, 0.99, 1.02, 1.024, 0.80, 1.01])
edges = np.linspace(0.0, 8.192, 9)

times = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 1.0, 1.01, 1.02, 1.03, 1.04])
    
class TestRebinningBinned(unittest.TestCase):
    
    def test_combine_by_factor(self): 
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        new_counts, new_exposure, new_edges = combine_by_factor(counts, exposure, edges, 4)
        assertArray(new_counts, np.array([10.0, 26.0]))
        assertArray(new_exposure, np.array([4.024, 3.854]))
        assertArray(new_edges, np.array([0.0, 4.096, 8.192]))
        
        self.assertRaises(AssertionError, combine_by_factor, counts, exposure, edges, 0)

    def test_rebin_by_time(self): 
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        new_counts, new_exposure, new_edges = rebin_by_time(counts, exposure, edges, 4.096)
        assertArray(new_counts, np.array([10.0, 26.0]))
        assertArray(new_exposure, np.array([4.024, 3.854]))
        assertArray(new_edges, np.array([0.0, 4.096, 8.192]))
        
        self.assertRaises(AssertionError, rebin_by_time, counts, exposure, edges, -1.0)

    def test_combine_into_one(self): 
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        new_counts, new_exposure, new_edges = combine_into_one(counts, exposure, edges)
        assertArray(new_counts, np.array([36.0]))
        assertArray(new_exposure, np.array([7.878]))
        assertArray(new_edges, np.array([0.0, 8.192]))

class TestBinningUnbinned(unittest.TestCase):
    
    def test_bin_by_time(self):
        from gbm.binning.unbinned import bin_by_time
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        
        new_counts, new_edges = bin_by_time(times, 1.0)
        assertArray(new_counts, np.array([6, 5]))
        assertArray(new_edges, np.array([0.0, 1.0, 2.0]))
        
        new_counts, new_edges = bin_by_time(times, 0.5, time_ref=1.0)
        assertArray(new_counts, np.array([6, 0, 5]))
        assertArray(new_edges, np.array([0.0, 0.5, 1.0, 1.5]))
        
        self.assertRaises(AssertionError, bin_by_time, times, -1.0)
        
    
    def test_combine_into_one(self):
        from gbm.binning.unbinned import combine_into_one
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        
        new_counts, new_edges = combine_into_one(times, 0.5, 2.0)
        assertArray(new_counts, np.array([5]))
        assertArray(new_edges, np.array([0.5, 2.0]))
        
        self.assertRaises(AssertionError, combine_into_one, times, 2.0, 0.5)
        self.assertRaises(ValueError, combine_into_one, times, 10.0, 20.0)
        
    
    def test_time_to_spill(self):
        from gbm.binning.unbinned import time_to_spill
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        
        new_counts, new_edges = time_to_spill(times, 5)
        assertArray(new_counts, np.array([5, 6]))
        assertArray(new_edges, np.array([0.0, 0.05, 1.04]))
        
        self.assertRaises(AssertionError, time_to_spill, times, -1)
        
        
