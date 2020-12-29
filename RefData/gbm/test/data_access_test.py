import unittest, os, shutil
import numpy as np

from gbm.data.phaii import Ctime, Cspec, TTE
from gbm.data.drm import RSP


class TestDataAccess(unittest.TestCase):
    
    def test_read_ctime(self):
        filename = 'glg_ctime_nb_bn120415958_v00.pha'
        ctime = Ctime(filename)
        self.assertEqual(ctime.filename, filename)
        self.assertEqual(len(ctime.header), 4)
        self.assertEqual(len(ctime.ebounds), 8)
        self.assertEqual(len(ctime.spectrum), 14429)
        self.assertEqual(len(ctime.gti), 1)
      
    def test_read_cspec(self):
        filename = 'glg_cspec_b0_bn120415958_v00.pha'
        cspec = Cspec(filename)
        self.assertEqual(cspec.filename, filename)
        self.assertEqual(len(cspec.header), 4)
        self.assertEqual(len(cspec.ebounds), 128)
        self.assertEqual(len(cspec.spectrum), 1993)
        self.assertEqual(len(cspec.gti), 2)

    def test_read_tte(self):
        filename = 'glg_tte_n9_bn090131090_v00.fit'
        tte = TTE(filename)
        self.assertEqual(tte.filename, filename)
        self.assertEqual(len(tte.header), 4)
        self.assertEqual(len(tte.ebounds), 128)
        self.assertEqual(len(tte.events), 422405)
        self.assertEqual(len(tte.gti), 1)

    def test_read_rsp2(self):
        filename = 'glg_cspec_n4_bn120415958_v00.rsp2'
        rsp = RSP(filename)
        self.assertEqual(rsp.filename, filename)
        self.assertEqual(len(rsp.header), 14)
        self.assertEqual(len(rsp.ebounds), 128)
        self.assertEqual(len(rsp.pbounds[0]), 140)
        self.assertEqual(rsp.num_drms, 12)
        self.assertEqual(len(rsp.drm_list), 12)
        self.assertEqual(len(rsp.times), 13)
        self.assertEqual(rsp.data_type, 'CSPEC')
        
        # Test to see that weightResponse picks out the correct matrix:
        times = [0., 1.]
        self.assertEqual(rsp.checkResponse(times), 'glg_cspec_n4_bn120415958_v00.rsp2{4}')
        matrix3 = rsp.drm_list[3]
        matrix2 = rsp.weightResponse(tint = times)
        self.assertEqual(matrix2.shape, (140, 128))
        self.assertTrue(np.array_equal(matrix3, matrix2))
        
        # Exercise the weighting function:
        hist = np.zeros(10, dtype=float) + .1
        lo = np.arange(0., 50., 5.)
        hi = lo + 5.
        times = np.array([lo, hi])
        matrix = rsp.weightResponse(norm_hist = hist, time_bins = times)
        matrix4 = rsp.drm_list[4]
        weightedRSP = matrix3 * 0.7 + matrix4 * 0.3
        self.assertTrue(np.allclose(matrix, weightedRSP))
        
if __name__ == '__main__':
    unittest.main()

