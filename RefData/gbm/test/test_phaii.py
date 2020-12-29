import unittest, os, shutil, sys
from gbm.data.phaii import RateHisto, Ctime, Cspec, TTE
from gbm.binning.unbinned import bin_by_time
from gbm.selection import Selection
import numpy as np

py_version = sys.version_info[0]

class TestPHAII(unittest.TestCase):
    
    def test_RateHisto_exposure(self):        
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        
        counts = np.array([5.0, 1.0, 0.0])
        bounds = (np.array([1.0, 10.0, 12.0]), np.array([10.0, 12.0, 15.0]))
        exp_float = 0.5
        exp_array = np.array([8.5, 1.4, 0.0])
        
        # no exposure, use binwidth
        rates1 = RateHisto(counts, bounds)
        assertArray(rates1.rate, np.array([5.0/9.0, 1.0/2.0, 0.0]))
        assertArray(rates1.uncertainty, np.array([np.sqrt(5.0)/9.0, 1.0/2.0, 0.0]))
        
        # scalar exposure
        rates2 = RateHisto(counts, bounds, exposure=exp_float)
        assertArray(rates2.rate, np.array([10.0, 2.0, 0.0]))
        assertArray(rates2.uncertainty, np.array([np.sqrt(5.0)/0.5, 2.0, 0.0]))
        
        # exposure array
        rates3 = RateHisto(counts, bounds, exposure=exp_array)
        assertArray(rates3.rate, np.array([5.0/8.5, 1.0/1.4, 0.0]))
        assertArray(rates3.uncertainty, np.array([np.sqrt(5.0)/8.5, 1.0/1.4, 0.0]))
    
    def test_RateHisto_sum(self):
        if py_version == 2:
            assertArray = self.assertItemsEqual
        elif py_version == 3:
            assertArray = self.assertCountEqual
        
        counts = np.array([5.0, 1.0, 0.0])
        bounds = (np.array([1.0, 10.0, 12.0]), \
                  np.array([10.0, 12.0, 15.0]))
        rates = RateHisto(counts, bounds)
        
        # adding rates
        rates_list = [rates, rates]
        rates_sum = RateHisto.sum_time(rates_list)
        assertArray(rates_sum.rate, np.array([10.0/9.0, 1.0, 0.0]))
        assertArray(rates_sum.uncertainty, \
                              np.array([np.sqrt(10.0)/9.0, np.sqrt(2.0)/2.0, 0.0]))
    
    def test_RateHisto_errors(self):
        counts = None
        bounds = (np.array([1.0, 10.0, 12.0]), np.array([10.0, 12.0, 15.0]))

        self.assertRaises(AssertionError, RateHisto, counts, bounds) 
        
        # changed from an error to a warning
        #counts = np.array([10.0, -1.0, 11.0])
        #self.assertRaises(ValueError, RateHisto, counts, bounds)
        
        exposure = '0.5'
        counts = np.array([5.0, 1.0, 0.0])
        self.assertRaises(TypeError, RateHisto, counts, bounds, exposure=exposure)

        exposure = np.array([0.5, 1.0])
        self.assertRaises(AssertionError, RateHisto, counts, bounds, exposure=exposure)
                
        bounds2 = (np.array([1.0, 11.0, 13.0]), np.array([11.0, 13.0, 15.0]))
        rates1 = RateHisto(counts, bounds)
        rates2 = RateHisto(counts, bounds2)
        self.assertRaises(AssertionError, RateHisto.sum_time, [rates1, rates2])
        

    
    def test_ctime_phaii(self):
        ctime = Ctime('glg_ctime_nb_bn120415958_v00.pha')
        time_range = ctime.time_range
        energy_range = ctime.energy_range
        numchans = ctime.num_chans
        self.assertAlmostEqual(time_range[0],  -899.342442035675, places=6)
        self.assertAlmostEqual(time_range[1], 1000.8578699827194, places=6)
        self.assertAlmostEqual(energy_range[0], 4.323754, places=6)
        self.assertAlmostEqual(energy_range[1], 2000.0, places=6)
        self.assertEqual(numchans, 8)
        self.assertAlmostEqual(ctime.trigtime, 356223561.133346, places=6)
        
        lightcurve = ctime.view(integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 14429)
        count_spectrum = ctime.view(integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 8)
        
        time_slice = (time_range[0]+100.0, time_range[1]-100.0)
        lightcurve = ctime.view(time_range=time_slice, integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 13649)
        count_spectrum = ctime.view(time_range=time_slice, integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 8)
        
        channel_slice = (3, 4)
        lightcurve = ctime.view(channel_range=channel_slice, integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 14429)
        count_spectrum = ctime.view(channel_range=channel_slice, integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 2)
        
        energy_slice = (50.0, 300.0)
        lightcurve = ctime.view(energy_range=energy_slice, integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 14429)
        count_spectrum = ctime.view(energy_range=energy_slice, integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 3)
        
        lightcurve = ctime.view(time_range=time_slice, energy_range=energy_slice, \
                                integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 13649)
        count_spectrum = ctime.view(time_range=time_slice, energy_range=energy_slice, \
                                    integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 3)
        
    def test_cspec_phaii(self):
        cspec = Cspec('glg_cspec_b0_bn120415958_v00.pha')
        time_range = cspec.time_range
        energy_range = cspec.energy_range
        numchans = cspec.num_chans
        self.assertAlmostEqual(time_range[0], -4001.872472047806, places=6)
        self.assertAlmostEqual(time_range[1], 3805.9237439632416, places=6)
        self.assertAlmostEqual(energy_range[0], 114.5423, places=4)
        self.assertAlmostEqual(energy_range[1], 50000.0, places=6)
        self.assertEqual(numchans, 128)
        self.assertAlmostEqual(cspec.trigtime, 356223561.133346, places=6)
        
        lightcurve = cspec.view(integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 404)
        count_spectrum = cspec.view(integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 128)
        
        time_slice = (time_range[0]+100.0, time_range[1]-100.0)
        lightcurve = cspec.view(time_range=time_slice, integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 380)
        count_spectrum = cspec.view(time_range=time_slice, integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 128)
        
        channel_slice = (3, 4)
        lightcurve = cspec.view(channel_range=channel_slice, integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 404)
        count_spectrum = cspec.view(channel_range=channel_slice, integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 2)
        
        energy_slice = (300.0, 1000.0)
        lightcurve = cspec.view(energy_range=energy_slice, integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 404)
        count_spectrum = cspec.view(energy_range=energy_slice, integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 19)
        
        lightcurve = cspec.view(time_range=time_slice, energy_range=energy_slice, \
                                integration_axis='energy')
        self.assertEqual(lightcurve[0].num_bins, 380)
        count_spectrum = cspec.view(time_range=time_slice, energy_range=energy_slice, \
                                    integration_axis='time')
        self.assertEqual(count_spectrum[0].num_bins, 19)

    
    def test_tte_phaii(self):
        tte = TTE('glg_tte_n9_bn090131090_v00.fit')
        time_range = tte.time_range
        energy_range = tte.energy_range
        numchans = tte.num_chans
        self.assertAlmostEqual(time_range[0], -25.49109798669815, places=6)
        self.assertAlmostEqual(time_range[1], 300.73552399873734, places=6)
        self.assertAlmostEqual(energy_range[0], 4.389729, places=6)
        self.assertAlmostEqual(energy_range[1], 2000.0, places=6)
        self.assertEqual(numchans, 128)
        self.assertAlmostEqual(tte.trigtime, 255060563.149072, places=6)
        
        select = Selection(tte)
        select.set_time_binning(bin_by_time, 0.256, time_ref=0.0)
        lightcurve = select.lightcurve
        self.assertEqual(lightcurve[0].num_bins, 1275)
        
        count_spectrum = select.spectrum
        self.assertEqual(count_spectrum[0].num_bins, 128)
        
        time_slice = (time_range[0]+10.0, time_range[1]-10.0)
        select2 = select.slice(time_selections=time_slice)
        lightcurve = select2.lightcurve
        self.assertEqual(lightcurve[0].num_bins, 1197)
        count_spectrum = select2.spectrum
        self.assertEqual(count_spectrum[0].num_bins, 128)
                
        energy_slice = (50.0, 300.0)
        select2 = select.slice(energy_selections=energy_slice)
        lightcurve = select2.lightcurve
        self.assertEqual(lightcurve[0].num_bins, 1275)
        count_spectrum = select2.spectrum
        self.assertEqual(count_spectrum[0].num_bins, 54)
        
        select2 = select.slice(time_selections=time_slice, energy_selections=energy_slice)
        lightcurve = select2.lightcurve
        self.assertEqual(lightcurve[0].num_bins, 1197)
        count_spectrum = select2.spectrum
        self.assertEqual(count_spectrum[0].num_bins, 54)
        
        select.set_time_binning(bin_by_time, 0.016, time_ref=0.0)
        lightcurve = select.lightcurve
        self.assertEqual(lightcurve[0].num_bins, 20381)
        select.set_time_binning(bin_by_time, 2.048, time_ref=0.0)
        lightcurve = select.lightcurve
        self.assertEqual(lightcurve[0].num_bins, 160)
        
    
    def test_phaii_errors(self):
        ctime = Ctime('glg_ctime_nb_bn120415958_v00.pha')
        time_range = ctime.time_range
        energy_range = ctime.energy_range
        numchans = ctime.num_chans
        
        self.assertRaises(ValueError, ctime.view, integration_axis='peekaboo')
        self.assertRaises(AssertionError, ctime.view, time_range=(time_range[1], time_range[0]))
        self.assertRaises(AssertionError, ctime.view, energy_range=(-10.0, 100.0))
        self.assertRaises(AssertionError, ctime.view, channel_range=(-10, -20))
               

if __name__ == '__main__':
    unittest.main()
