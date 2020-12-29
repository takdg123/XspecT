# Copyright (C) 2016 William Cleveland, Universities Space Research Association (USRA)
# for NASA's Fermi Gamma-ray Burst Monitor (GBM).

import unittest
from gbm.detectors import *
from copy import copy


class TestDetectors(unittest.TestCase):

    expected_nai = ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb']
    expected_bgo = ['b0', 'b1']
    expected_all = expected_nai + expected_bgo

    def test_detector_list(self):
        detectors = copy(self.expected_all)
        for d in Detector:
            if d.short_name in detectors:
                detectors.remove(d.short_name)
            else:
                self.fail("Detector was not expected")
        if detectors:
            self.fail("All of the detectors weren't removed")

    def test_is_nai(self):
        for d in self.expected_nai:
            self.assertTrue(Detector.is_nai(Detector.from_str(d)))
        for d in self.expected_bgo:
            self.assertFalse(Detector.is_nai(Detector.from_str(d)))

    def test_is_bgo(self):
        for d in self.expected_bgo:
            self.assertTrue(Detector.is_bgo(Detector.from_str(d)))
        for d in self.expected_nai:
            self.assertFalse(Detector.is_bgo(Detector.from_str(d)))

    def test_nai(self):
        l = Detector.nai()
        self.assertEquals(len(l), len(self.expected_nai))
        for d in self.expected_nai:
            self.assertIsNotNone(Detector.from_str(d))

    def test_bgo(self):
        l = Detector.bgo()
        self.assertEquals(len(l), len(self.expected_bgo))
        for d in self.expected_bgo:
            self.assertIsNotNone(Detector.from_str(d))

    def test_from_number(self):
        num = 0
        for d in self.expected_all:
            self.assertEquals(d, Detector.from_num(num).short_name)
            num += 1

    def test_meta_data(self):
        expected_data = [
            ('N0', 'NAI_00', 0, 45.89, 20.58),
            ('N1', 'NAI_01', 1, 45.11, 45.31),
            ('N2', 'NAI_02', 2, 58.44, 90.21),
            ('N3', 'NAI_03', 3, 314.87, 45.24),
            ('N4', 'NAI_04', 4, 303.15, 90.27),
            ('N5', 'NAI_05', 5, 3.35, 89.79),
            ('N6', 'NAI_06', 6, 224.93, 20.43),
            ('N7', 'NAI_07', 7, 224.62, 46.18),
            ('N8', 'NAI_08', 8, 236.61, 89.97),
            ('N9', 'NAI_09', 9, 135.19, 45.55),
            ('NA', 'NAI_10', 10, 123.73, 90.42),
            ('NB', 'NAI_11', 11, 183.74, 90.32),
            ('B0', 'BGO_00', 12, 0.00, 90.00),
            ('B1', 'BGO_01', 13, 180.00, 90.00),
        ]
        for d in expected_data:
            det = Detector.from_num(d[2])
            self.assertEquals(det.name,  d[0])
            self.assertEquals(det.long_name, d[1])
            self.assertEquals(det.azimuth, d[3])
            self.assertEquals(det.zenith, d[4])
            self.assertEquals(det.pointing, (d[3], d[4]))
            self.assertEquals(det.__repr__(), "Detector(\"{}\", \"{}\", {})".format(d[0], d[1], d[2]))

    def test_invalid_str(self):
        self.assertIsNone(Detector.from_str("FAKE"))

    def test_invalid_num(self):
        self.assertIsNone(Detector.from_num(20))

