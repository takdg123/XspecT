# detectors.py: GBM Detector module
# Copyright (C) 2016 William Cleveland, Universities Space Research Association (USRA)
# for NASA's Fermi Gamma-ray Burst Monitor (GBM).
# April 20, 2016
from enum import Enum


class Detector(Enum):
    N0 = ('NAI_00', 0, 45.89, 20.58)
    N1 = ('NAI_01', 1, 45.11, 45.31)
    N2 = ('NAI_02', 2, 58.44, 90.21)
    N3 = ('NAI_03', 3, 314.87, 45.24)
    N4 = ('NAI_04', 4, 303.15, 90.27)
    N5 = ('NAI_05', 5, 3.35, 89.79)
    N6 = ('NAI_06', 6, 224.93, 20.43)
    N7 = ('NAI_07', 7, 224.62, 46.18)
    N8 = ('NAI_08', 8, 236.61, 89.97)
    N9 = ('NAI_09', 9, 135.19, 45.55)
    NA = ('NAI_10', 10, 123.73, 90.42)
    NB = ('NAI_11', 11, 183.74, 90.32)
    B0 = ('BGO_00', 12, 0.00, 90.00)
    B1 = ('BGO_01', 13, 180.00, 90.00)
    LL = ('LLE', 14, 0, 0)

    def __init__(self, long_name, number, azimuth, zenith):
        self.long_name = long_name
        self.number = number
        self.azimuth = azimuth
        self.zenith = zenith

    def __repr__(self):
        return "Detector(\"{}\", \"{}\", {})".format(self.name, self.long_name, self.number)

    def __str__(self):
        return str.lower(self.name)

    @property
    def short_name(self):
        return self.__str__()

    @property
    def pointing(self):
        return self.azimuth, self.zenith

    def is_nai(self):
        return self.name[0] == 'N'

    def is_bgo(self):
        return self.name[0] == 'B'

    @classmethod
    def from_str(cls, value):
        if value.upper() in cls.__members__:
            return cls[value.upper()]
        # TODO: Reconsider returning None
        return None

    @classmethod
    def from_num(cls, num):
        for d in cls:
            if d.number == num:
                return d
        # TODO: Reconsider returning None
        return None

    @classmethod
    def nai(cls):
        return [x for x in cls if x.is_nai()]

    @classmethod
    def bgo(cls):
        return [x for x in Detector if x.is_bgo()]
