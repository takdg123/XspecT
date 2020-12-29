# XspecT

This is an analysis tool for Fermi GBM and LAT data. This tool utilizes 

- PyXspec (2.0.3, https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/index.html)
- gspec (0.9.1, https://fermi.gsfc.nasa.gov/ssc/data/analysis/rmfit/)
- Fermitools (2.0.0, https://github.com/fermi-lat/Fermitools-conda/wiki). 

##### Note that you need to update Xspec functions (/RefData/XSFunctions) beforehand.
- The original step function (xsstep.f) is changed to a cutoff powerlaw function. The functional form is same to the one defined in cutoffPowerLaw.cxx, but the approximation method is different (This uses two-point approximation, whereas cutoffPowerLaw.cxx uses the gamma function).
- The power law (powerLaw.cxx) and cutoff power-law functions (cutoffPowerLaw.cxx) are normalized at 100 keV, instead of 1 keV (original).
- A parameter in the Band function (xsgrbm.f) is changed. Now, a folded energy (E_f) is replaced with the peak energy (E_p), E_p = E_f/(alpha+2).

##### Note that you need to update gspec modules (/RefData/gbm) for generating pha and bak files for LLE.
- LLE data can be obtained from https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/triggers, or using XspecT.Preparation.download.DownloadData with lle=True.
- The energy band for LLE is from 31 MeV to 100 MeV.
