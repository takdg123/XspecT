# XspecT

This is an analysis tool for Fermi GBM and LAT data. This tool utilizes 

- Xspec (12.11.1, https://heasarc.gsfc.nasa.gov/xanadu/xspec/)
- gspec (0.9.1, https://fermi.gsfc.nasa.gov/ssc/data/analysis/rmfit/)
- Fermitools (2.0.0, https://github.com/fermi-lat/Fermitools-conda/wiki). 

Note that you need to update Xspec functions (/RefData/XSFunctions) before using this tool.
- The original step function (xsstep.f) is changed to a cutoff powerlaw function (the functional form is same to the one defined in cutoffPowerLaw.cxx, but the approximation method is different).
- The power law (powerLaw.cxx) and cutoff power-law functions (cutoffPowerLaw.cxx) are normalized at 100 keV, instead of 1 keV (original)
- A parameter in the Band function (xsgrbm.f) is changed. Now, a folded energy (E_f) is replaced with the peak energy (E_p), E_p = E_f/(\alpha+2))
