// Hydrogen model atmosphere model by V. Suleimanov
// with Comptonization taken into account using
// Kompaneets equation.
// Implemented by D. Klochkov 2015
// email: klochkov@astro.uni-tuebingen.de
// modified by kaa 5/12/17 to use a single FITS file for the input table
// and to share code with carbatm.cxx.
//
// parameter[0] - effective (unredshifted) temperature of the 
//                neutron star surface (in MK)
//                T = 2.0..10.0
// parameter[1] - neutron star gravitational mass (in solar mass)
// parameter[2] - neutron star radius (in km)
//
//////////////////////////////////////////////////

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>
#include "XSFunctions/Utilities/FunctionUtility.h"
#include "BinarySearch.h"
#include <CCfits/CCfits>

// this routine can be found in carbatm.cxx
void calcNSatm(const RealArray& energy, const RealArray& parameter, 
	       const string inputFilename, RealArray& flux);


void hatm(const RealArray& energy, const RealArray& parameter, int spectrum, 
	  RealArray& flux, RealArray& fluxError, const string& init)
{

  //..path to the directory with model tables:
  // if variable HATM is set, use that path. If not, use
  // $HEADAS/../spectral/modelData
  string directory = FunctionUtility::getModelString("HATM");
  if(directory.compare("$$NOT$$") == 0) { // HATM variable is not set
    directory = FunctionUtility::modelDataPath();
  }
  //..initialize filename
  const string inputFilename(directory+"spH.fits");

  calcNSatm(energy, parameter, inputFilename, flux);
  fluxError.resize(0);
  return;
}
