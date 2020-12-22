// XSPEC model subroutine for redshifted "partial covering 
// absorption". Assumes that part of the emitter is covered by 
// the given absorption and the rest of the emitter is unobscured. 
//---
// number of model parameters: 4
//      0      ANH      Hydrogen column density (in units of 10**22
//                      atoms per square centimeter
//      1      xi
//      2      FRACT    Covering fraction (0 implies no absorption,
//                      1 implies the emitter is all absorbed with
//                      the indicated column ANH.
//      3      REDSHIFT

#include <xsTypes.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

// prototype from tableInterpolate.cxx

void tableInterpolate(const RealArray& energyArray, const RealArray& params, 
		      string fileName, int spectrumNumber, RealArray& fluxArray, 
		      RealArray& fluxErrArray, const string& initString, 
		      const string& tableType, const bool readFull);


void zxipcf(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& trans, RealArray& transErr, 
	    const string& initString)
{
  RealArray eparams(3);
  eparams[0] = params[0]*1.e22;   //  nh in units of 1e22 for param(0)      
  eparams[1] = params[1];         //  logxi
  eparams[2] = params[3];         //  z

  // find the path to the mtable file required

  string pname = "ZXIPCF_DIR";
  string DirName(FunctionUtility::getModelString(pname));
  if ( DirName.length() == 0 || DirName == FunctionUtility::NOT_A_KEY() )
    DirName = FunctionUtility::modelDataPath();

  string fileName = DirName + "zxipcf_mtable.fits";

// interpolate on the mtable

  tableInterpolate(energyArray, eparams, fileName, spectrumNumber, trans, 
		   transErr, initString, "mul", true);
  if ( trans.size() == 0 ) {
    FunctionUtility::xsWrite("Failed to read "+fileName+" use xset ZXIPCF_DIR to set directory containing file", 5);
    return;
  }

  // now modify for the partial covering

  trans = (1.0-params[2]) + params[2]*trans;
  if (transErr.size() > 0 ) transErr = (1.0-params[2]) + params[2]*transErr;

  return;
}

