// An emission line from an accretion disk around a BH including GR effects
// Laor ApJ 376, 90
// Parameters
//   0     line energy (keV)
//   1     alpha, power law index for emissivity (R^-alpha)
//   2     inner radius in units of GM/c^2
//   3     outer radius in units of GM/c^2
//   4     inclination (degrees)

#include <xsTypes.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <CCfits/CCfits>
#include <sstream>


void readLaorFile(RealArray& laorRadiusArray, RealArray& laorMidRadiusArray,
		  RealArray& laorAngleArray, RealArray& laorEnergyArray,
		  std::vector<RealArray>& laorTransferFnArray);


void laor(const RealArray& energyArray, const RealArray& params,
	  int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray, 
	  const string& initString)
{
  using namespace Numerics;
  
  static bool first(true);
  static RealArray laorRadiusArray, laorMidRadiusArray, laorAngleArray;
  static RealArray laorEnergyArray;
  static std::vector<RealArray> laorTransferFnArray;
  
  const size_t nE = energyArray.size();
  const size_t nF = nE - 1;

  fluxArray.resize(nF);
  fluxErrArray.resize(0);

  // if this is the first time through read in the model file

  if ( first ) {
    readLaorFile(laorRadiusArray, laorMidRadiusArray, laorAngleArray,
		 laorEnergyArray, laorTransferFnArray);
    first = false;
  }

  // Set up energy array for the line. Note the use of the high energy bin
  // and the low energy bin to suppress error messages from the interpolation
  // program.

  const size_t nlaorE = laorEnergyArray.size();
  RealArray lineEnergyArray(nlaorE+3);

  lineEnergyArray[0] = 0.0;
  lineEnergyArray[1] = (3*laorEnergyArray[0]-laorEnergyArray[1]) * params[0] / 2.;
  for (size_t i=2; i<=nlaorE; i++) {
    lineEnergyArray[i] = (laorEnergyArray[i-2]+laorEnergyArray[i-1]) * params[0] / 2.;
  }
  lineEnergyArray[nlaorE+1] = (3*laorEnergyArray[nlaorE-1]-laorEnergyArray[nlaorE-2]) * params[0] / 2.;
  lineEnergyArray[nlaorE+2] = 1.0e6;

  // Initialize the flux array
  
  RealArray lineFluxArray(nlaorE+2);
  lineFluxArray = 0.0;

  // selecting the integer values of the radial bins which are at the inner and
  // outer disk boundary (according to inner and outer radii)

  const size_t nlaorRad = laorRadiusArray.size();
  size_t iro = 0;
  size_t iri = nlaorRad-1;
  for (size_t i=0; i<nlaorRad-1; i++) {
    if ( params[3] <= laorRadiusArray[i] && params[3] > laorRadiusArray[i+1] ) iro = i;
    if ( params[2] <= laorRadiusArray[i] && params[2] > laorRadiusArray[i+1] ) iri = i;
  }

  // Set the cos(inclination) variable.
  
  Real inclin = cos(params[4]*3.14159/180.);

  //  find the tabulated inclinations above and below that requested
  //  and the associated weights

  const size_t nlaorAngles = laorAngleArray.size();

  size_t inclow = 0;
  while ( inclin > laorAngleArray[inclow] && inclow < nlaorAngles-1 ) inclow++;
  if ( inclow > 0 ) inclow--;
  Real angleBinSize = laorAngleArray[inclow+1] - laorAngleArray[inclow];
  Real wlow = (laorAngleArray[inclow+1]-inclin)/angleBinSize;
  Real whgh = (inclin-laorAngleArray[inclow])/angleBinSize;

  // Summing the profiles of individual rings according to a given
  // weight function

  for (size_t j=iro; j<=iri; j++) {

    Real r = laorRadiusArray[j];
    Real area = 0.0;
    if ( j == iro ) {
      r = laorMidRadiusArray[iro];
      area = params[3]*params[3] - r*r;
    } else if ( j == iri ) {
      r = laorMidRadiusArray[iri-1];
      area = r*r - params[2]*params[2];
    } else {
      area = laorMidRadiusArray[j-1]*laorMidRadiusArray[j-1] -
	     laorMidRadiusArray[j]*laorMidRadiusArray[j];
    }

    // ems-the radial weight function

    Real ems = area * pow(r, -params[1]);

    // sum in the required transfer functions

    size_t ioff = j*nlaorAngles + inclow;

    for (size_t ie=0; ie<nlaorE; ie++) {
      lineFluxArray[ie+1] += wlow*laorTransferFnArray[ioff][ie]*ems;
    }

    ioff ++;
    for (size_t ie=0; ie<nlaorE; ie++) {
      lineFluxArray[ie+1] += whgh*laorTransferFnArray[ioff][ie]*ems;
    }

    // end loop over radii
  }

  // Perform 3-pt smoothing to reduce the small wiggles which result
  // from the finite number of radial bins

  RealArray smo(3);

  smo[0] = lineFluxArray[1];
  smo[1] = lineFluxArray[2];
  size_t ismo = 2;
  for (size_t i=1; i<nlaorE-2; i++) {
    smo[ismo] = lineFluxArray[i+2];
    ismo++;
    if ( ismo == 3 ) ismo = 0;
    lineFluxArray[i+1] = 0.;
    for (size_t j=0; j<3; j++) lineFluxArray[i+1] += smo[j];
    lineFluxArray[i+1] /= 3.0;
  }

  // Rebin onto passed energies

  size_t inputBin;
  size_t outputBin;
  IntegerVector startBin(fluxArray.size());
  IntegerVector endBin(fluxArray.size());
  RealArray startWeight(fluxArray.size());
  RealArray endWeight(fluxArray.size());

  const Real FUZZY = 1.0e-6;

  Rebin::findFirstBins(lineEnergyArray, energyArray, FUZZY, inputBin, outputBin);
  Rebin::initializeBins(lineEnergyArray, energyArray, FUZZY, inputBin, outputBin,
				  startBin, endBin, startWeight, endWeight);

  Rebin::rebin(lineFluxArray, startBin, endBin, startWeight, endWeight, fluxArray);

  // normalise values to total

  Real spm = fluxArray.sum();
  if ( spm > 0.0 && spm != 1.0 ) fluxArray /= spm;

}
  
void readLaorFile(RealArray& laorRadiusArray, RealArray& laorMidRadiusArray,
		  RealArray& laorAngleArray, RealArray& laorEnergyArray,
		  std::vector<RealArray>& laorTransferFnArray)
{
  using namespace CCfits;

  string modelFileName(FunctionUtility::modelDataPath());
  modelFileName += "ari.mod";

  std::unique_ptr<FITS> pInfile((FITS*)0);

  try {
    pInfile.reset(new FITS(modelFileName));
  } catch(...) {
    string msg = "laor: failed to open " + modelFileName + "\n";
    throw YellowAlert(msg);
  }

  try {
    ExtHDU& radiusExt = pInfile->extension("RADII");
    long nRows = (long)radiusExt.rows();
    radiusExt.column("RADIUS").read(laorRadiusArray, (long)1, nRows);
  } catch(...) {
    string msg = "laor: failed to read RADIUS column from RADII extension.\n";
    throw YellowAlert(msg);
  }
      

  laorMidRadiusArray.resize(laorRadiusArray.size()-1);
  for (size_t i=0; i<laorMidRadiusArray.size(); i++) {
    laorMidRadiusArray[i] = (laorRadiusArray[i+1] + laorRadiusArray[i])/2.0;
  }

  try {
    ExtHDU& angleExt = pInfile->extension("ANGLES");
    long nRows = (long)angleExt.rows();
    angleExt.column("ANGLE").read(laorAngleArray, (long)1, nRows);
  } catch(...) {
    string msg = "laor: failed to read ANGLE column from ANGLES extension.\n";
    throw YellowAlert(msg);
  }
      
  try {
    ExtHDU& energyExt = pInfile->extension("ENERGIES");
    long nRows = (long)energyExt.rows();
    energyExt.column("ENERGY").read(laorEnergyArray, (long)1, nRows);
  } catch(...) {
    string msg = "laor: failed to read ENERGY column from ENERGIES extension.\n";
    throw YellowAlert(msg);
  }
      

  try {
    ExtHDU& transferExt = pInfile->extension("TRANSFER FUNCTIONS");
    long nRows = (long)transferExt.rows();
    laorTransferFnArray.resize(nRows);
    transferExt.column("TRANSFER").readArrays(laorTransferFnArray, (long)1, nRows);
  } catch(...) {
    string msg = "laor: failed to read TRANSFER column from TRANSFER FUNCTIONS extension.\n";
    throw YellowAlert(msg);
  }
      
}




