#include <MZCompRefl.h>
#include <XSFunctions/Utilities/FunctionUtility.h>

// prototypes for functions in this file

void calcXilRefl(const string modelName, const RealArray& energyArray,
		 RealArray& tableParams, Real index, RealArray& Spref);

// prototypes for functions in rfxconv.cxx

void calcReflFlux(const string modelName, RealArray& tableParams, Real Scale, 
		  Real cosIncl, Real Abund, Real cutoff, Real zshift, 
		  const RealArray& energyArray, RealArray& flux);

// prototype from tableInterpolate.cxx

void tableInterpolate(const RealArray& energyArray, const RealArray& params, 
		      string fileName, int spectrumNumber, RealArray& fluxArray, 
		      RealArray& fluxErrArray, const string& initString,
		      const string& tableType, const bool readFull);



static const Real PI = 4.0 * atan(1.0);

// filename for ionization model

#define XILLVER_FILE "xillver-a-Ec3.fits"


// Main model routine

void xilconv(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
   //  Reflection from an ionized disk - combines the Magdziarz & Zdziarski
   //  Compton reflection code at high energies with tabulated models of
   //  ionized disks at lower energies. 
   //
   //  See Magdziarz & Zdziarski, 1995, MNRAS for Comptonization model.
   //
   //  The output spectrum is the sum of the input (the original contents
   //  of Photar) and the reflection component.
   //  The reflection component alone can be obtained
   //  for scale (see below) = rel_refl < 0. Then the actual
   //  reflection normalization is |scale|. Note that you need to
   //  change the limits of rel_refl. The range of rel_refl in that case
   //  should exclude zero (as then the direct component appears).
   //
   //  number of model parameters:6
   //  1: scale, scaling factor for reflection; if <0, no direct component
   //     (scale=1 for isotropic source above disk)
   //  2: redshift, z
   //  3: iron abundance relative to the solar abundances
   //  4: cosine of inclination angle
   //  5: log10(Xi) for ionization
   //  6: Exponential cut-off energy
   //
   // also uses the following values which can be set using xset
   //   XILCONV_MAX_E      maximum energy for which to calculate the output spectrum
   //   XILCONV_PRECISION  precision used in adaptive Gauss-Kronrod quadrature of Greens'
   //                      function integral
   //   XILCONV_DIR        directory to use for ionization model filenames instead of
   //                      the standard modelData directory
   //  Based on rfxconv.cxx
   
   using namespace XSutility;

   static bool isFirst = true;
   if (isFirst)
   {
      string msg("Compton reflection from ionized medium.");
      msg += "\nIf you use results from this model in a paper";
      msg += "\nplease refer to Chris Done";
      xs_write(const_cast<char*>(msg.c_str()), 5);
      isFirst = false;
   }

   // parameters

   Real Scale = params[0];
   Real zshift = 1.0 + params[1];
   Real Abund = 1.0;
   Real cosIncl = params[3];
   if ( cosIncl < 0.05 ) cosIncl = 0.05;
   if ( cosIncl > 0.95 ) cosIncl = 0.95;
   Real cutoff = params[5];

   // parameters that will be used to get the spectrum from the xillver table model

   RealArray tableParams(6);
   tableParams[1] = params[2];
   tableParams[2] = params[4];
   tableParams[3] = cutoff;
   tableParams[4] = acos(cosIncl) * 180.0 / PI;
   tableParams[5] = 0.0;

   // calculate the Compton reflection (non-relativistic and relativistic) using 
   // the XILLVER table model for the ionization.

   calcReflFlux("XILCONV", tableParams, Scale, cosIncl, Abund, cutoff, zshift,
		energyArray, flux);

}

// *****************************************************************************
// Use the XILLVER ionized disk table model to calculate reflected flux

void calcXilRefl(const string modelName, const RealArray& energyArray,
		 RealArray& tableParams, Real index, RealArray& Spref)
{

  Spref.resize(energyArray.size()-1);
  Spref = 0.0;

  tableParams[0] = index;

  // Table model is only tabulated for slopes of 1.2 to 3.4 so if calculated slope
  // lies outside this range then use extrema and write warning message.

  if ( tableParams[0] < 1.2 ) {
    stringstream msg;
    msg << modelName << ": A 2-10 keV slope of " << tableParams[0] << " is outside the tabulated range - setting to 1.2" << endl;
    FunctionUtility::xsWrite(msg.str(), 5);
    tableParams[0] = 1.2;
  } else if ( tableParams[0] > 3.4 ) {
    stringstream msg;
    msg << modelName << ": A 2-10 keV slope of " << tableParams[0] << " is outside the tabulated range - setting to 3.4" << endl;
    FunctionUtility::xsWrite(msg.str(), 5);
    tableParams[0] = 3.4;
  }

  // Check whether a directory has been set for the input table files.
  // If not use the standard location

  string pname = modelName+"_DIR";
  string DirName(FunctionUtility::getModelString(pname));
  
  if ( DirName.length() == 0 || DirName == FunctionUtility::NOT_A_KEY() )
    DirName =  FunctionUtility::modelDataPath();

  string FileName =  DirName + XILLVER_FILE;

  RealArray fluxErrArray;

  int ifl = 0;

  tableInterpolate(energyArray, tableParams, FileName, ifl, Spref, fluxErrArray, " ", 
		   "add", true);
  if ( Spref.size() == 0 ) {
    stringstream msg;
    msg << modelName << ": Failed to read " << FileName << endl;
    FunctionUtility::xsWrite(msg.str(),5);
  }

}
