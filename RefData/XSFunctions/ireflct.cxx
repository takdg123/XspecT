#include "MZCompRefl.h"

// definition of wrap-up routine

void doreflect(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString, Real Emax);


// Main model routine

void ireflct(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
  //  Driver for angle-dependent reflection
  //
  //  See Magdziarz & Zdziarski, 1995, MNRAS.
  //  See Zdziarski et al., 1995, ApJL 438, L63 for description of
  //  calculation of ionization (based on Done et al. 1992, ApJ 395, 275).
  //  The abundances are defined by the command abund
  //
  //  The output spectrum is the sum of the input (the original contents
  //  of Photar) and the reflection component.
  //  The reflection component alone can be obtained
  //  for scale (see below) = rel_refl < 0. Then the actual
  //  reflection normalization is |scale|. Note that you need to
  //  change the limits of rel_refl. The range of rel_refl in that case
  //  should exclude zero (as then the direct component appears).
  //
  //  number of model parameters:7
  //  1: scale, scaling factor for reflection; if <0, no direct component
  //     (scale=1 for isotropic source above disk)
  //  2: redshift, z
  //  3: abundance of elements heavier than He relative to
  //     the solar abundances
  //  4: iron abundance relative to the solar abundances
  //  5: cosine of inclination angle
  //  6: disk temperature in K
  //  7: disk ionization parameter = L/nR^2
  //
  // also uses the following values which can be set using xset
  //   IREFLECT_MAX_E      maximum energy for which to calculate the output spectrum
  //   IREFLECT_PRECISION  fractional precision for Greens' fn. adaptive integration
  //  Based on reflct.cxx

  using namespace XSutility;

  // find whether a model name has been give in the initString.

  string modelName = "IREFLECT";
  if (initString.find_first_not_of(" ") != string::npos) modelName = initString;

  // check whether the user has defined a maximum energy that they want calculated

  Real Emax(-1.0);
  string pname = modelName + "_MAX_E";
  string pvalue = FunctionUtility::getModelString(pname);
  if ( pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY() ) {
    std::istringstream iss(pvalue);
    if ( !(iss >> Emax) || !iss.eof() ) {
      std::ostringstream err;
      err << "Invalid " << modelName << "_MAX_E value: " << pvalue;
      xs_write(const_cast<char*>(err.str().c_str()), 10);
    }
  }
  

  doreflect(energyArray, params, spectrumNumber, flux, fluxErr, modelName, Emax);

}


void doreflect(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& modelName, Real Emax)
{
  
  using namespace XSutility;

  Real Xmax = Emax/511.0;

  static bool isFirst = true;
  if (isFirst) {
    if ( modelName == "IREFLECT" || modelName == "REFLECT" ) {
      string msg("Compton reflection. See help for details.");
      msg += "\nIf you use results of this model in a paper,";
      msg += "\nplease refer to Magdziarz & Zdziarski 1995 MNRAS, 273, 837";
      xs_write(const_cast<char*>(msg.c_str()), 5);
    }
    isFirst = false;
  }
  
  static int nesave = -1;
  const int Ne = static_cast<int>(energyArray.size()) - 1;
  static RealArray Sptot;
  static RealArray Spinc;
  static RealArray X;

  if (Ne != nesave) {
    Sptot.resize(Ne);
    Spinc.resize(Ne);
    X.resize(Ne);
    nesave = Ne;
  }

  // parameters

  Real Scale = params[0];
  Real zshift = 1.0 + params[1];
  Real Abund = log10(params[2]);
  Real FeAbund  = log10(params[3]);
  Real cosIncl = params[4];
  if ( cosIncl < 0.05 ) cosIncl = 0.05;
  if ( cosIncl > 0.95 ) cosIncl = 0.95;
  Real Temp = params[5];
  Real Xi   = params[6];

  //  x is in the source frame and in units of m_e c^2.

  Real factor = zshift/2.0/511.0;
  for (int i=0; i<Ne; ++i)
    X[i] = (energyArray[i+1]+energyArray[i]) * factor;

  // Generate the spinc array from the input photar array.

  for (int i=0; i<Ne; ++i) {
    if (energyArray[i+1] != energyArray[i]) {
      Real avgE = (energyArray[i+1]+energyArray[i])/2.;
      Spinc[i] = flux[i]*avgE*avgE/(energyArray[i+1]-energyArray[i]);
    }
  }

  //     x is the energy array (units m_e c^2)
  //     spinc is the input spectrum array (E F_E)
  //     spref is the reflected spectrum array (E F_E)
  //     sptot is the total spectrum array (E F_E), = spinc if no reflection
  //     all dimensions = Ne
  

  // calculate the total flux including direct and Compton reflection (non-relativistic
  // and relativistic).

  calcCompReflTotalFlux(modelName, Scale, cosIncl, Abund, FeAbund, Xi, Temp, Xmax,
			X, Spinc, Sptot);

  // convert back to xspec internal units
  
  for (int i=0; i<Ne; ++i) {
    Real avgE = (energyArray[i+1]+energyArray[i])/2.;
    flux[i] = Sptot[i]*(energyArray[i+1]-energyArray[i])/(avgE*avgE);
  }

}
