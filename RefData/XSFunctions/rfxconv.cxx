#include <MZCompRefl.h>

// prototypes for functions in this file

void calcReflFlux(const string modelName, RealArray& tableParams, Real Scale, 
		  Real cosIncl, Real Abund, Real cutoff, Real zshift, 
		  const RealArray& energyArray, RealArray& flux);

void calcRossRefl(const string modelName, const RealArray& energyArray,
		  RealArray& tableParams, Real index, RealArray& Spref);

Real estimatePowerLawSlope(Real Xlow, Real Xhigh, RealArray& X, RealArray& Y, bool& good);

// prototype for functions in xilconv.cxx

void calcXilRefl(const string modelName, const RealArray& energyArray,
		 RealArray& tableParams, Real index, RealArray& Spref);

// prototype from cutoffPowerLaw.cxx

Real calcCutoffPowerLaw(const RealArray& energyArray, const Real& index, 
			const Real& cutoff, bool isRenorm, RealArray& fluxArray);

// prototype from tableInterpolate.cxx

void tableInterpolate(const RealArray& energyArray, const RealArray& params, 
		      string fileName, int spectrumNumber, RealArray& fluxArray, 
		      RealArray& fluxErrArray, const string& initString, 
		      const string& tableType, const bool readFull);

// handy energies in m_e c^2 units

#define X2       2.0/511.0
#define X10     10.0/511.0
#define X12     12.0/511.0
#define X14     14.0/511.0

// filename for ionization model

#define REFLION_FILE "reflionx.mod"


// Main model routine

void rfxconv(const RealArray& energyArray, const RealArray& params,
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
   //  number of model parameters:5
   //  1: scale, scaling factor for reflection; if <0, no direct component
   //     (scale=1 for isotropic source above disk)
   //  2: redshift, z
   //  3: iron abundance relative to the solar abundances
   //  4: cosine of inclination angle
   //  5: log10(Xi) for ionization
   //
   // also uses the following values which can be set using xset
   //   RFXCONV_MAX_E      maximum energy for which to calculate the output spectrum
   //   RFXCONV_PRECISION  precision used in adaptive Gauss-Kronrod quadrature of Greens'
   //                      function integral
   //   RFXCONV_DIR        directory to use for ionization model filenames instead of
   //                      the standard modelData directory
   //  Based on reflct.cxx
   
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
   Real FeAbund  = params[2];
   Real cosIncl = params[3];
   if ( cosIncl < 0.05 ) cosIncl = 0.05;
   if ( cosIncl > 0.95 ) cosIncl = 0.95;
   Real Xi   = pow(10.0,params[4]);

   // parameters that will be used to get the spectrum from the reflion table model

   RealArray tableParams(4);
   tableParams[0] = FeAbund;
   tableParams[2] = Xi;
   tableParams[3] = 0.0;

   // calculate the Compton reflection (non-relativistic and relativistic) using 
   // the Ross & Fabian table models for the ionization.

   calcReflFlux("RFXCONV", tableParams, Scale, cosIncl, Abund, 0.0, zshift, 
		energyArray, flux);
                   
}

// *********************************************************************************
// Routine to calculate the reflected flux.

void calcReflFlux(const string modelName, RealArray& tableParams, Real Scale, 
		  Real cosIncl, Real Abund, Real cutoff, Real zshift, 
		  const RealArray& energyArray, RealArray& flux)
{

  static Real xnor(pow(10.0/511.0, 3.0));

  size_t Ne = energyArray.size() - 1;
  RealArray X(Ne);

  //  x is in the source frame and in units of m_e c^2.

  Real factor = zshift/2.0/511.0;
  for (size_t i=0; i<Ne; ++i) X[i] = (energyArray[i+1]+energyArray[i]) * factor;

  // This model will not work and may seg fault if the input energy range does not include
  // 2 - 14 keV in the source frame.

  if ( X[0] > X2 || X[Ne-1] < X14 ) {
    std::ostringstream err;
    err << modelName << " requires that the energy range include at least 2 - 14 keV in the source frame.\n"
	<< "The current energy range is " << energyArray[0] << " to " << energyArray[Ne-1] 
	<< " and the redshift is " << zshift-1 << ".\n"
	<< "No reflection will be applied. Use the energies command to extend the energy range.";
    xs_write(const_cast<char*>(err.str().c_str()), 10);
    return;
  }

  // Generate the spinc array from the input photar array.
   
  RealArray Spinc(Ne);
  for (size_t i=0; i<Ne; ++i) {
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

  // check whether the user has defined a maximum X that they want calculated

  Real Xmax(-1.0);
  string pname = modelName + "_MAX_E";
  string pvalue = FunctionUtility::getModelString(pname);
  if ( pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY() ) {
    std::istringstream iss(pvalue);
    if ( !(iss >> Xmax) || !iss.eof() ) {
      std::ostringstream err;
      err << "Invalid " << modelName + "_MAX_E value: " << pvalue;
      xs_write(const_cast<char*>(err.str().c_str()), 10);
    }
  }
  Xmax /= 511.0;

  // Include the direct flux

  RealArray Spref(Ne);
  RealArray Sptot(Ne);
  Sptot = 0.0;

  if ( Scale >= 0.0 ) Sptot += Spinc;
   
  if ( Scale != 0.0 ) {

    // Find average index of the input spectrum between 2 and 10 keV
    // Since we test above for the input spectrum containing the 2-10 keV
    // interval this should not produce a problem

    bool good;
    Real index = estimatePowerLawSlope(X2, X10, X, Spinc, good);
    if ( !good ) return;

    // Calculate reflection spectrum from the reflionx.mod or xillver table models

    if ( modelName == "RFXCONV" ) {

      calcRossRefl(modelName, energyArray, tableParams, index, Spref);
      cutoff = 300.0;

    } else if ( modelName == "XILCONV" ) {

      calcXilRefl(modelName, energyArray, tableParams, index, Spref);

    }
    if ( Spref.size() == 0 ) Spref = Spinc;

    // calculate the input continuum array which was used to generate this
    // reflected spectrum

    RealArray contArray;
    calcCutoffPowerLaw(energyArray, index, cutoff, false, contArray);

    // generate the reflected spectrum

    Spref *= Spinc/contArray;

    // set the maximum X value required in the reflected component

    Real xrefmax = 1/(YMIN+1-sqrt(1-(cosIncl-.05)*(cosIncl-.05)));
    if ( Xmax < 0.0 ) {
      Xmax = xrefmax;
    } else {
      if ( xrefmax < Xmax ) Xmax = xrefmax;
    }

    // Find the index corresponding to X14

    size_t j14 = 0;
    while ( X[j14] < X14 && j14 < X.size() ) j14++;
    j14--;
    
    // Estimate the average photon index of the reflection between
    // X12 and X14. Since we test above for the input spectrum containing the 12-14 keV
    // interval this should not produce a problem

    Real gamma_ref = estimatePowerLawSlope(X12, X14, X, Spref, good);
    if ( !good ) return;

    // Now iterate over reflection models over 12-14 keV finding one with 
    // the same slope. Use Newton-Raphson since there is a nice monotonic
    // relation between xnor and slope.

    RealArray Sprefp(Ne);

    size_t count = 0;

    MZCompRefl greenir;

    // Initial calculation of the slope based on the current value of xnor

    Real xnor1 = xnor;

    // Greens' fn calculation for this albedo

    greenir.CalcReflection(modelName, cosIncl, xnor1, X14, X, Spinc, Sprefp);
	
    // calculate the 12-14 keV index. Since we test above for the input spectrum 
    // containing the 12-14 keV interval this should not produce a problem


    Real gamma_pex1 = estimatePowerLawSlope(X12, X14, X, Sprefp, good);
    if ( !good ) return;

    // search until difference between reference and calculated slope is
    // less than 0.1% relative.
    
    Real gamma_pex = gamma_pex1;
    
    while ( abs(gamma_ref-gamma_pex)/abs(gamma_ref) > 1.0e-3 && count < 40 ) {
      
      count++;

      // calculate derivative for current value of xnor

      xnor1 = xnor + 1.0e-7;
      greenir.CalcReflection(modelName, cosIncl, xnor1, X14, X, Spinc, Sprefp);
      gamma_pex1 = estimatePowerLawSlope(X12, X14, X, Sprefp, good);
      if ( !good ) return;

      Real dgamma_dxnor = 1.0e7 * (gamma_pex1-gamma_pex);

      // calculate new value of xnor and evaluate slope
      // ensure that xnor > 0.0

      xnor -= (gamma_pex-gamma_ref)/dgamma_dxnor;
      while ( xnor < 0.0 ) {
	xnor += 0.1*(gamma_pex-gamma_ref)/dgamma_dxnor;
      }

      greenir.CalcReflection(modelName, cosIncl, xnor, X14, X, Spinc, Sprefp);
      gamma_pex = estimatePowerLawSlope(X12, X14, X, Sprefp, good);
      if ( !good ) return;

    }
    
    // Run reflection for this albedo over entire energy range
    
    greenir.CalcReflection(modelName, cosIncl, xnor, Xmax, X, Spinc, Sprefp);

    // Renormalize below 14 keV by the difference at 14 keV

    Real factor = Sprefp[j14]/Spref[j14];

    if ( factor != factor ) {
      ostringstream msg;
      msg << "NaN detected : " << factor << " " << Sprefp[j14] << " " << Spref[j14] 
	    << " " << gamma_ref << " " << gamma_pex << " " << X[0]*511.0 << " "
	    << X[Ne-1]*511.0 << " " << xnor << endl;
      FunctionUtility::xsWrite(msg.str(),5);
    }
    
    for (size_t j=0; j<=j14; j++) Spref[j] *= factor;

    // Fill up the rest of the Spref array from the CalcReflection output

    for (size_t j=j14+1; j<Ne; j++) Spref[j] = Sprefp[j];

    // Add in reflected spectrum

    Sptot += abs(Scale) * Spref;

  }

   // convert back to xspec internal units
   
   for (size_t i=0; i<Ne; ++i)
   {
     Real avgE = (energyArray[i+1]+energyArray[i])/2.;
     flux[i] = Sptot[i]*(energyArray[i+1]-energyArray[i])/(avgE*avgE);
   }

  return;

}

// *****************************************************************************
// Use the Ross & Fabian reflionx ionized disk table model to calculate
// reflected flux

void calcRossRefl(const string modelName, const RealArray& energyArray,
		  RealArray& tableParams, Real index, RealArray& Spref)
{

  Spref.resize(energyArray.size()-1);
  Spref = 0.0;

  // Average index of the input spectrum between 2 and 10 keV

  tableParams[1] = index;

  // Table model is only tabulated for slopes of 1.4 to 3.3 so if calculated slope
  // lies outside this range then use extrema and write warning message.

  if ( tableParams[1] < 1.4 ) {
    stringstream msg;
    msg << modelName << ": A 2-10 keV slope of " << tableParams[1] << " is outside the tabulated range - setting to 1.4" << endl;
    FunctionUtility::xsWrite(msg.str(),5);
    tableParams[1] = 1.4;
  } else if ( tableParams[1] > 3.3 ) {
    stringstream msg;
    msg << modelName << ": A 2-10 keV slope of " << tableParams[1] << " is outside the tabulated range - setting to 3.3" << endl;
    FunctionUtility::xsWrite(msg.str(),5);
    tableParams[1] = 3.3;
  }

  // Check whether a directory has been set for the input table files.
  // If not use the standard location

  string pname = modelName + "_DIR";
  string DirName(FunctionUtility::getModelString(pname));
  
  if ( DirName.length() == 0 || DirName == FunctionUtility::NOT_A_KEY() )
    DirName =  FunctionUtility::modelDataPath();

  string FileName =  DirName + REFLION_FILE;

  RealArray fluxErrArray;

  int ifl = 0;
  tableInterpolate(energyArray, tableParams, FileName, ifl, Spref, fluxErrArray, " ", "add", true);
  if ( Spref.size() == 0 ) {
    stringstream msg;
    msg << modelName << ": Failed to read " << FileName << endl;
    FunctionUtility::xsWrite(msg.str(),5);
  }

}

// *****************************************************************************
// Estimate the power-law slope between Xlow and Xhigh

Real estimatePowerLawSlope(Real Xlow, Real Xhigh, RealArray& X, RealArray& Y, bool& good)
{
  Real slope(0.0);
  good = true;

  if ( Xhigh > X[0] && Xlow < X[X.size()-1] ) {

    size_t jlow = 0;
    while ( jlow < X.size() && X[jlow] < Xlow ) jlow++;

    size_t jhigh = jlow;
    while ( jhigh < X.size() && X[jhigh] < Xhigh ) jhigh++;
    jhigh--;

    if ( jhigh > jlow ) {
      slope = log10(Y[jhigh]/(X[jhigh]*X[jhigh])) - log10(Y[jlow]/(X[jlow]*X[jlow]));
      slope /= log10(X[jlow]) - log10(X[jhigh]);
    } else {
      good = false;
    }

  } else {

    good = false;

  }

  return(slope);
}
