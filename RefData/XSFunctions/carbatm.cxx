// Carbon model atmosphere model by V. Suleimanov
// (Suleimanov et al. A&A, 2013)
// Implemented by D. Klochkov 2013
// Updated by D. Klochkov in 2016 to use new tables
// (Suleimanov et al. A&A 2016)
// email: klochkov@astro.uni-tuebingen.de
// modified by kaa 5/9/17 to use a single FITS file for the input table
//
// parameter[0] - effective (unredshifted) temperature of the 
//                neutron star surface (in MK)
//                T = 1.0..4.0
// parameter[1] - neutron star gravitational mass (in solar mass)
// parameter[2] - neutron star radius (in km)
//
//////////////////////////////////////////////////



#include "xsTypes.h"
#include <XSFunctions/functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include "BinarySearch.h"
#include <CCfits/CCfits>
#include <memory>

void calcNSatm(const RealArray& energy, const RealArray& parameter, 
	       const string inputFilename, RealArray& flux);


void carbatm(const RealArray& energy, const RealArray& parameter, int spectrum, 
	     RealArray& flux, RealArray& fluxError, const string& init)
{

  //..path to the directory with model tables:
  // if variable CARBATM is set, use that path. If not, use
  // $HEADAS/../spectral/modelData
  string directory = FunctionUtility::getModelString("CARBATM");
  if(directory.compare("$$NOT$$") == 0) { // CARBATM variable is not set
    directory = FunctionUtility::modelDataPath();
  }
  //..initialize filename
  const string inputFilename(directory+"spC.fits");

  calcNSatm(energy, parameter, inputFilename, flux);
  fluxError.resize(0);
  return;
}

void calcNSatm(const RealArray& energy, const RealArray& parameter, 
	       const string inputFilename, RealArray& flux)
{
  using namespace std;
  using namespace CCfits;

  Real T            =parameter[0];
  const Real NSmass =parameter[1];
  const Real NSrad  =parameter[2];

  static int Ntemps;     //..number of tabulated temperatures
  static RealArray Ttab; //..tabulated temperatures in MK
  static int Nenergies;   //..number of energy bins in tabulated spectra
  static RealArray Elow,Ehigh;//..model energy bins
  static int Nlogg;       //..number of tabulated logg
  static RealArray loggtab;//..tabulated values of log g (in cgs)
  static vector<vector<RealArray> >  Flux;

  Real logg;             //..common logarithm of surface gravity (in cgs)
  Real zfactor;          //..1 + gravitational redshift at NS surface

  Real lowlogg, highlogg;//..bracketing values of logg, i.e. the closest
                         // values of tabulated logg around the current value
  int  ilogg;            //..index of the lower bracketing logg file
  Real lowT, highT;      //..bracketing values of T
  int  iT;               //..index of the lower bracketing T
  Real checkFactor=1.0;  //..The final flux is multiplied by this factor, 
                         // which is normally 1.0 but goes up and generates
                         // a warning message if log(g) is
                         // outside the allowed range. This is to avoid the
                         // fit routine to go deep outside the available log(g)
                         // range by "steppar" or "error" commands)
  static string lastInputFilename(" ");//..to check whether the input filename has changed
                                       // will automatically catch the first invocation


  //..calculate (1+z) factor and redshift and logg
  if(NSmass/(0.337*NSrad) < 0.99) {
    zfactor = 1./sqrt(1. - NSmass/(0.337*NSrad));
  } else {
    zfactor = 1./sqrt(1. - 0.99);
    cout<<"!! Warning: the calculated Z is above the limit (9.0)."<<endl;
    cout<<"The results are not trustable for this M and R ("<<NSmass<<", "<<NSrad<<")!!"<<endl;
  }
  logg = 16.125 + log10(NSmass*zfactor/(NSrad*NSrad));
  
  //cout<<"logg= "<<logg<<", (1+z)= "<<zfactor<<", M= "<<NSmass<<", R= "<<NSrad<<", T= "<<T<<endl;


  // if first time through or filename has changed then load the arrays from the input file

  if ( inputFilename != lastInputFilename ) {

    FunctionUtility::xsWrite("Reading NS atmosphere input data from "+inputFilename,25);

    // open the input FITS file

    std::unique_ptr<FITS> pInfile((FITS*)0);
    try {
      pInfile.reset(new FITS(inputFilename, Read));
    } catch(...) {
      FunctionUtility::xsWrite("Failed to read "+inputFilename,5);
      return;
    }

    // read the extension with energies

    try {
      ExtHDU& Etable = pInfile->extension("ENERGIES");
      Nenergies = Etable.rows();
      Etable.column("ELOW").read(Elow,1,Nenergies);
      Etable.column("EHIGH").read(Ehigh,1,Nenergies);
    } catch(...) {
      FunctionUtility::xsWrite("Failed to find ENERGIES extension in "+inputFilename,5);
      return;
    }

    // read the extension with the logg values and extension names
    vector<string> extname;
    try {
      ExtHDU& LGtable = pInfile->extension("LOGGVALS");
      Nlogg = LGtable.rows();
      LGtable.column("LOGG").read(loggtab,1,Nlogg);
      LGtable.column("LOGGNAME").read(extname,1,Nlogg);
    } catch(...) {
      FunctionUtility::xsWrite("Failed to find LOGGVALS extension in "+inputFilename,5);
      return;
    }

    // read the extension with the temperature values
    try {
      ExtHDU& Ttable = pInfile->extension("TEMPVALS");
      Ntemps = Ttable.rows();
      Ttable.column("TEMP").read(Ttab,1,Ntemps);
    } catch(...) {
      FunctionUtility::xsWrite("Failed to find TEMPVALS extension in "+inputFilename,5);
      return;
    }

    // loop round the extensions with different log g loading the fluxes

    Flux.resize(Nlogg);
    for (size_t ig=0; ig<(size_t)Nlogg; ig++) {
      try {
	ExtHDU& Ftable = pInfile->extension(extname[ig]);
	Ftable.column("FLUX").readArrays(Flux[ig],1,Ntemps);
      } catch(...) {
	FunctionUtility::xsWrite("Failed to find "+extname[ig]+" extension in "+inputFilename,5);
	return;
      }
    }

    lastInputFilename = inputFilename;

  }


  //..find bracketing logg and corresponding index
  ilogg = Numerics::BinarySearch(loggtab, logg);
  if(ilogg < 0) {
    cout<<"!! Warning: log(g)="<<logg<<" is outside the limits for which model atmospheres are calculated ("<< loggtab[0] << "..." << loggtab[Nlogg-1] << ")"<<endl;
    cout<<"The results are not trustable for this M and R ("<<NSmass<<", "<<NSrad<<")!!"<<endl; 
    if(ilogg == -1) {
      checkFactor = 1.+1.e4*(loggtab[0]-logg)/logg;
      logg = loggtab[0];
      ilogg = 0;
    } else if(ilogg == -2) {
      checkFactor = 1.+1.e4*(logg-loggtab[Nlogg-1])/loggtab[Nlogg-1];
      logg = loggtab[Nlogg-1];
      ilogg = Nlogg-2;
    }
  }
  lowlogg = loggtab[ilogg];
  highlogg= loggtab[ilogg+1];

  //..find bracketing T and corresponding index
  iT = Numerics::BinarySearch(Ttab, T);
  if(iT < 0) {
    cout<<"!! Warning: T="<<T<<" is outside the limits for which model atmospheres are calculated ("<< Ttab[0] << "..." << Ttab[Ntemps-1] << ")"<<endl;
    cout<<"The results are not trustable for this temperature!!"<<endl;
    if(iT == -1) {
      checkFactor = 1.+1.e4*(Ttab[0]-T)/T;
      T = Ttab[0];
      iT = 0;
    } else if(iT == -2) {
      checkFactor = 1.+1.e4*(T-Ttab[Ntemps-1])/Ttab[Ntemps-1];
      T = Ttab[Ntemps-1];
      iT = Ntemps-2;
    }
  }
  lowT = Ttab[iT];
  highT= Ttab[iT+1];

  // and interpolate between bracketing spectra to calculate flux
  Real A = 1./((highlogg-lowlogg)*(highT-lowT));
  Real B = (highlogg-logg)*(highT-T);
  Real C = (logg -lowlogg)*(highT-T);
  Real D = (highlogg-logg)*(T- lowT);
  Real E = (logg -lowlogg)*(T- lowT);
  RealArray fluxInterp(Nenergies);  //..bilinearly interpolated model flux
  fluxInterp = A*( Flux[ilogg][iT]*B + Flux[ilogg+1][iT]*C+
		   Flux[ilogg][iT+1]*D + Flux[ilogg+1][iT+1]*E );
  
  //..correct interpolated model spectrum for redshift
  RealArray zElow = Elow/zfactor;
  RealArray zEhigh = Ehigh/zfactor;
  fluxInterp /= zfactor;
  
  //..calculate flux in provided energy bins
  size_t N(energy.size()-1);
  flux.resize(N);
  size_t i=0;
  for(size_t k=0; k<N; ++k) { //..loop over energy bins in observed spectrum
    if(energy[k+1] < zEhigh[Nenergies-1]) {
      while(zEhigh[i] < energy[k]) {i++;}
      if(zEhigh[i] < energy[k+1]) {
	flux[k] = fluxInterp[i]*(zEhigh[i]-energy[k])/(zEhigh[i]-zElow[i]);
	i++;
	while(zEhigh[i] < energy[k+1]) {
	  flux[k]+=fluxInterp[i];
	  i++;
	}
	flux[k] += fluxInterp[i]*(energy[k+1]-zElow[i])/(zEhigh[i]-zElow[i]);
      } else {
	flux[k] = fluxInterp[i]*(energy[k+1]-energy[k])/(zEhigh[i]-zElow[i]);
      }
    } else { flux[k] = 0.; }
  }

  //..multiply flux by NSrad^2
  flux *= NSrad*NSrad*checkFactor;

  return;
}
