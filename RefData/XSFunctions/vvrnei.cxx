//    Vvrnei model
//	Recombining plasma model. This model requires NEIVERS 3.0 and above.
//		   parameters (param):
//              (0):    current temperature t(keV)
//              (1):    initial temperature t_i(keV)
//              (2):    hydrogen abundance (switches on and off free-free
//                      continuum
//              (3..31): abundances of all elements from Z=2 to Z=30
//                        with respect to solar values
//              (32):   ionization timescale tau (cm^-3 s)
//              (33):   redshift z
//
//              Parameter (2) should be usually set to one. Its only use
//              in the code is to set hydrogen abundance to zero; in this
//              case this parameter should be set to zero. It is not the
//              density parameter.


#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/XSstream.h>
#include <IonBalNei.h>
#include <Aped.h>
#include <CCfits/CCfits>

using namespace std;
using namespace XSutility;
using namespace CCfits;


void vvrnei(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
  // for speed reasons set Zinput and abundance arrays only for elements
  // with non-zero abundance

  int nZ = 0;
  for (size_t i=0; i<TOTZ; i++) {
    if ( params[i+2] > 0.0 ) nZ++;
  }
  IntegerVector Zinput(nZ);
  RealArray abundance(nZ);
  int iZ = 0;
  for (size_t i=0; i<TOTZ; i++) {
    if ( params[i+2] > 0.0 ) {
      Zinput[iZ] = i+1;
      abundance[iZ] = params[i+2];
      iZ++;
    }
  }

  // set the ion fractions for the initial equilibrium temperature

  vector<RealArray> InitialIonFrac;
  calcCEIfractions(params[1], Zinput, InitialIonFrac);

  //debug
  //    tcout << "Initial ion fractions:" << std::endl;
  //    for (size_t i=0; i<InitialIonFrac.size(); i++) {
  //      tcout << Zinput[i] << ": ";
  //      for (size_t j=0; j<InitialIonFrac[i].size(); j++) {
  //	tcout << InitialIonFrac[i][j] << ", ";
  //      }
  //      tcout << std::endl;
  //    }

  // calculate the ion fractions given the current temperature, ionization parameter
  // and initial ion fractions

  vector<RealArray> IonFrac;
  calcNEIfractions(params[0], params[32], Zinput, InitialIonFrac, IonFrac);

  // calculate the spectrum for the calculated ion fractions

  Real Redshift(params[33]);
  Real Tinput(params[0]);

  int status = calcNEISpectrum(energyArray, Zinput, abundance, Redshift, Tinput, 
				IonFrac, false, 0.0, flux, fluxErr);
  if ( status != 0 ) {
    std::stringstream oss;
    oss << "calcNEISpectrum failed in vvrnei: status = " << status << "\n";
    oss << "set chatter to 25 and retry to get additional diagnostics.\n";
    xs_write(const_cast<char*>(oss.str().c_str()),10);
  }

  return;
  
}



