//    Mar 09   Translated from xseq.f (CG)
//    May 02    added Ar
//    Aug 23,99 subroutine hunt used instead of ibine for distribution of
//              lines into bins.
//    feb 95 Equlibrium ionization models 
//		  additive model for XSPEC to use with general models.
//		  parameters (param):
//             [0]:     temperature t(keV)
//             [1..12]: abundances of He,C,N,O,Ne,Mg,Si,S,Ar,Ca,Fe,Ni
//                              with respect to solar values
//             [13]:    redshift z

#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <CCfits/CCfits>
#include <Aped.h>
#include <IonBalNei.h>

using namespace std;
using namespace XSutility;


void vequil(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
  const size_t nZ = 13;
  const int Zinfo[nZ] = {1, 2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28};
  const IntegerVector Zinput(&Zinfo[0],&Zinfo[0]+13);

  vector<RealArray> IonFrac;
  calcCEIfractions(params[0], Zinput, IonFrac);

  // set up Zinput and abundance arrays

  RealArray abundance(nZ);
  abundance[0] = 1.0;
  for (size_t i=1; i<nZ; i++) abundance[i] = params[i];

  // calculate output spectrum

  Real Redshift(params[13]);

  int status = calcNEISpectrum(energyArray, Zinput, abundance, Redshift, 
			       params[0], IonFrac, false, 0.0, flux, fluxErr);
  if ( status != 0 ) {
    std::stringstream oss;
    oss << "calcNEISpectrum failed in vequil: status = " << status << "\n";
    oss << "set chatter to 25 and retry to get additional diagnostics.\n";
    xs_write(const_cast<char*>(oss.str().c_str()),10);
  }

  return;
}
