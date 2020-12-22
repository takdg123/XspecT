//     Aug 99    Plane-parallel shock models with unequal ion and electron 
//               temperatures. May be used to model a shock section.
//		    additive model for XSPEC to use with general models.
//		    parameters (param):
//               (0):    postshock temperature t(keV)
//               (1):    postshock electron temperature t_e(keV)
//               (2):    hydrogen abundance (switches on and off free-free
//                       continuum
//               (3..14): abundances of He,C,N,O,Ne,Mg,Si,S,Ar,Ca,Fe,Ni
//                                  with respect to solar values
//               (15):   ionization parameter taul(cm^-3 s)
//               (16):   ionization parameter tauu(cm^-3 s)
//               (17):   redshift z
//
//               Parameter (2) should be usually set to one. Its only use
//               in the code is to set hydrogen abundance to zero; in this
//               case this parameter should be set to zero. It is not the
//               density parameter.

#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/XSstream.h>
#include <IonBalNei.h>
#include <Aped.h>

using namespace std;
using namespace XSutility;

// prototype for routine to get the trace abundance

Real getNEItraceAbund(const IntegerVector& Z, const RealArray& abundance, 
		      const string& callr);



void vnpshock(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{
  const size_t nZin = 12;
  const int Zinfo[12] = {2,6,7,8,10,12,14,16,18,20,26,28};
  IntegerVector Zinput(&Zinfo[0],&Zinfo[0]+12);

  RealArray abundance(nZin);
  for (size_t i=0; i<nZin; i++) abundance[i] = params[i+3];
  Real TraceAbund = getNEItraceAbund(Zinput, abundance, "vnpshock");

  // set parameters to pass through to vvnpshock

  RealArray vvpar(35);
  vvpar[0] = params[0];
  vvpar[1] = params[1];
  vvpar[2] = params[2];
  for (size_t i=3; i<32; i++) vvpar[i] = TraceAbund;
  for (size_t i=0; i<nZin; i++) vvpar[Zinput[i]+1] = abundance[i];
  vvpar[32] = params[15];
  vvpar[33] = params[16];
  vvpar[34] = params[17];

  vvnpshock(energyArray, vvpar, spectrumNumber, flux, fluxErr, initString);

  return;
}
