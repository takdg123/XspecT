//c     Aug 9,99  Nonequlibrium ionization models with constant electron
//c               temperature, with a distribution of ionization timescales
//c               appropriate for a plane parallel shock or its section.
//c               additive model for XSPEC to use with general models.
//c               parameters (param):
//c               (1):    temperature t(keV)
//c               (2):    hydrogen abundance (switches on and off free-free
//c                       continuum
//c               (3..14): abundances of He,C,N,O,Ne,Mg,Si,S,Ar,Ca,Fe,Ni
//c                                  with respect to solar values
//c                     ionization range bounded by:
//c               (15):   ionization parameter taul(cm^-3 s)
//c               (16):   ionization parameter tauu(cm^-3 s)
//c               (17):   redshift z
//c
//c               Parameter (2) should be usually set to one. Its only use
//c               in the code is to set hydrogen abundance to zero; in this
//c               case this parameter should be set to zero. It is not the
//c               density parameter.

// based on doxsneqs.f code

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



void vpshock(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
  const size_t nZin=12;
  const int Zinfo[nZin] = {2,6,7,8,10,12,14,16,18,20,26,28};
  IntegerVector Zinput(&Zinfo[0], &Zinfo[0]+12);

  RealArray abundance(nZin);
  for (size_t i=0; i<12; i++) abundance[i] = params[i+2];
  Real TraceAbund = getNEItraceAbund(Zinput, abundance, "vpshock");

  // set parameters to pass through to vvpshock

  RealArray vvpar(34);
  vvpar[0] = params[0];
  vvpar[1] = params[1];
  for (size_t i=2; i<31; i++) vvpar[i] = TraceAbund;
  for (size_t i=0; i<nZin; i++) vvpar[Zinput[i]] = params[i+2];
  vvpar[31] = params[14];
  vvpar[32] = params[15];
  vvpar[33] = params[16];

  vvpshock(energyArray, vvpar, spectrumNumber, flux, fluxErr, initString);

  return;
}
