// This version of the GNEI model using the Aped class.

#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/XSstream.h>
#include <IonBalNei.h>
#include <Aped.h>

using namespace std;
using namespace XSutility;

// prototype for routine to get the trace abundance

Real getNEItraceAbund(const IntegerVector& Z, const RealArray& abundance, 
		      const string& callr);




void vgnei(const RealArray& energyArray, const RealArray& params,
	   int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	   const string& initString)
{
   //
   //  aug 95   Non-equlibrium ionization models with two electron
   //           temperatures (final and ionization-timescale averaged)
   //		additive model for XSPEC to use with general models.
   //		parameters (params):
   //           (0):    temperature t(keV)
   //           (1):    hydrogen abundance (switches on and off free-free
   //                   continuum
   //		(2..13): abundances of He,C,N,O,Ne,Mg,Si,S,Ar,Ca,Fe,Ni
   //                                  with respect to solar values
   //           (14):   ionization timescale tau (cm^-3 s^-1)
   //		(15):	tau-averaged temperature tav(keV)
   //		(16):	redshift z
   //
   //           Parameter (1) should be usually set to one. Its only use
   //           in the code is to set hydrogen abundance to zero; in this
   //           case this parameter should be set to zero. It is not the
   //           density parameter.
   //
   //           Calls vvgnei after setting the trace element abundances
   //           appropriately.

  const size_t nZin=12;
  const int Zinfo[nZin] = {2, 6, 7, 8, 10, 12, 14, 16, 18, 20, 26, 28};
  IntegerVector Zinput(&Zinfo[0],&Zinfo[0]+12);

  RealArray abundance(nZin);
  for (size_t i=0; i<nZin; i++) abundance[i] = params[i+2];
  Real TraceAbund = getNEItraceAbund(Zinput, abundance, "vgnei");

  // set parameters to pass through to vvgnei

  RealArray vvpar(34);
  vvpar[0] = params[0];
  vvpar[1] = params[1];
  for (size_t i=2; i<31; i++) vvpar[i] = TraceAbund;
  for (size_t i=0; i<nZin; i++) vvpar[Zinput[i]] = abundance[i];
  vvpar[31] = params[14];
  vvpar[32] = params[15];
  vvpar[33] = params[16];

  vvgnei(energyArray, vvpar, spectrumNumber, flux, fluxErr, initString);

  return;
}

//---------------------------------------------------------------------------

Real getNEItraceAbund(const IntegerVector& Z, const RealArray& abundance, 
		      const string& callr)
{
  string pname("NEI_TRACE_ABUND");
  string pvalue(FunctionUtility::getModelString(pname));
  Real TraceAbund(1.0);

  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    int ielt=-1;
    string celt = XSutility::lowerCase(pvalue);
    for (size_t i=0; i<TOTZ; ++i) if (celt == atomName[i]) ielt = i;

    if ( ielt != -1 ) {
      bool good = false;
      for (size_t i=0; i<Z.size(); i++) {
	if ( ielt == (Z[i]-1) ) good = true;
      }
      if ( good ) {
	TraceAbund = abundance[ielt];
      } else {
	string message = "Trace abundance cannot be set to that of " 
	  + pvalue 
	  + " since the abundance of that element is not a parameter in the model.";
	xs_write(const_cast<char*>(message.c_str()),10);
      }
    } else {
      std::istringstream TraceAbundstr(pvalue);
      if (!(TraceAbundstr >> TraceAbund) || !TraceAbundstr.eof()) {
	string message = "Failed to read value from " + TraceAbundstr.str() + ": setting trace abundances to Solar";
	xs_write(const_cast<char*>(message.c_str()),10);
      }
    }
  }

  std::ostringstream oss;
  oss << callr << " : TraceAbund = " << TraceAbund << std::endl;
  xs_write(const_cast<char *>(oss.str().c_str()), 25);

  return TraceAbund;

}
