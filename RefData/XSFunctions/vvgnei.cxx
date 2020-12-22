// This a version of the GNEI model using the Aped
// class.

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

// prototype for routine which calls Kazik Borkowski's code for old versions
// of the NEI models

int KBcalcNEISpectrum(const RealArray& energyArray, const IntegerVector& Zinput,
		      const RealArray& abundance, const Real Redshift,
		      const Real& Tinput, 
		      const vector<RealArray>& IonFrac,
		      const bool qtherm, const Real velocity,
		      RealArray& fluxArray, RealArray& fluxErrArray);




void vvgnei(const RealArray& energyArray, const RealArray& params,
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
   //		(2..30): abundances of all elements from Z=2 to Z=30
   //                                  with respect to solar values
   //           (31):   ionization timescale tau (cm^-3 s^-1)
   //		(32):	tau-averaged timescale tav(keV)
   //		(33):	redshift z
   //
   //           Parameter (1) should be usually set to one. Its only use
   //           in the code is to set hydrogen abundance to zero; in this
   //           case this parameter should be set to zero. It is not the
   //           density parameter.
   //


   vector<RealArray> IonFrac;
   RealArray tavkeV(2);
   RealArray tau(1);
   RealArray weight(1);
   weight[0] = 1.0;
   tau[0] = params[31];
   tavkeV[0] = params[32];
   tavkeV[1] = params[32];

   // for speed reasons set Zinput and abundance arrays only for elements
   // with non-zero abundance

   int nZ = 0;
   for (size_t i=0; i<TOTZ; i++) {
     if ( params[i+1] > 0.0 ) nZ++;
   }
   IntegerVector Zinput(nZ);
   RealArray abundance(nZ);
   int iZ = 0;
   for (size_t i=0; i<TOTZ; i++) {
     if ( params[i+1] > 0.0 ) {
       Zinput[iZ] = i+1;
       abundance[iZ] = params[i+1];
       iZ++;
     }
   }

   // calculate the ion fractions

   calcNEIfractions(tavkeV, tau, weight, Zinput, IonFrac);

   //   string outstr = writeIonFrac(14, Zinput, IonFrac);
   //   xs_write(const_cast<char*>(outstr.c_str()),15);       

   // calculate output spectrum

   Real Redshift(params[33]);

   int status = calcNEISpectrum(energyArray, Zinput, abundance, 
				Redshift, params[0], IonFrac, false, 0.0, 
				flux, fluxErr);
   if ( status != 0 ) {
     std::stringstream oss;
     oss << "calcNEISpectrum failed in vvgnei: status = " << status << "\n";
     oss << "set chatter to 25 and retry to get additional diagnostics.\n";
     xs_write(const_cast<char*>(oss.str().c_str()),10);       
   }

   return;

}
