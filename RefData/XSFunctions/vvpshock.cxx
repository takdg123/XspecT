//c     Aug 9,99  Nonequlibrium ionization models with constant electron
//c               temperature, with a distribution of ionization timescales
//c               appropriate for a plane parallel shock or its section.
//c               additive model for XSPEC to use with general models.
//c               parameters (param):
//c               (0):    temperature t(keV)
//c               (1):    hydrogen abundance (switches on and off free-free
//c                       continuum
//c               (2..30): abundances of all elements from Z=2 to Z=30
//c                                  with respect to solar values
//c                     ionization range bounded by:
//c               (31):   ionization parameter taul(cm^-3 s)
//c               (32):   ionization parameter tauu(cm^-3 s)
//c               (33):   redshift z
//c
//c               Parameter (1) should be usually set to one. Its only use
//c               in the code is to set hydrogen abundance to zero; in this
//c               case this parameter should be set to zero. It is not the
//c               density parameter.

// based on doxsneqs.f code

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

void vvpshock(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
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

   static bool isFirst = true;
   static RealArray tauSave;
   static RealArray weightSave;

   // if necessary calculate the tau and weight arrays

   RealArray tau;
   RealArray weight;

   if ( !identicalArrays(tau,tauSave) || !identicalArrays(weight,weightSave) ||
	isFirst ) {

     Real tauLow = params[31];
     Real tauHigh = params[32];

     // loop over zones setting up taus and weights

     if ( tauLow == 0.0 ) {

       size_t nzones = 200;
       tau.resize(nzones);
       weight.resize(nzones);
       
       Real taumin = 1.0e8;
       Real deltau = log10(tauHigh/taumin)/(Real)(nzones-1);
       Real tauiml = 0.0;
       for (size_t i=0; i<nzones; i++) {
	 Real taui = taumin * pow(10.0, i*deltau);
	 weight[i] = (taui - tauiml)/tauHigh;
	 tau[i] = 0.5*(taui+tauiml);
	 tauiml = taui;
       }

     } else if ( tauHigh == tauLow ) {

       tau.resize(1);
       weight.resize(1);
       tau[0] = tauHigh;
       weight[0] = 1.0;

     } else {

       size_t nzones = 40;
       tau.resize(nzones);
       weight.resize(nzones);

       Real deltau = (tauHigh - tauLow)/(Real)nzones;
       Real tauiml = tauLow;
       for (size_t i=0; i<nzones; i++) {
	 Real taui = tauLow + i * deltau;
	 weight[i] = 0.025;
	 tau[i] = 0.5*(taui+tauiml);
	 tauiml = taui;
       }

     }

     tauSave.resize(tau.size());
     tauSave = tau;
     weightSave.resize(weight.size());
     weightSave = weight;

   }

   RealArray tkeV(2);
   tkeV[0] = params[0];
   tkeV[1] = params[1];

   // calculate the ion fractions

   vector<RealArray> IonFrac;
   calcNEIfractions(tkeV, tau, weight, Zinput, IonFrac);

   // calculate output spectrum

   Real Redshift(params[33]);

   int status = calcNEISpectrum(energyArray, Zinput, abundance, Redshift, 
				params[0], IonFrac, false, 0.0, flux, fluxErr);
   if ( status != 0 ) {
     std::stringstream oss;
     oss << "calcNEISpectrum failed in vvpshock: status = " << status << "\n";
     oss << "set chatter to 25 and retry to get additional diagnostics.\n";
     xs_write(const_cast<char*>(oss.str().c_str()),10);
   }

   isFirst = false;

   return;

}

