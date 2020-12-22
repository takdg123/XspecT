//c ------------------------------------------------------------------------
//c
//c  ** NOTE - FOR DETECTORS WITH LIMITED RESPONSES, THE ARRAY OF
//c     ENERGIES SHOULD BE EXTENDED WITH THE 'ENERGIES' COMMAND. 
//c
//c   FOR EXAMPLE, WITH RXTE:
//c   >    energies 0.2 50 500 lin 
//c   
//c
//c   This routine is designed to give a physically motivated
//c   prescription to model Comptonization of an input seed spectrum.
//c   Photon number is preserved in the output.  
//c
//c     USER INPUT:
//c
//c     param(1) - Photon index of the Comptonized component in the spectrum. 
//c     param(2) - Normalization of the spectrum. 
//c     param(3) - Flag - >0 for up-scattering only (SIMPL-1), <=0 for both
//c                 up-scattering and down-scattering (SIMPL-2).
//c------------------------------------------------------------------------
//c     XSPEC INPUT
//c     
//c     2 Vectors:
//c     -----------------
//c     1) Energy
//c     2) Flux (Energy)
//c
//c
//c     EAR(0:NE) - array of observed energy bins in keV
//c     NE        - # observed energy bins
//c     PARAM     - list of parameters
//c     PHOTAR() - an array of photon flux in each bin of energy PHOTAR(NE)
//c     PHOTER() - formal output for errors
//c
//c
//c
//c
//c  EXAMPLE:
//c      model phabs(simpl(diskbb))
//c
//c==========================================================================

#include <XSFunctions/functionMap.h>
#include "xsTypes.h"

void simpl (const RealArray& energyArray, 
            const RealArray& params, 
            int spectrumNumber,
            RealArray& fluxArray, 
            RealArray& fluxErrArray,
            const string& initString)
{


  //      integer i,j,k
  //      double precision norm,gamma,gnormUP,gnormDN,gamma1,gamma2
  //      real ear(0:ne),param(3),photar(ne),photer(ne)
  //      real tmparr(ne),enavgam1(ne),enavgam2(ne)
  //      real engam1(0:ne),engam2(0:ne)

  Real norm (params[1]);
  Real gamma (params[0]);
  Real scatteringflag (params[2]);
  const size_t nE = energyArray.size();
  const size_t nBins = nE-1;

  //  WORKAROUND to avoid pole at gamma=1

  if ( gamma == 1 ) gamma = 1.001;

  Real gamma1 (gamma-1.);
  Real gamma2 (gamma+2.);

  Real gnormUP ((gamma+2.)/(1.+2.*gamma));
  Real gnormDN ((gamma-1.)/(1.+2.*gamma));

  //  Initialize arrays

  RealArray engam1 (nE);
  RealArray enavgam1 (nBins);
  RealArray tmparr (nBins);

  engam1 = pow(energyArray,-gamma1);
  for (size_t i=0; i<nBins; i++) {
    tmparr[i] = 0.0;
    enavgam1[i] = pow(0.5*(energyArray[i]+energyArray[i+1]),gamma1);
  }

  if ( scatteringflag > 0 ) {

    // !UP-SCATTERING ONLY

    for (size_t i=0; i<nBins; i++) {

      //  Calculate scattered contribution to the same bin
      tmparr[i] += fluxArray[i]*(1.-enavgam1[i]*engam1[i+1]);

      //  Loop over all higher bins and calculate scattered contributions
      for (size_t j=i+1; j<nBins; j++) {
	tmparr[j] += enavgam1[i]*(engam1[j]-engam1[j+1])*fluxArray[i];
      }

    }

  } else {

    // UP-SCATTERING & DOWN-SCATTERING 

    RealArray engam2 (nE);
    RealArray enavgam2 (nBins);

    engam2 = pow(energyArray,gamma2);
    for (size_t i=0; i<nBins; i++) {
      enavgam2[i]=pow(0.5*(energyArray[i]+energyArray[i+1]),-gamma2);
    }


    //  Loop over energy bins

    for (size_t i=0; i<nBins; i++) {

      // Calculate scattered contribution to the same bin
      tmparr[i] += fluxArray[i]*(gnormUP*(1.-enavgam1[i]*engam1[i+1])
			       + gnormDN*(1.-enavgam2[i]*engam2[i+1]));

      // Loop over all bins and calculate scattered contributions

      for (size_t j=0; j<nBins; j++) {

	if ( j < i ) {
	  tmparr[j] += enavgam2[i]*(engam2[j+1]-engam2[j])*fluxArray[i]*gnormDN;
	} else if ( j > i ) {
	  tmparr[j] += enavgam1[i]*(engam1[j]-engam1[j+1])*fluxArray[i]*gnormUP;
	}

      }

    }

  }

  fluxArray = (1.-norm)*fluxArray + norm*tmparr;

}
