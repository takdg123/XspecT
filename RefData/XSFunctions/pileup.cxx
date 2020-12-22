#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSUtil/Utils/XSutility.h>
#include <cfortran.h>
#include <cmath>

void pileup(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{

// John Davis' forward-folding algorithm for dealing with pile-up.
// This assumes that the input photar array has been multiplied by
// the effective areas. The model parameters are the frame time (sec)
// and the number of pile-ups to consider.

// Model parameters are :
//   1     frame time
//   2     max number of piled up photons
//   3     g0 - good event fraction
//   4     alpha - grade migration factor
//   5     psffrac - fraction of the counts to pile-up
//   6     nregions - number of independent regions
//   7     fracexpo - Chandra FRACEXPO keyword value

// Arguments :
//    energyArray     r        i: Energy ranges
//    params          r        i: Model parameters
//    spectrumNumber  i        i: Data set
//    flux            r      i/r: Model flux
//    fluxErr         r      i/r: Model flux errors
//    initString      s        i: Unused

// Calculates the expression

//    exp(-tau/g0 Int dE' S'(E') ) [(exp(tau F)-1)/(tau F)] o S'(E)

// where tau is the frame time, S' is the input, and F is the functional 
// defined by

//   F o f(E) = alpha Int_0^E dE' S'(E') f(E-E')

// The expression [(exp(tau F)-1)/(tau F)] is expanded in an infinite
// series

//   Sum_p  (tau F)^(p-1)/p!

  using namespace XSutility;

  int nE = flux.size();

  Real frame = params[0];
  int npiled  = static_cast<int>(floor(params[1]+0.5));
  Real g0 = params[2];
  Real alpha = params[3];
  Real psffrac = params[4];
  int nregions = static_cast<int>(floor(params[5]+0.5));
  Real fracexpo = params[6];

// This is the fraction of the input counts that are used to calculate the pile-up

  Real fract = psffrac/nregions/fracexpo;

// This is an offset we will require in the convolution in case ear(0) is
// non-zero. Note that we assume the energy array is linear and starts at
// a multiple of the energy bin size.

  int ioff = -energyArray[0]/(energyArray[1]-energyArray[0]);

// Set the output array and the temporary arrays for the single photon case
// and find the sum of the input array

  Real pnorm = 0.0;
  RealArray outar(nE);
  RealArray tmpar(nE);

  outar = flux * fract;
  tmpar = outar;
  pnorm = outar.sum();

  Real incounts = pnorm * frame;

// Calculate the exponential factor. 

  Real expfac = exp(-1 * frame * pnorm / g0);

  RealArray pfrac(npiled);

  pfrac[0] = 0.0;
  pfrac[1] = pnorm * frame * expfac;

// Loop round pile-ups

  Real factor = 1.0;

  for (int ipile=2; ipile<npiled; ipile++) {

// Do the convolution (should probably do an fft here but use BFI for now).

    RealArray conv(0.0,nE);

    for (int ie=0; ie<nE; ie++) {
      for (int je=0; je<ie+ioff; je++) {
	conv[ie] += flux[je]*fract*tmpar[ie-je+ioff];
      }
    }

// Add the latest iteration of the temporary array into the summation

    factor *= alpha * frame / ipile;

    outar += conv*factor;
    pfrac[ipile] = conv.sum()*factor*frame*expfac;

    tmpar = conv;

  }

// Now multiply by the exponential term and add in the unpiled fraction

//  Real piledflux = outar.sum()*expfac;
  flux *= (1.0 - psffrac);
  flux += nregions*fracexpo*expfac*outar;

// Write out diagnostic information if required. The second number is the
// probability of seeing that number of photons piled up in a frame and the
// third number is the probability of an observed event being due to that
// number of piled up photons.

  Real sum = pfrac.sum();
  if ( sum > 0.0 ) {

    pnorm = expfac;
    for (int ipile=1; ipile<=npiled; ipile++) {
      pnorm *= incounts / ipile;
    } 

  }


}
