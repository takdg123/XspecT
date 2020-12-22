// C++ conversion of xsdili.f for the diskline model
// model to calculate the line shape for a rotating accretion
// disk. does not include GR effects. note that if param[1] is
// set to 10 then do the special case of an accretion disk
// emissivity.
// Formulae from the Appendix to Fabian et al. 1989, MNRAS 238, 729
// Parameters :
//    0        line energy
//    1        power law index for emissivity (10 for disk)
//    2        inner radius (GM/c**2)
//    3        outer radius (GM/c**2)
//    4        inclination  (degrees)

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Numerics/BinarySearch.h>
#include <iostream>

void diskline(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	      const string& initString)
{
  const int numberRadii(1000);
  const int numberAngles(90);
  const Real dera(57.295779);

  //  this model does not calculate errors
  const size_t nE(energyArray.size());
  fluxArray.resize(nE-1,0.0);
  fluxErrArray.resize(0);

  // convert input parameters - inclination set to radians.

  Real enel = params[0];
  Real alp = params[1];
  Real ri = params[2];
  Real rilog10 = log10(ri);
  Real ro = params[3];

  Real sinIncl = sin(params[4]/dera);
  Real sinIncl2 = sinIncl*sinIncl;
  Real cosIncl = cos(params[4]/dera);
  Real cosIncl2 = cosIncl*cosIncl;

  // trap case where inner radius is greater than outer

  if ( ri >= ro ) {
    FunctionUtility::xsWrite("Inner radius > outer radius  -  model is invalid", 5);
    return;
  }

  // trap case of inner radius being less than the last stable orbit.

  if (ri < 6.0) {
    FunctionUtility::xsWrite("Inner radius < 6 R_g  -  model is invalid", 5);
    return;
  }

  // calculate log step in radius

  Real dlogr = (log10(ro)-rilog10)/((Real)(numberRadii-1));

  // calculate radii

  RealArray radius(numberRadii);
  for (size_t i=0; i<radius.size(); i++) radius[i] = pow(10.0,rilog10+i*dlogr);

  // set up azimuthal angle quantities so we don't have to recalculate them
  // for every radius. phi is the azimuthal disk coordinate with phi=0 being the
  // line of nodes where the disk crosses the plane of the sky.
  
  RealArray phi(2*numberAngles+1);
  for (size_t i=0; i<2*numberAngles+1; i++) phi[i] = i * (90.0/numberAngles) / dera;

  // tanxi2 is defined by equation A7 and cosbeta by equation A3 in Fabian et al.
  // xi + pi/2 is the angle between direction of emission of the photon and the
  // line connecting the emitting point to the BH. beta is the angle between the
  // disk plane and the plane defined by the BH, the emitting point, and the observer.
  
  RealArray tanxi2(2*numberAngles+1), cosbeta(2*numberAngles+1);
  tanxi2 = sinIncl2*sin(phi)*sin(phi) / (1.0 - sinIncl2*sin(phi)*sin(phi));
  cosbeta = sinIncl*cos(phi) / sqrt(cosIncl2 + sinIncl2*cos(phi)*cos(phi));

  // its also useful to have an array of energy bin sizes

  RealArray energyArrayBinSize(nE-1);
  for (size_t i=0; i<nE-1; i++) energyArrayBinSize[i] = energyArray[i+1]-energyArray[i];

  // big loop for radii. note that the radius cannot reach 3 else
  // the metric goes singular

  for (size_t irad=0; irad<radius.size()-1; irad++) {

    Real ra = (radius[irad]+radius[irad+1])/2.;
    Real dra = radius[irad+1] - radius[irad];

    Real rafact = sqrt(1.-3./ra);

    // if power-law index is less than ten use to calculate emissivity
    // else use the accretion disk emissivity law.

    Real fra;
    if ( alp < 9.9 ) {
      fra = pow(ra,alp);
    } else {
      fra = (1.-sqrt(6./ra))/(ra*ra*ra);
    }

    // calculate the radius dependent part of the flux density
    
    Real fluxfactor = fra*ra*dra*2./dera;

    // loop over disk azimuthal angles. this is in steps of 2 degrees from
    // 0 (line of nodes) to 178.
    
    for (size_t kang = 0; kang<numberAngles; kang++) {

      // calculate mean redshift (1+z = zpo) for the bin using equation A2
      // in Fabian et al.

      Real zpo = (1.+cosbeta[2*kang+1]/sqrt(ra*(1.+tanxi2[2*kang+1])-2.))/rafact;

      // and the low and high redshifts for the bin. note the traps for
      // the case of an inclination of 90 degrees.

      Real zpol;
      if ( params[4] > 89.9 && kang == (numberAngles/2) ) {
	zpol = 1./rafact;
      } else {
	zpol = (1.+cosbeta[2*kang]/sqrt(ra*(1.+tanxi2[2*kang])-2.))/rafact;
      }

      Real zpoh;
      if ( params[4] >  89.9 && kang == (numberAngles/2 - 1) ) {
	zpoh = 1./rafact;
      } else {
	zpoh = (1.+cosbeta[2*kang+2]/sqrt(ra*(1.+tanxi2[2*kang+2])-2.))/rafact;
      }

      // enobl and enobh are the lower and upper observed energy from
      // this azimuthal and radial bin

      Real enobl = fmin(enel/zpol, enel/zpoh);
      Real enobh = fmax(enel/zpol, enel/zpoh);
      Real total = enobh - enobl;

      // calculate flux density from this bin using equation A4 in Fabian et al.

      Real fluxDensity(0.0);
      if ( total > 0.0 ) fluxDensity = fluxfactor/(zpo*zpo*zpo)/total;

      // find fractions of emission from this bin to place in each energy
      // range. Find the index in energyArray immediately below enobl and
      // enobh

      int ienobl = Numerics::BinarySearch(energyArray, enobl);
      int ienobh = Numerics::BinarySearch(energyArray, enobh);

      if (ienobl == -2 || ienobh == -1 || total <= 0.0) {
	// do nothing because (enobl,enobh) lies outside energyArray
	// or enobl=enobh so there is no contribution
      } else if ( ienobl == -1 && ienobh == -2 ) {
	// (enobl,enobh) encompasses the entire energyArray
	for (size_t ie=0; ie<nE-1; ie++) {
	  fluxArray[ie] += fluxDensity*energyArrayBinSize[ie];
	}
      } else if ( ienobl == -1 && ienobh > 0 ) {
	// enobl below energyArray but ienobh within it
	for (size_t ie=0; ie<(size_t)ienobh; ie++) {
	  fluxArray[ie] += fluxDensity*energyArrayBinSize[ie];
	}
	fluxArray[ienobh] += fluxDensity*(enobh-energyArray[ienobh]);	     
      } else if ( ienobl > 0 && ienobh == -2 ) {
	// enobl within energyArray but ienobh above it
	fluxArray[ienobl] += fluxDensity*(energyArray[ienobl+1]-enobl);
	for (size_t ie=ienobl+1; ie<nE-1; ie++) {
	  fluxArray[ie] += fluxDensity*energyArrayBinSize[ie];
	}
      } else {
	// both enobl and enobh within energyArray
	if ( ienobl == ienobh ) {
	  fluxArray[ienobl] += fluxDensity*total;
	} else {
	  fluxArray[ienobl] += fluxDensity*(energyArray[ienobl+1]-enobl);
	  for (size_t ie=(size_t)ienobl+1; ie<(size_t)ienobh; ie++) {
	    fluxArray[ie] += fluxDensity*energyArrayBinSize[ie];
	  }
	  fluxArray[ienobh] += fluxDensity*(enobh-energyArray[ienobh]);
	}
      }

      // end loop over angles	
    }
    // end loop over radii
  }

  // normalise values to total

  Real spm = 0.;
  for (size_t ie=0; ie<nE-1; ie++) spm += fluxArray[ie];
  if ( spm != 0.0 ) {
    for (size_t ie=0; ie<nE-1; ie++) fluxArray[ie] /= spm;;
  }    

  return;
}
