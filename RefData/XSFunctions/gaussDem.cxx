#include <functionMap.h>
#include <xsTypes.h>

// prototype for call to vgadem model in vgaussDem.cxx.

void vgaussDem(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString);

// XSPEC model subroutine to calculate collisional plasma with a gaussian DEM
// 
// Parameters:
//    param(1) = Temperature mean
//    param(2) = Temperature sigma
//    param(3) = nH (cm^-3)  Fixed at 1 for most applications
//    param(4) = abundance
//    param(5) = redshift
//    param(6) = switch(0=calculate MEKAL model, 1=interpolate MEKAL model,
//                       2=AtomDB model)


void gaussDem(const RealArray& energyArray, const RealArray& params,
	      int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	      const string& initString)
{

  RealArray pparams(19);

  for (int i=0; i<3; i++) pparams[i] = params[i];
  pparams[3] = 1.0;
  for (int i=4; i<17; i++) pparams[i] = params[3];
  pparams[17] = params[4];
  pparams[18] = params[5];

  vgaussDem(energyArray, pparams, spectrumNumber, flux, fluxErr, 
	    initString);

  return;

}

