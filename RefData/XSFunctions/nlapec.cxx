// XSPEC model for an APEC spectrum with no lines
// parameters :
//       0:    kT (keV)
//       1:    Abundances (H & He assumed Solar)
//       2:    Redshift
// Norm = (4 pi 1e14)^-1 Int n^2 dV / D^2

#include <xsTypes.h>
#include <functionMap.h>
#include "Aped.h"
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <sstream>

void nlapec(const RealArray& energyArray, const RealArray& params,
	    int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	    const string& initString)
{

  bool qtherm = false;
  bool noLines = true;
  Real velocity = 0.0;
  Real dem = 1.0;

  IntegerVector Zinput(30);
  RealArray abundance(30);
  for (size_t i=0; i<30; i++) Zinput[i] = i+1;
  abundance[0] = 1.0;
  abundance[1] = 1.0;
  for (size_t i=2; i<30; i++) abundance[i] = params[1];

  int status = calcCEISpectrum(energyArray, Zinput, abundance, params[2], params[0], dem,
			       qtherm, velocity, noLines, flux, fluxErr);

  if ( status != 0 ) {
    std::ostringstream message;
    message << "***nlapec : failure in calcCEISpectrum, status = " << status;
    FunctionUtility::xsWrite(message.str(), 5);
  }

  return;

}
