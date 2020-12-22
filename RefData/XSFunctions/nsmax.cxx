//-----------------------------------------------------------------------------
//     Model neutron star spectrum for magnetized, partially ionized atmosphere
//     (see Ho, WCG, Potekhin, AY, Chabrier, G 2008, ApJS, submitted)
//     Converted to C++ and to use a single input FITS file by kaa on 5/12/17
//
//     parameter[0] = log(unredshifted effective temperature, in K)
//     parameter[1] = redshift, 1+zg = 1/(1-2GM/Rc^2)^1/2
//     parameter[2] = model to use, see nsmax.dat
//     normalization, (radius/distance)^2, with radius in km
//           and distance in kpc (automatically added by XSPEC)
//-----------------------------------------------------------------------------

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include "BinarySearch.h"
#include <CCfits/CCfits>
#include <memory>


void nsmax(const RealArray& energy, const RealArray& parameter, int spectrum, 
	   RealArray& flux, RealArray& fluxError, const string& init)
{
  using namespace std;
  using namespace CCfits;

  static vector<RealArray> tabulatedFlux;
  static string lastFilename(" ");
  static int lastModel(-1);

  Real temperature(parameter[0]);
  const Real zshift(parameter[1]);
  const int imodel(parameter[2]);


  // check for NSMAX_DIR and set the filename
  string directory = FunctionUtility::getModelString("NSMAX_DIR");
  if (directory.compare("$$NOT$$") == 0) {
    directory = FunctionUtility::modelDataPath();
  }
  const string filename(directory+"nsmax.fits");

  // if the filename has changed or the model number has changed then need
  // to read the input file.

  static int Nmodels;
  static IntegerVector modelNumber;
  static vector<string> modelName;
  static vector<RealArray> tabulatedTemperature;
  static vector<RealArray> tabulatedEnergy;
  static int Ntemps;
  static int Nenergies;

  static RealArray modelTemperature;
  static RealArray modelEnergy;

  if ( filename != lastFilename || imodel != lastModel ) {

    std::unique_ptr<FITS> pInfile((FITS*)0);
    try {
      pInfile.reset(new FITS(filename, Read));
    } catch(...) {
      FunctionUtility::xsWrite("Failed to read"+filename,5);
      return;
    }

    // if the filename has changed then need to reread MODELINFO extension

    if ( filename != lastFilename ) {

      try {
	ExtHDU& Mtable = pInfile->extension("MODELINFO");
	Nmodels = Mtable.rows();
	Mtable.column("MODNUM").read(modelNumber,1,Nmodels);
	Mtable.column("MODNAME").read(modelName,1,Nmodels);
	Mtable.column("TEMPS").readArrays(tabulatedTemperature,1,Nmodels);
	Mtable.column("ENERGIES").readArrays(tabulatedEnergy,1,Nmodels);
      } catch(...) {
	FunctionUtility::xsWrite("Failed to find MODELINFO extension in "+filename,5);
	return;
      }

      // convert tabulatedEnergy array to log10
      for (size_t i=0; i<(size_t)Nmodels; i++) tabulatedEnergy[i] = log10(tabulatedEnergy[i]);
      

    }

    // read the tabulatedFlux information for the current model

    size_t itarg;
    for (itarg=0; itarg<(size_t)Nmodels; itarg++) {
      if ( imodel == modelNumber[itarg] ) break;
    }

    Ntemps = tabulatedTemperature[itarg].size();
    Nenergies = tabulatedEnergy[itarg].size();

    modelTemperature.resize(Ntemps);
    modelTemperature = tabulatedTemperature[itarg];
    modelEnergy.resize(Nenergies);
    modelEnergy = tabulatedEnergy[itarg];

    try {
      ExtHDU& Ftable = pInfile->extension(modelName[itarg]);
      Ftable.column("FLUX").readArrays(tabulatedFlux, 1, Ntemps);
    } catch(...) {
      FunctionUtility::xsWrite("Failed to find "+modelName[itarg]+" extension in "+filename,5);
      return;
    }

    // convert tabulated flux into log10(flux)

    for (size_t i=0; i<(size_t)Ntemps; i++) tabulatedFlux[i] = log10(tabulatedFlux[i]);

    lastFilename = filename;
    lastModel = imodel;
  }

  Real minT = modelTemperature[0];
  Real maxT = modelTemperature[Ntemps-1];

  // find the bracketing temperatures and interpolate the tabulatedFlux array

  int iT = Numerics::BinarySearch(modelTemperature,temperature);
  Real checkFactor(1.0);
  if(iT < 0) {
    cout<<"!! Warning: T="<<temperature<<" is outside the limits for which model atmospheres are calculated ("<< minT << "..." << maxT << ")"<<std::endl;
    cout<<"The results are not trustable for this temperature!!"<<std::endl;
    if ( iT == -1 ) {
      checkFactor = 1.+1.e4*(minT-temperature)/temperature;
      temperature = minT;
      iT = 0;
    } else if(iT == -2) {
      checkFactor = 1.+1.e4*(temperature-maxT)/maxT;
      temperature = maxT;
      iT = Ntemps-2;
    }
  }
  Real lowT = modelTemperature[iT];
  Real highT= modelTemperature[iT+1];

  RealArray fluxInterp(Nenergies);
  fluxInterp = ( (highT-temperature)*tabulatedFlux[iT] + 
		 (temperature-lowT)*tabulatedFlux[iT+1] )/(highT-lowT);

  // interpolate onto the output flux array - notice using the upper energy of each
  // bin in the interpolation (this is not really correct but is what was supplied
  // for this model).

  flux.resize(energy.size()-1);
  for (size_t ie=0; ie<energy.size()-1; ie++) {
    Real zEnergy = log10(energy[ie+1] * zshift);
    int ie1 = Numerics::BinarySearch(modelEnergy,zEnergy);
    Real ce1, ce2;
    if ( ie1 == -1 ) {
      ie1 = 0;
      ce1 = 1.0;
      ce2 = 0.0;
    } else if ( ie1 == -2 ) {
      ie1 = Nenergies - 2;
      ce1 = 0.0;
      ce2 = 1.0;
    } else {
      ce1 = (modelEnergy[ie1+1] - zEnergy) / (modelEnergy[ie1+1] - modelEnergy[ie1]);
      ce2 = 1.0 - ce1;
    }
    flux[ie] = ( fluxInterp[ie1]*ce1 + fluxInterp[ie1+1]*ce2 );

    // Scale by redshift
    flux[ie] -= log10(zshift);
    // Scale by (radius/distance)^2, where radius (in km) and distance (in kpc)
    // PC=3.0856775807e18 cm
    flux[ie] -= 32.978701090;
    // Convert from ergs/(s cm^2 Hz) to counts/(s cm^2 keV)
    // HH=6.6260693e-27 ergs s
    flux[ie] += 26.178744 - log10(energy[ie+1]);
    flux[ie] = pow(10.0, flux[ie]);
    // Convert counts/(s cm^2 bin)
    flux[ie] *= (energy[ie+1] - energy[ie]);

  }


  return;
}
