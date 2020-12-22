//-----------------------------------------------------------------------------
//    Model neutron star spectrum for magnetized, partially ionized atmosphere
//    (see Ho, WCG, Potekhin, AY, Chabrier, G 2008, ApJS, 178, 102)
//     Converted to C++ and to use a single input FITS file by kaa on 5/16/17
//
//    parameter[0] = log(unredshifted effective temperature, in K)
//    parameter[1] = neutron star mass, in solar masses
//    parameter[2] = neutron star radius, in km
//    parameter[3] = distance to neutron star, in kpc
//    parameter[4] = model to use, see MODELINFO extension of nsmaxg.fits
//-----------------------------------------------------------------------------

#include "xsTypes.h"
#include <XSFunctions/functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include "BinarySearch.h"
#include <CCfits/CCfits>
#include <memory>

// routine to do the interpolation over temperature and surface gravity
// returns tableFlux.
void calcHoNSatmos(Real& temperature, Real& gg, const string filename,
		   const int imodel, RealArray& tableEnergy, 
		   RealArray& tableFlux);


void nsmaxg(const RealArray& energy, const RealArray& parameter, int spectrum, 
	    RealArray& flux, RealArray& fluxError, const string& init)
{
  using namespace std;

  Real temperature(parameter[0]);
  const Real redshift(1.0/sqrt(1.0-2.95316*(parameter[1]/parameter[2])));
  Real gg(log10(1.3271e16*((parameter[1]/(parameter[2]*parameter[2]))*redshift)));
  const int imodel(parameter[4]);

  // check for NSMAXG_DIR and set the filename
  string directory = FunctionUtility::getModelString("NSMAXG_DIR");
  if (directory.compare("$$NOT$$") == 0) {
    directory = FunctionUtility::modelDataPath();
  }
  const string filename(directory+"nsmaxg.fits");

  RealArray tableEnergy, tableFlux;
  calcHoNSatmos(temperature, gg, filename, imodel, tableEnergy, tableFlux);

  // interpolate onto the output flux array - notice using the upper energy of each
  // bin in the interpolation (this is not really correct but is what was supplied
  // for this model). Using log interpolation so convert tableEnergy array to log10

  tableEnergy = log10(tableEnergy);

  flux.resize(energy.size()-1);
  for (size_t ie=0; ie<energy.size()-1; ie++) {
    Real zEnergy = log10(energy[ie+1] * redshift);
    int ie1 = Numerics::BinarySearch(tableEnergy,zEnergy);
    Real ce1, ce2;
    if ( ie1 == -1 ) {
      ie1 = 0;
      ce1 = 1.0;
      ce2 = 0.0;
    } else if ( ie1 == -2 ) {

      ie1 = tableEnergy.size() - 2;
      ce1 = 0.0;
      ce2 = 1.0;
    } else {
      ce1 = (tableEnergy[ie1+1] - zEnergy) / (tableEnergy[ie1+1] - tableEnergy[ie1]);
      ce2 = 1.0 - ce1;
    }
    flux[ie] = ( tableFlux[ie1]*ce1 + tableFlux[ie1+1]*ce2 );

    // Scale by redshift
    flux[ie] -= log10(redshift);
    // Scale by (radius/distance)^2, where radius (in km) and distance (in kpc)
    // PC=3.0856775807e18 cm
    flux[ie] += -32.978701090 + 2.0*log10(parameter[2]/parameter[3]);
    // Convert from ergs/(s cm^2 Hz) to counts/(s cm^2 keV)
    // HH=6.6260693e-27 ergs s
    flux[ie] += 26.178744 - log10(energy[ie+1]);
    flux[ie] = pow(10.0, flux[ie]);
    // Convert counts/(s cm^2 bin)
    flux[ie] *= (energy[ie+1] - energy[ie]);

  }

  return;
}

void calcHoNSatmos(Real& temperature, Real& gg, const string filename,
		   const int imodel, RealArray& tableEnergy, 
		   RealArray& tableFlux)
{
  using namespace std;
  using namespace CCfits;

  static vector<vector<RealArray> > tabulatedFlux;
  static string lastFilename(" ");
  static int lastModel(-1);


  // if the filename has changed or the model number has changed then need
  // to read the input file.

  static int Nmodels;
  static IntegerVector modelNumber;
  static vector<string> modelName;
  static vector<RealArray> tabulatedTemperature;
  static vector<RealArray> tabulatedGrav;
  static vector<RealArray> tabulatedEnergy;
  static int Ntemps;
  static int Ngravs;
  static int Nenergies;

  static RealArray modelTemperature;
  static RealArray modelGrav;
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
	Mtable.column("GRAVS").readArrays(tabulatedGrav,1,Nmodels);
	Mtable.column("ENERGIES").readArrays(tabulatedEnergy,1,Nmodels);
      } catch(...) {
	FunctionUtility::xsWrite("Failed to find MODELINFO extension in "+filename,5);
	return;
      }

    }

    // find the index for the requested model and create the arrays for this model

    size_t itarg;
    for (itarg=0; itarg<(size_t)Nmodels; itarg++) {
      if ( imodel == modelNumber[itarg] ) break;
    }

    Ntemps = tabulatedTemperature[itarg].size();
    Ngravs = tabulatedGrav[itarg].size();
    Nenergies = tabulatedEnergy[itarg].size();

    modelTemperature.resize(Ntemps);
    modelTemperature = tabulatedTemperature[itarg];
    modelGrav.resize(Ngravs);
    modelGrav = tabulatedGrav[itarg];
    modelEnergy.resize(Nenergies);
    modelEnergy = tabulatedEnergy[itarg];

    // read the tabulatedFlux information for the current model

    vector<RealArray> readFlux;
    try {
      ExtHDU& Ftable = pInfile->extension(modelName[itarg]);
      Ftable.column("FLUX").readArrays(readFlux, 1, Ntemps*Ngravs);
    } catch(...) {
      FunctionUtility::xsWrite("Failed to find "+modelName[itarg]+" extension in "+filename,5);
      return;
    }

    // map readFlux onto tabulatedFlux to separate out the looping over temperature
    // and g the convert to log10(flux)

    tabulatedFlux.resize(Ntemps);
    size_t ipt=0;
    for (size_t itemp=0; itemp<(size_t)Ntemps; itemp++) {
      tabulatedFlux[itemp].resize(Ngravs);
      for (size_t igrav=0; igrav<(size_t)Ngravs; igrav++) {
	tabulatedFlux[itemp][igrav].resize(Nenergies);
	tabulatedFlux[itemp][igrav] = log10(readFlux[ipt]);
	ipt++;
      }
    }

    lastFilename = filename;
    lastModel = imodel;
  }

  Real minT = modelTemperature[0];
  Real maxT = modelTemperature[Ntemps-1];
  Real ming = modelGrav[0];
  Real maxg = modelGrav[Ngravs-1];

  // find the bracketing temperatures and gravs and interpolate the tabulatedFlux array

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
  Real ct1 = (modelTemperature[iT+1]-temperature)/(modelTemperature[iT+1]-modelTemperature[iT]);
  Real ct2 = 1.0 - ct1;

  Real cg1(1.0), cg2(0.0);
  int ig = 0;
  if ( Ngravs > 1 ) {
    ig = Numerics::BinarySearch(modelGrav,gg);
    if(ig < 0) {
      cout<<"!! Warning: g="<<gg<<" is outside the limits for which model atmospheres are calculated ("<< ming << "..." << maxg << ")"<<std::endl;
      cout<<"The results are not trustable for this g!!"<<std::endl;
      if ( ig== -1 ) {
	checkFactor += 1.e4*(ming-gg)/gg;
	gg = ming;
	ig = 0;
      } else if (ig == -2) {
	checkFactor += 1.e4*(gg-maxg)/maxg;
	gg = maxg;
	ig = Ngravs-2;
      }
    }
    cg1 = (modelGrav[ig+1]-gg)/(modelGrav[ig+1]-modelGrav[ig]);
    cg2 = 1.0 - cg1;
  }

  // interpolate on the temperature and g values

  tableFlux.resize(Nenergies);
  tableFlux = cg1 * ( ct1*tabulatedFlux[iT][ig] + ct2*tabulatedFlux[iT+1][ig] );
  if ( Ngravs > 1 ) {
    tableFlux += cg2 * ( ct1*tabulatedFlux[iT][ig+1] + ct2*tabulatedFlux[iT+1][ig+1] );
  }

  // place modelEnergy contents in tableEnergy output array

  tableEnergy.resize(modelEnergy.size());
  tableEnergy = modelEnergy;

  return;
}
