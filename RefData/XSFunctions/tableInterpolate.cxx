// routine to return the fluxes and errors from the table model file fileName using
// the parameter values in params. This routine handles additive, multiplicative and
// exponential multiplicative models.

// This uses classes defined in the heasp library
#include "table.h"
#include "SPutils.h"

#include <XSFunctions/Utilities/FunctionUtility.h>
#include <Numerics/LinearInterp.h>

void tableInterpolate(const RealArray& energyArray, 
		      const RealArray& params, 
		      string fileName,
		      int spectrumNumber,
		      RealArray& fluxArray, 
		      RealArray& fluxErrArray,
		      const string& initString,
		      const string& tableType,
		      const bool readFull)
{
  using namespace std;
  using namespace Numerics;
  using namespace Rebin;

  if ( tableType != "add" && tableType != "mul" && tableType != "exp" ) {
    FunctionUtility::xsWrite(tableType+" is not a valid table type", 5);
    return;
  }

  static vector<pair<string,table> > saveTables;

  // check whether we already have the table from the input filename

  bool found = false;
  size_t itable;
  for (size_t i=0; i<saveTables.size(); i++) {
    if ( fileName == saveTables[i].first ) {
      found = true;
      itable = i;
      break;
    }
  }

  if ( !found ) {
    table readTable;
    int status = readTable.read(fileName, readFull);
    if ( status != 0 ) {
      FunctionUtility::xsWrite("Failed to read "+fileName, 5);
      return;
    }
    pair<string,table> savePair;
    savePair = make_pair(fileName, readTable);
    saveTables.push_back(savePair);
    itable = saveTables.size()-1;
  }

  table& inTable = saveTables[itable].second;

  size_t nE(energyArray.size());
  RealArray tableEnergyBins, tableValues, tableErrors;
  Real minEnergy = energyArray[0];
  Real maxEnergy = energyArray[nE-1];

  // interpolate on the table model grid

  int status = inTable.getValues(params, minEnergy, maxEnergy, tableEnergyBins,
				 tableValues, tableErrors);
  if ( status != 0 ) {
    FunctionUtility::xsWrite(SPgetErrorStack(), 5);
    SPclearErrorStack();
    return;
  }

  // get the values to be used for energies below or above those tabulated in the
  // table

  Real loLimit = inTable.getLowEnergyLimit();
  Real hiLimit = inTable.getHighEnergyLimit();
  
  // now remap onto energyArray. for the additive model we rebin, for the others
  // we interpolate

  size_t inputBin;
  size_t outputBin;

  IntegerVector startBin(nE-1), endBin(nE-1);
  RealArray startWeight(nE-1), endWeight(nE-1);

  findFirstBins(tableEnergyBins, energyArray, FUZZY, inputBin, outputBin);
  initializeBins(tableEnergyBins, energyArray, FUZZY, inputBin, outputBin,
		 startBin, endBin, startWeight, endWeight);

  fluxArray.resize(nE-1);
  if ( tableType == "add" ) {
    rebin(tableValues, startBin, endBin, startWeight, endWeight, fluxArray, loLimit,
	  hiLimit);
  } else if ( tableType == "mul" ) {
    interpolate(tableValues, startBin, endBin, startWeight, endWeight, fluxArray, false, loLimit, hiLimit);
  } else {
    interpolate(tableValues, startBin, endBin, startWeight, endWeight, fluxArray, true, loLimit, hiLimit);
  }

  if ( tableErrors.size() > 0 ) {
    fluxErrArray.resize(energyArray.size()-1);
    if ( tableType == "add" ) {
      rebin(tableErrors, startBin, endBin, startWeight, endWeight, fluxErrArray);
    } else if ( tableType == "mul" ) {
      interpolate(tableErrors, startBin, endBin, startWeight, endWeight, fluxErrArray, false);
    } else {
      interpolate(tableErrors, startBin, endBin, startWeight, endWeight, fluxErrArray, true);
      fluxErrArray *= fluxArray;
    }
  }

  return;
}


// C routine interface to tableInterpolate

extern "C" void tabint(float* ear, int ne, float* param, int npar, const char* filenm, int ifl, 
            const char* tabtyp, float* photar, float* photer)
{
   string fileName(filenm);
   string initString("");
   string tableType(tabtyp);

   RealArray energyArray(.0,ne+1);
   RealArray parameters(npar);
   RealArray fluxArray, fluxErrArray;

   for (size_t i=0; i<(size_t)ne+1; i++) energyArray[i] = ear[i];
   for (size_t i=0; i<(size_t)npar; i++) parameters[i] = param[i];

   tableInterpolate(energyArray, parameters, fileName, ifl, fluxArray, 
		    fluxErrArray, initString, tableType, true);

   for (size_t i=0; i<(size_t)ne; i++) photar[i] = fluxArray[i];
   if ( fluxErrArray.size() > 0 ) {
     for (size_t i=0; i<(size_t)ne; i++) photer[i] = fluxErrArray[i];
   } else {
     for (size_t i=0; i<(size_t)ne; i++) photer[i] = 0.0;
   }

   return;
}

// Fortran wrapper

#include <cfortran.h>

FCALLSCSUB9(tabint,TABINT,tabint,FLOATV,INT,FLOATV,INT,STRING,INT,STRING,FLOATV,FLOATV)

