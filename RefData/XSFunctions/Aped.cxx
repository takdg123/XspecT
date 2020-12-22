// Code for Aped, ApedTemperatureRecord, and ApedElementRecord classes


#include <XSFunctions/IonBalNei.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Numerics/BinarySearch.h>
#include <XSUtil/Numerics/Numerics.h>
#include <sstream>
#include <fstream>
#include <cmath>
#include <CCfits/CCfits>
#include <memory>

#include <time.h>

#include "Aped.h"

using namespace CCfits;
using namespace Numerics;

// conversion from keV to Kelvin
#define keVtoK 1.1604505e7
// km in cm
#define km_cm 1e5



// declaration for the calcManyLines function in calcLines.cxx.

void calcManyLines(const RealArray& energyArray, const RealArray& ecenterArray,
		   const std::vector<RealArray>& lineParamsArray,
		   const RealArray& linefluxArray, const Real crtLevel,
		   const int lineShape, const bool qspeedy, RealArray& fluxArray);

// Methods for Aped class. Most of code is to load from the cont and line files

// default constructor

Aped::Aped()
{
  TemperatureRecord.resize(0);
  Temperatures.resize(0);
  NelementsLine.resize(0);
  NelementsCoco.resize(0);
  Nline.resize(0);
  Ncont.resize(0);
  coconame = "";
  linename = "";
  noLines = false;
  noResonanceLines = false;
  thermalBroadening = false;
  velocityBroadening = 0.0;
  logTempInterpolation = false;
}

// destructor

Aped::~Aped()
{
  TemperatureRecord.resize(0);
  Temperatures.resize(0);
  NelementsLine.resize(0);
  NelementsCoco.resize(0);
  Nline.resize(0);
  Ncont.resize(0);
}

// Set the noLines, noResonanceLines, thermalBroadening, velocityBroadening, and
// logTempInterpolation by checking the APECNOLINES, APECNORES, APECTHERMAL,  
// APECVELOCITY, and APECLOGTINTERP values

void Aped::SetNoLines(const bool qno)
{
  string pname = "APECNOLINES";
  string pvalue = FunctionUtility::getModelString(pname);
  noLines = qno;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    pvalue = XSutility::lowerCase(pvalue);
    if (pvalue.substr(0,1) == "t" || pvalue.substr(0,1) == "y" || 
	pvalue.substr(0,2) == "on") {
      noLines = true;
    } else {
      noLines = false;
    }
  }
  return;
}

void Aped::SetNoResonanceLines(const bool qnores)
{
  string pname = "APECNORES";
  string pvalue = FunctionUtility::getModelString(pname);
  noResonanceLines = qnores;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    pvalue = XSutility::lowerCase(pvalue);
    if (pvalue.substr(0,1) == "t" || pvalue.substr(0,1) == "y" || 
	pvalue.substr(0,2) == "on") {
      noResonanceLines = true;
    } else {
      noResonanceLines = false;
    }
  }
  return;
}

void Aped::SetThermalBroadening(const bool qtherm)
{
  string pname = "APECTHERMAL";
  string pvalue(FunctionUtility::getModelString(pname));
  thermalBroadening = qtherm;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    pvalue = XSutility::lowerCase(pvalue);
    if (pvalue.substr(0,1) == "t" || pvalue.substr(0,1) == "y" || 
	pvalue.substr(0,2) == "on") {
      thermalBroadening = true;
    } else {
      thermalBroadening = false;
    }
  }
  return;
}

void Aped::SetVelocityBroadening(const Real velocity)
{
  string pname = "APECVELOCITY";
  string pvalue = FunctionUtility::getModelString(pname);
  velocityBroadening = velocity;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    std::istringstream iss(pvalue);
    Real tmpVel = 0.0;
    if (!(iss >> tmpVel) || !iss.eof()) {
      std::ostringstream oss;
      oss << "Failed to read value from APECVELOCITY - assuming "
	  << velocity << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
    } else {
      if (velocity > 0. && fabs(tmpVel-velocity) > 1e-5) {
	string msg = "Warning: inconsistency between velocity broadening set by model and APECVELOCITY variable - using that from model.";
	FunctionUtility::xsWrite(msg, 10);
      } else
	velocityBroadening = tmpVel;
    }
  }
  return;
}

void Aped::SetMinimumLinefluxForBroadening(const Real flux)
{
  string pname = "APECMINFLUX";
  string pvalue = FunctionUtility::getModelString(pname);
  minimumLinefluxForBroadening = flux;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    std::istringstream iss(pvalue);
    Real tmpFlux = 0.0;
    if (!(iss >> tmpFlux) || !iss.eof()) {
      std::ostringstream oss;
      oss << "Failed to read value from APECMINFLUX - assuming "
	  << flux << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
    } else {
      minimumLinefluxForBroadening = tmpFlux;
    }
  }
  return;
}

void Aped::SetBroadenPseudoContinuum(const bool qpbroad)
{
  string pname = "APECBROADPSEUDO";
  string pvalue = FunctionUtility::getModelString(pname);
  broadenPseudoContinuum = qpbroad;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    pvalue = XSutility::lowerCase(pvalue);
    if (pvalue.substr(0,1) == "t" || pvalue.substr(0,1) == "y" || 
	pvalue.substr(0,2) == "on") {
      broadenPseudoContinuum = true;
    } else {
      broadenPseudoContinuum = false;
    }
  }
  return;
}

void Aped::SetLogTempInterpolation(const bool qltinterp)
{
  string pname = "APECLOGTINTERP";
  string pvalue = FunctionUtility::getModelString(pname);
  logTempInterpolation = qltinterp;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    pvalue = XSutility::lowerCase(pvalue);
    if (pvalue.substr(0,1) == "t" || pvalue.substr(0,1) == "y" || 
	pvalue.substr(0,2) == "on") {
      logTempInterpolation = true;
    } else {
      logTempInterpolation = false;
    }
  }
  return;
}

// Wrap-up routine to check the APECNOLINES, APECNORES, APECTHERMAL, APECVELOCITY, 
// APECMINFLUX, APECBROADPSEUDO, and APECLOGTINTERP values and set the object

void Aped::SetLineSpecs(const bool qno, const bool qnores, const bool qtherm, 
			const Real velocity, const Real flux, const bool qpbroad,
			const bool qltinterp)
{

  SetNoLines(qno);
  SetNoResonanceLines(qnores);
  SetThermalBroadening(qtherm);
  SetVelocityBroadening(velocity);
  SetMinimumLinefluxForBroadening(flux);
  SetBroadenPseudoContinuum(qpbroad);
  SetLogTempInterpolation(qltinterp);

    // write out diagnostic information at chatter level 20

  std::ostringstream oss;
  if ( noLines ) {
    oss << "No lines will be calculated" << std::endl;
  } else {
    if ( noResonanceLines ) {
      oss << "No resonance lines will be calculated" << std::endl;
    }
    oss << "Thermal broadening : ";
    thermalBroadening ? oss << "T" : oss << "F";
    oss << std::endl;
    oss << "Velocity broadening : " << velocityBroadening << std::endl;;
    if ( thermalBroadening || velocityBroadening > 0.0 ) {
      if ( minimumLinefluxForBroadening > 0.0 ) {
	oss << "Thermal and/or velocity broadening will be applied only for lines with flux > " << minimumLinefluxForBroadening << std::endl;
      }
      if ( broadenPseudoContinuum ) {
	oss << "Thermal and/or velocity broadening will be applied to the pseudo-continuum" << std::endl;
      }
      if ( logTempInterpolation ) {
	oss << "Logarithmic interpolation on tabulated temperatures" << std::endl;
      } else {
	oss << "Linear interpolation on tabulated temperatures" << std::endl;
      }
    }
  }
  FunctionUtility::xsWrite(oss.str(), 20);
  return;
}

// read from files

int Aped::Read(string cocofilename, string linefilename)
{

  const string hduName("PARAMETERS");

  const vector<string> hduKeys;
  const vector<string> primaryKeys;

  // open files
  std::unique_ptr<FITS> pCocofile((FITS*)0);
  std::unique_ptr<FITS> pLinefile((FITS*)0);

  string msg = "Reading continuum data from " + cocofilename;
  FunctionUtility::xsWrite(msg,25);
  msg = "Reading line data from " + linefilename;
  FunctionUtility::xsWrite(msg,25);

  try {
    pCocofile.reset(new FITS(cocofilename,CCfits::Read,hduName,false,hduKeys,primaryKeys,(int)1));
  } catch(...) {
    msg = "Failed to read continuum data from " + cocofilename;
    FunctionUtility::xsWrite(msg,5);
    return(1);
  }
  try {
    pLinefile.reset(new FITS(linefilename,CCfits::Read,hduName,false,hduKeys,primaryKeys,(int)1));
  } catch(...) {
    msg = "Failed to read line data from " + linefilename;
    FunctionUtility::xsWrite(msg,5);
    return(2);
  }

  // store the filenames we are using

  coconame = cocofilename;
  linename = linefilename;

  // Reset the temperature records in case they have been set from a different
  // pair of files

  TemperatureRecord.clear();  

  ExtHDU& cocoExt = pCocofile->currentExtension();
  ExtHDU& lineExt = pLinefile->currentExtension();

  // get the number of temperatures from the number of rows in the PARAMETERS extension

  size_t NumberTemperatures = (size_t)cocoExt.rows();

  // Load the temperatures from the PARAMETERS extension

  Temperatures.resize(NumberTemperatures);

  try {
    cocoExt.column("kT").read(Temperatures,1,NumberTemperatures);
  } catch(...) {
    msg = "Failed to read kT column in PARAMETERS extension of " + cocofilename;
    FunctionUtility::xsWrite(msg,5);
    return(3);
  }

  // check for consistency between the two files

  size_t LineNumberTemperatures = (size_t)lineExt.rows();

  if ( NumberTemperatures != LineNumberTemperatures ) return(3);

  RealArray LineTemperatures(LineNumberTemperatures);
  try {
    lineExt.column("kT").read(LineTemperatures,1,LineNumberTemperatures);
  } catch(...) {
    msg = "Failed to read kT column in PARAMETERS extension of " + linefilename;
    FunctionUtility::xsWrite(msg,5);
    return(4);
  }

  // loop round temperatures
  
  for (size_t iTemp=0; iTemp<NumberTemperatures; iTemp++) {
    if ( Temperatures[iTemp] != LineTemperatures[iTemp] ) {
      msg = "Inconsistency in tabulated temperatures between " + 
	cocofilename + " and " + linefilename;
      FunctionUtility::xsWrite(msg,5);
      return(5);
    }
  }

  // store other information from the PARAMETERS extension

  NelementsLine.resize(NumberTemperatures);
  NelementsCoco.resize(NumberTemperatures);
  Nline.resize(NumberTemperatures);
  Ncont.resize(NumberTemperatures);

  try {
    lineExt.column("Nelement").read(NelementsLine,1,NumberTemperatures);
  } catch(...) {
    msg = "Failed to read Nelement column in PARAMETERS extension of " + linefilename;
    FunctionUtility::xsWrite(msg,5);
    return(6);
  }
  try {
    cocoExt.column("NElement").read(NelementsCoco,1,NumberTemperatures);
  } catch(...) {
    msg = "Failed to read NElement column in PARAMETERS extension of " + cocofilename;
    FunctionUtility::xsWrite(msg,5);
    return(7);
  }
  try {
    lineExt.column("Nline").read(Nline,1,NumberTemperatures);
  } catch(...) {
    msg = "Failed to read Nline column in PARAMETERS extension of " + linefilename;
    FunctionUtility::xsWrite(msg,5);
    return(7);
  }
  try {
    cocoExt.column("NCont").read(Ncont,1,NumberTemperatures);
  } catch(...) {
    msg = "Failed to read NCont column in PARAMETERS extension of " + cocofilename;
    FunctionUtility::xsWrite(msg,5);
    return(8);
  }

  // size the TemperatureRecord array

  TemperatureRecord.resize(NumberTemperatures);

  return(0);

}

// read in a temperature record

int Aped::ReadTemperature(const int TemperatureIndex)
{
  vector<int> TemperatureIndexArray(1);
  TemperatureIndexArray[0] = TemperatureIndex;
  return ReadTemperature(TemperatureIndexArray);
}

// read in temperature records

int Aped::ReadTemperature(const vector<int>& TemperatureIndex)
{
  const string hduName("PARAMETERS");

  const vector<string> hduKeys;
  const vector<string> primaryKeys;

  // open files
  std::unique_ptr<FITS> pCocofile((FITS*)0);
  std::unique_ptr<FITS> pLinefile((FITS*)0);

  try {
    pCocofile.reset(new FITS(coconame,CCfits::Read,hduName,false,hduKeys,primaryKeys,(int)1));
  } catch(...) {
    string msg = "Failed to read " + coconame;
    FunctionUtility::xsWrite(msg, 5);
    return(1);
  }
  try {
    pLinefile.reset(new FITS(linename,CCfits::Read,hduName,false,hduKeys,primaryKeys,(int)1));
  } catch(...) {
    string msg = "Failed to read " + linename;
    FunctionUtility::xsWrite(msg, 5);
    return(2);
  }


  // loop round requested temperature indices
  
  for (size_t iTemp=0; iTemp<TemperatureIndex.size(); iTemp++) {


    if ( TemperatureRecord[TemperatureIndex[iTemp]].ElementRecord.size() == 0 ) {

      int extNum = TemperatureIndex[iTemp] + 2;
      stringstream tst;
      tst << extNum;
      string sextNum = tst.str();

      try {
	pCocofile->read(extNum, true, hduKeys);
      } catch(...) {
	string msg = "Failed to read extension " + sextNum + " from " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(3);
      }
      try {
	pLinefile->read(extNum, true, hduKeys);
      } catch(...) {
	string msg = "Failed to read extension " + sextNum + " from " + linename;
	FunctionUtility::xsWrite(msg, 5);
	return(4);
      }

      ExtHDU& cocoExt = pCocofile->currentExtension();
      ExtHDU& lineExt = pLinefile->currentExtension();

      // set up temperature record and read temperature for this pair of
      // extensions

      ApedTemperatureRecord Trecord;

      Real Temp;
      try {
	cocoExt.readKey("TEMPERATURE", Temp);
      } catch(...) {
	string msg = "Failed to read TEMPERATURE keyword from extension " 
	  + sextNum + " from " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(5);
      }
      Trecord.Temperature = Temp / keVtoK;

      std::ostringstream oss;
      oss << "Reading data for temperature " << TemperatureIndex[iTemp] << " "
	  << Temperatures[TemperatureIndex[iTemp]] << " " << Trecord.Temperature
	  << "\n";
      FunctionUtility::xsWrite(oss.str(),25);

      // Read the line data into arrays.

      size_t NumberLines = (size_t)lineExt.rows();

      RealArray Lambda(NumberLines);
      RealArray LambdaError(NumberLines);
      RealArray Epsilon(NumberLines);
      RealArray EpsilonError(NumberLines);
      IntegerVector Element(NumberLines);
      IntegerVector Ion(NumberLines);
      IntegerVector ElementDriver(NumberLines);
      IntegerVector IonDriver(NumberLines);
      IntegerVector UpperLev(NumberLines);
      IntegerVector LowerLev(NumberLines);

      if ( NumberLines > 0 ) {

	try {
	  lineExt.column("Lambda").read(Lambda,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read Lambda column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("Lambda_Err").read(LambdaError,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read Lambda_Err column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("Epsilon").read(Epsilon,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read Epsilon column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("Epsilon_Err").read(EpsilonError,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read Epsilon_Err column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("Element").read(Element,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read Element column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("Ion").read(Ion,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read Ion column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("UpperLev").read(UpperLev,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read UpperLev column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("LowerLev").read(LowerLev,1,NumberLines);
	} catch(...) {
	  string msg = "Failed to read LowerLev column from extension " + sextNum + " of " + linename;
	  FunctionUtility::xsWrite(msg, 5);
	  return(6);
	}
	try {
	  lineExt.column("Elem_drv").read(ElementDriver,1,NumberLines);
	} catch(...) {
	  try {
	    lineExt.column("Element_drv").read(ElementDriver,1,NumberLines);
	  } catch(...) {
	    ElementDriver = Element;
	  }
	}
	try {
	  lineExt.column("Ion_drv").read(IonDriver,1,NumberLines);
	} catch(...) {
	  IonDriver = Ion;
	}

      }

      // Get the number of rows in the continuum extension

      size_t NumberRows = (size_t)cocoExt.rows();

      // First put the scalar columns into arrays

      IntegerVector Z(NumberRows);
      IntegerVector rmJ(NumberRows);
      IntegerVector Ncont(NumberRows);
      IntegerVector Npseudo(NumberRows);

      try {
	cocoExt.column("Z").read(Z, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read Z column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("rmJ").read(rmJ, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read rmJ column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("N_Cont").read(Ncont, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read N_cont column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("N_Pseudo").read(Npseudo, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read N_Pseudo column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }

      // and the vector columns into vectors of arrays

      vector<RealArray> E_Cont(NumberRows);
      vector<RealArray> Continuum(NumberRows);
      vector<RealArray> Cont_Err(NumberRows);
      vector<RealArray> E_Pseudo(NumberRows);
      vector<RealArray> Pseudo(NumberRows);
      vector<RealArray> Pseudo_Err(NumberRows);

      try {
	cocoExt.column("E_Cont").readArrays(E_Cont, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read E_Cont column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("Continuum").readArrays(Continuum, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read Continuum column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("Cont_Err").readArrays(Cont_Err, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read Cont_Err column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("E_Pseudo").readArrays(E_Pseudo, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read E_Pseudo column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("Pseudo").readArrays(Pseudo, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read Pseudo column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }
      try {
	cocoExt.column("Pseudo_Err").readArrays(Pseudo_Err, (long)1, (long)NumberRows);
      } catch(...) {
	string msg = "Failed to read Pseudo_Err column from extension " + sextNum + " of " + coconame;
	FunctionUtility::xsWrite(msg, 5);
	return(7);
      }

      // Set up a vector with unique elements included.

      vector<int> UniqueElements;
      for (size_t i=0; i<NumberLines; i++) {
	bool found(false);
	for (size_t j=0; j<UniqueElements.size(); j++) {
	  if ( UniqueElements[j] == Element[i] ) found = true;
	}
	if ( !found ) UniqueElements.push_back(Element[i]);
      }
      for (size_t i=0; i<NumberRows; i++) {
	bool found(false);
	for (size_t j=0; j<UniqueElements.size(); j++) {
	  if ( UniqueElements[j] == Z[i] ) found = true;
	}
	if ( !found ) UniqueElements.push_back(Z[i]);
      }

      oss.clear();
      oss << "Read " << NumberLines << " element lines and " << NumberRows << " "
	  << " continuum points for " << UniqueElements.size() << " unique elements.\n";
      FunctionUtility::xsWrite(oss.str(),25);

      // Now loop over the elemets

      for (size_t iElt=0; iElt<UniqueElements.size(); iElt++) {

	// Set up record for this element

	ApedElementRecord Erecord;
	Erecord.AtomicNumber = UniqueElements[iElt];

	// Set up vector of ions for this element. Need to run through both
	// the continuum and line data to generate a list of ions

	vector<int> Ions;
	for (size_t i=0; i<NumberRows; i++) {
	  if ( Z[i] == UniqueElements[iElt] ) {
	    Ions.push_back(rmJ[i]);
	  }
	}
	for (size_t i=0; i<NumberLines; i++) {
	  if ( Element[i] == UniqueElements[iElt] ) {
	    bool found = false;
	    for (size_t j=0; j<Ions.size(); j++) {
	      if ( Ion[i] == Ions[j] ) found = true;
	    }
	    if ( !found ) {
	      Ions.push_back(Ion[i]);
	    }
	  }
	}

	// Loop over ions for this element

	for (size_t iIon=0; iIon<Ions.size(); iIon++) {

	  // Set up record for this ion

	  ApedIonRecord IonRecord;
	  IonRecord.Ion = Ions[iIon];

	  // Initialize this ion in case we have not read in any data for it

	  IonRecord.ContinuumEnergy.resize(0);
	  IonRecord.ContinuumFlux.resize(0);
	  IonRecord.ContinuumFluxError.resize(0);
	  IonRecord.PseudoContinuumEnergy.resize(0);
	  IonRecord.PseudoContinuumFlux.resize(0);
	  IonRecord.PseudoContinuumFluxError.resize(0);

	  // Loop round the continuum and pseudo arrays loading the
	  // appropriate ones into IonRecord.

	  for (size_t i=0; i<NumberRows; i++) {
	    if ( Z[i] == UniqueElements[iElt] && rmJ[i] == Ions[iIon] ) {

	      size_t N = (size_t)Ncont[i];
	      IonRecord.ContinuumEnergy.resize(N);
	      const RealArray& E_Cont_i = E_Cont[i];
	      for (size_t j=0; j<N; ++j)
		IonRecord.ContinuumEnergy[j] = E_Cont_i[j];
	      IonRecord.ContinuumFlux.resize(N);
	      const RealArray& Continuum_i = Continuum[i];
	      for (size_t j=0; j<N; ++j)
		IonRecord.ContinuumFlux[j] = Continuum_i[j];
	      IonRecord.ContinuumFluxError.resize(N);
	      const RealArray& Cont_Err_i = Cont_Err[i];
	      for (size_t j=0; j<N; ++j)
		IonRecord.ContinuumFluxError[j] = Cont_Err_i[j];

	      N = (size_t)Npseudo[i];
	      IonRecord.PseudoContinuumEnergy.resize(N);
	      const RealArray& E_Pseudo_i = E_Pseudo[i];
	      for (size_t j=0; j<N; ++j)
		IonRecord.PseudoContinuumEnergy[j] = E_Pseudo_i[j];
	      IonRecord.PseudoContinuumFlux.resize(N);
	      const RealArray& Pseudo_i = Pseudo[i];
	      for (size_t j=0; j<N; ++j)
		IonRecord.PseudoContinuumFlux[j] = Pseudo_i[j];
	      IonRecord.PseudoContinuumFluxError.resize(N);
	      const RealArray& Pseudo_Err_i = Pseudo_Err[i];
	      for (size_t j=0; j<N; ++j)
		IonRecord.PseudoContinuumFluxError[j] = Pseudo_Err_i[j];
	      
	    }
	  }

	  // now find and load the emission lines for this element and ion
	  // first need to find the number of lines and size the arrays

	  size_t N = 0;
	  for (size_t i=0; i<NumberLines; i++) {
	    if ( Element[i] == Erecord.AtomicNumber && Ion[i] == IonRecord.Ion ) N++;
	  }
	  IonRecord.LineEnergy.resize(N);
	  IonRecord.LineEnergyError.resize(N);
	  IonRecord.LineEmissivity.resize(N);
	  IonRecord.LineEmissivityError.resize(N);
	  IonRecord.ElementDriver.resize(N);
	  IonRecord.IonDriver.resize(N);
	  IonRecord.UpperLevel.resize(N);
	  IonRecord.LowerLevel.resize(N);

	  // now loop through the lines

	  size_t iLine(0);
	  for (size_t i=0; i<NumberLines; i++) {
	    if ( Element[i] == Erecord.AtomicNumber && Ion[i] == IonRecord.Ion ) {
	      IonRecord.LineEnergy[iLine] = KEVTOA / Lambda[i];
	      IonRecord.LineEnergyError[iLine] = KEVTOA * LambdaError[i] / Lambda[i] / Lambda[i];
	      IonRecord.LineEmissivity[iLine] = Epsilon[i];
	      IonRecord.LineEmissivityError[iLine] = EpsilonError[i];
	      IonRecord.ElementDriver[iLine] = ElementDriver[i];
	      IonRecord.IonDriver[iLine] = IonDriver[i];
	      IonRecord.UpperLevel[iLine] = UpperLev[i];
	      IonRecord.LowerLevel[iLine] = LowerLev[i];
	      // test for funnies - energy negative or not increasing
	      if ( IonRecord.LineEnergy[iLine] < 0.0 ) {
		std::stringstream strmsg;
		strmsg << "WARNING: Wavelength negative for line " << iLine+1 << " of element " 
		       << Erecord.AtomicNumber << " and ion " << IonRecord.Ion << ".\n"
		       << "         Do not continue analysis without correcting HDU " 
		       << sextNum << " in " << linename << std::endl;
		FunctionUtility::xsWrite(strmsg.str(), 5);
		return(7);
	      }
	      if (iLine > 0 && 
		  IonRecord.LineEnergy[iLine] < IonRecord.LineEnergy[iLine-1]) {
		std::stringstream strmsg;
		strmsg << "WARNING: Wavelength for line " << iLine+1 << " (" 
		       << KEVTOA/IonRecord.LineEnergy[iLine] << ") is > than that for line "
		       << iLine << " (" << KEVTOA/IonRecord.LineEnergy[iLine-1] << ") of element " 
		       << Erecord.AtomicNumber << " and ion " << IonRecord.Ion << ".\n"
		       << "         Do not continue analysis without correcting row "
		       << i+1 << " of HDU " << sextNum << " in " << linename << std::endl;
		FunctionUtility::xsWrite(strmsg.str(), 5);
		return(7);
	      }
	      iLine++;
	    }
	  }


	  // load this ion record into the element record

	  Erecord.LoadIonRecord(IonRecord);

	}

	// load this element record into the temperature record

	Trecord.LoadElementRecord(Erecord);

      }

      // load this temperature record into the Aped object

      TemperatureRecord[TemperatureIndex[iTemp]] = Trecord;

    }

  }


  return(0);

}

// return the number of temperatures

int Aped::NumberTemperatures()
{
  return (int)Temperatures.size();
}

// return the number of elements

int Aped::NumberElements()
{
  return NelementsCoco[0];
}

// return the number of ions for the input element. Not all ions may be present
// in the data for all temperatures so need to loop through temperature records

int Aped::NumberIons(int Z)
{
  vector<int> Ions;
  for (size_t iT=0; iT<TemperatureRecord.size(); iT++) {
    ApedTemperatureRecord& thisT = TemperatureRecord[iT];
    for (size_t iElt=0; iElt<thisT.ElementRecord.size(); iElt++) {
      ApedElementRecord& thisElt = thisT.ElementRecord[iElt];
      if ( Z == thisElt.AtomicNumber ) {
	for (size_t iIon=0; iIon<thisElt.IonRecord.size(); iIon++) {
	  ApedIonRecord& thisIon = thisElt.IonRecord[iIon];
	  bool found = false;
	  for (size_t i=0; i<Ions.size(); i++) {
	    if ( Ions[i] == thisIon.Ion ) found = true;
	  }
	  if ( !found ) Ions.push_back(thisIon.Ion);
	}
      }
    }
  }

  return (int)Ions.size();
}

// return true if the data has been loaded for the requested temperature index

bool Aped::IsTemperatureLoaded(const int TemperatureIndex)
{
  bool loaded = false;
  if ( TemperatureRecord[TemperatureIndex].ElementRecord.size() > 0 ) loaded = true;
  return loaded;
}

// deep copy

Aped& Aped::operator=(const Aped& beta)
{
  TemperatureRecord.resize(beta.TemperatureRecord.size());
  for (size_t i=0; i<TemperatureRecord.size(); i++) {
    TemperatureRecord[i] = beta.TemperatureRecord[i];
  }
  Temperatures.resize(beta.Temperatures.size());
  Temperatures = beta.Temperatures;
  NelementsLine.resize(beta.NelementsLine.size());
  NelementsLine = beta.NelementsLine;
  NelementsCoco.resize(beta.NelementsCoco.size());
  NelementsCoco = beta.NelementsCoco;
  Nline.resize(beta.Nline.size());
  Nline = beta.Nline;
  Ncont.resize(beta.Ncont.size());
  Ncont = beta.Ncont;
  coconame = beta.coconame;
  linename = beta.linename;
  return *this;
}

// method to get spectrum for a single temperature and DEM. It
// assumes ionization equilibrium so the object must have been loaded from
// the ionization equilibrium input files

void Aped::SumEqSpectra(const RealArray& energyArray, 
			const IntegerVector& Zinput, const RealArray& abunZ,
			const Real Redshift, const Real& Tinput,
			const Real& Dem, RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  DemArr[0] = Dem;

  this->SumEqSpectra(energyArray, Zinput, abunZ, Redshift, TinputArr, DemArr,
		     fluxArray, fluxErrArray);

  return;
}

// method to sum spectra for the temperatures and associated DEMs. It
// assumes ionization equilibrium so the object must have been loaded from
// the ionization equilibrium input files

void Aped::SumEqSpectra(const RealArray& energyArray, 
			const IntegerVector& Zinput, const RealArray& abunZ,
			const Real Redshift, const RealArray& Tinput,
			const RealArray& Dem, RealArray& fluxArray, RealArray& fluxErrArray)
{
  // This routine loops over the temperatures requested, for each temperature
  // interpolating between the nearest temperatures tabulated. For the continuum 
  // and pseudo-continuum it also loops over elements, weighting by the 
  // abundances. The lines are added in with the line emissivity multiplied
  // by the abundance. 

  // Each member of the abundance array abunZ is the abundance for the
  // element whose atomic number is given by the matching member of the
  // Zinput array.

  // Calculate by calling SumSpectra with the IonFrac array set to
  // the Dem for each element.We are assuming that this Aped object is generated
  // for an equilibium ionization plasma so the ionization fractions are
  // implicitly included.

  vector<vector<RealArray> > IonFrac(Tinput.size());
  for (size_t it=0; it<Tinput.size(); it++) {
    IonFrac[it].resize(Zinput.size());
    for (size_t ielt=0; ielt<Zinput.size(); ielt++) {
      IonFrac[it][ielt].resize(Zinput[ielt]+2);
      for (size_t iion=0; iion<(size_t)Zinput[ielt]+2; iion++) {
        IonFrac[it][ielt][iion] = Dem[it];
      }
    }
  }

  RealArray Tbinput(Tinput.size());
  Tbinput = Tinput;

  this->SumSpectra(energyArray, Zinput, abunZ, Redshift, Tinput, Tbinput,
		   IonFrac, true, fluxArray, fluxErrArray);

  return;

}

// Cases for the temperature used for the thermal broadening differing from that
// for the ionization.

void Aped::SumEqSpectra(const RealArray& energyArray, 
			const IntegerVector& Zinput, const RealArray& abunZ,
			const Real Redshift, const Real& Tinput, const Real& Tbinput,
			const Real& Dem, RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray TbinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  TbinputArr[0] = Tbinput;
  DemArr[0] = Dem;

  this->SumEqSpectra(energyArray, Zinput, abunZ, Redshift, TinputArr, TbinputArr,
		     DemArr, fluxArray, fluxErrArray);

  return;
}


void Aped::SumEqSpectra(const RealArray& energyArray, 
			const IntegerVector& Zinput, const RealArray& abunZ,
			const Real Redshift, const RealArray& Tinput,
			const RealArray& Tbinput, const RealArray& Dem, 
			RealArray& fluxArray, RealArray& fluxErrArray)
{
  // Calculate by calling SumSpectra with the IonFrac array set to
  // the Dem for each element.We are assuming that this Aped object is generated
  // for an equilibrium ionization plasma so the ionization fractions are
  // implicitly included.

  vector<vector<RealArray> > IonFrac(Tinput.size());
  for (size_t it=0; it<Tinput.size(); it++) {
    IonFrac[it].resize(Zinput.size());
    for (size_t ielt=0; ielt<Zinput.size(); ielt++) {
      IonFrac[it][ielt].resize(Zinput[ielt]+2);
      for (size_t iion=0; iion<(size_t)Zinput[ielt]+2; iion++) {
        IonFrac[it][ielt][iion] = Dem[it];
      }
    }
  }

  this->SumSpectra(energyArray, Zinput, abunZ, Redshift, Tinput, Tbinput,
		   IonFrac, true, fluxArray, fluxErrArray);

  return;

}

// Now the NEI methods

void Aped::SumNeqSpectra(const RealArray& energyArray, 
			 const IntegerVector& Zinput, const RealArray& abunZ,
			 const Real Redshift, const Real& Tinput,
			 const vector<RealArray>& IonFrac,
			 RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  vector<vector<RealArray> > IonFracArr(1);
  TinputArr[0] = Tinput;
  IonFracArr[0].resize(IonFrac.size());
  for (size_t i=0; i<IonFrac.size(); i++) {
    IonFracArr[0][i].resize(IonFrac[i].size());
    for (size_t j=0; j<IonFrac[i].size(); j++) {
      IonFracArr[0][i][j] = IonFrac[i][j];
    }
  }

  this->SumNeqSpectra(energyArray, Zinput, abunZ, Redshift, TinputArr, IonFracArr,
		      fluxArray, fluxErrArray);

  return;
}

void Aped::SumNeqSpectra(const RealArray& energyArray, 
			 const IntegerVector& Zinput, const RealArray& abunZ,
			 const Real Redshift, const RealArray& Tinput,
			 const vector<vector<RealArray> >& IonFrac, 
			 RealArray& fluxArray, RealArray& fluxErrArray)
{
  // This routine loops over the temperatures requested, for each temperature
  // interpolating between the nearest temperatures tabulated. Within each
  // temperature it loops over elements and for each element over the ion stages
  // requested. Continua, pseudo-continua and lines are summed in weighted by
  // the abundance and ion fraction.

  // The abunZ input array is the abundance required for each element listed
  // in Zinput.

  // The IonFrac input variable is a vector of size Tinput.size() with
  // each element corresponding to a temperature in Tinput. Each element of
  // the vector is itself a vector of size Zinput.size() with each element
  // corresponding to an atomic number in Zinput. Each element of this inner
  // vector is itself a RealArray with each element being the ion fraction
  // and the index assumed to be the ion stage.

  RealArray Tbinput(Tinput.size());
  Tbinput = Tinput;

  this->SumSpectra(energyArray, Zinput, abunZ, Redshift, Tinput, Tbinput,
		   IonFrac, false, fluxArray, fluxErrArray);

  return;
}

// Cases for the temperature used for the thermal broadening differing from that
// for the ionization.


void Aped::SumNeqSpectra(const RealArray& energyArray, 
			 const IntegerVector& Zinput, const RealArray& abunZ,
			 const Real Redshift, const Real& Tinput, 
			 const Real& Tbinput,
			 const vector<RealArray>& IonFrac,
			 RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray TbinputArr(1);
  vector<vector<RealArray> > IonFracArr(1);
  TinputArr[0] = Tinput;
  TbinputArr[0] = Tbinput;
  IonFracArr[0].resize(IonFrac.size());
  for (size_t i=0; i<IonFrac.size(); i++) {
    IonFracArr[0][i].resize(IonFrac[i].size());
    for (size_t j=0; j<IonFrac[i].size(); j++) {
      IonFracArr[0][i][j] = IonFrac[i][j];
    }
  }

  this->SumNeqSpectra(energyArray, Zinput, abunZ, Redshift, TinputArr, TbinputArr,
		      IonFracArr, fluxArray, fluxErrArray);

  return;
}

void Aped::SumNeqSpectra(const RealArray& energyArray, 
			 const IntegerVector& Zinput, const RealArray& abunZ,
			 const Real Redshift, const RealArray& Tinput,
			 const RealArray& Tbinput,
			 const vector<vector<RealArray> >& IonFrac, 
			 RealArray& fluxArray, RealArray& fluxErrArray)
{
  // The abunZ input array is the abundance required for each element listed
  // in Zinput.

  // The IonFrac input variable is a vector of size Tinput.size() with
  // each element corresponding to a temperature in Tinput. Each element of
  // the vector is itself a vector of size Zinput.size() with each element
  // corresponding to an atomic number in Zinput. Each element of this inner
  // vector is itself a RealArray with each element being the ion fraction
  // and the index assumed to be the ion stage.

  this->SumSpectra(energyArray, Zinput, abunZ, Redshift, Tinput, Tbinput,
		   IonFrac, false, fluxArray, fluxErrArray);

  return;
}

void Aped::SumSpectra(const RealArray& energyArray, 
		      const IntegerVector& Zinput, const RealArray& abundance,
		      const Real Redshift, const RealArray& Tinput,
		      const RealArray& Tbinput, 
		      const vector<vector<RealArray> >& IonFrac, const bool isCEI,
		      RealArray& fluxArray, RealArray& fluxErrArray)
{
  // Routine to handle both CEI and NEI plasmas. In the CEI case the IonFrac
  // structure contains the DEM. For a CEI Aped object the continuum fluxes
  // for each element are contained in an object with Ion=0.

  // This routine loops over the temperatures requested, for each temperature
  // interpolating between the nearest temperatures tabulated. Within each
  // temperature it loops over elements and for each element over the ion stages
  // requested. Continua, pseudo-continua and lines are summed in weighted by
  // the abundance and ion fraction.

  // The abunZ input array is the abundance required for each element listed
  // in Zinput. This is relative to the abundances set by the abund command.
  // Inside SumSpectra these will be corrected to those relative to Anders &
  // Grevesse.

  // The IonFrac input variable is a vector of size Tinput.size() with
  // each element corresponding to a temperature in Tinput. Each element of
  // the vector is itself a vector of size Zinput.size() with each element
  // corresponding to an atomic number in Zinput. Each element of this inner
  // vector is itself a RealArray with each element being the ion fraction
  // and the index assumed to be the ion stage.

  // Note that IonFrac is assumed to be the same for Tinput and Tbinput. This is
  // fine for CEI but not for NEI. At the moment there are no NEI models which
  // have differing Tinput and Tbinput.

  // we have to convert from abundances relative to the currently defined
  // set to Anders & Grevesse, which are used internally.

  RealArray abunZ(Zinput.size());
  for (size_t i=0; i<Zinput.size(); i++) {
    abunZ[i] = abundance[i] * FunctionUtility::getAbundance(Zinput[i]) / 
      FunctionUtility::getAbundance("angr", Zinput[i]);
  }

  // Initialize output arrays

  fluxArray.resize(energyArray.size()-1);
  fluxErrArray.resize(0);
  fluxArray = 0.0;

  // Precalculate for speed
  Real z1 = 1. + Redshift;

  // Store the normalization factor to apply to the output flux
  Real normfactor = 1.0e14;
  if ( Redshift > -1.0 ) normfactor /= z1;

  // Divide the stored minimum line flux by the normalization factor
  // since we will compare to line fluxes before normalization
  Real minLineflux(minimumLinefluxForBroadening/normfactor);

  RealArray sourceFrameEnergy(energyArray.size());
  sourceFrameEnergy = energyArray * z1;

  // Calculate the interpolation of IonFrac onto the tabulated temperatures

  vector<vector<RealArray> > interpIonFrac;
  vector<vector<RealArray> > interpbIonFrac;
  vector<Real> coeffSum;
  vector<Real> coeffbSum;
  vector<int> interpTindex;
  vector<int> interpTbindex;
  RealArray tval = Temperatures;

  // Loop over input (electron) temperatures

  for (size_t it=0; it<Tinput.size(); it++) {

    int nmtval = tval.size();

    int itemp[2];

    // Find tabulated temperatures to be used in the interpolation and
    // calculate interpolation coefficients

    if ( Tinput[it] < tval[0] ) {
      std::ostringstream oss;
      oss << " Warning: temperature " << Tinput[it] << "< minimum in model file " << tval[0] << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
      itemp[0] = 0;
      itemp[1] = 0;
    } else if ( Tinput[it] > tval[nmtval-1] ) {
      std::ostringstream oss;
      oss << " Warning: temperature " << Tinput[it] << "> maximum in model file " << tval[nmtval-1] << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
      itemp[0] = nmtval-1;
      itemp[1] = nmtval-1;
    } else {
      itemp[1] = 0;
      while ( Tinput[it] >= tval[itemp[1]] ) itemp[1]++;
      itemp[0] = itemp[1] - 1;
      if ( Tinput[it] == tval[itemp[0]] ) itemp[1] = itemp[0];
    }

    Real tdelta;
    if ( logTempInterpolation ) {
      tdelta = log(tval[itemp[1]]) - log(tval[itemp[0]]);
    } else {
      tdelta = tval[itemp[1]] - tval[itemp[0]];
    }
    Real tcoeff[2];
    if ( tdelta >  0.0 ) {
      if ( logTempInterpolation ) {
	tcoeff[0] = (log(tval[itemp[1]]) - log(Tinput[it]))/tdelta;
	tcoeff[1] = (log(Tinput[it]) - log(tval[itemp[0]]))/tdelta;
      } else {
	tcoeff[0] = (tval[itemp[1]] - Tinput[it])/tdelta;
	tcoeff[1] = (Tinput[it] - tval[itemp[0]])/tdelta;
      }
    } else {
      tcoeff[0] = 1.0;
      tcoeff[1] = 0.0;
    }

    std::ostringstream oss;
    oss  << "For continuum interpolating between " << tval[itemp[0]]
	 << " and " << tval[itemp[1]] << " for " << Tinput[it] << "\n"
	 << "with coefficients " << tcoeff[0] << " and " << tcoeff[1] << "\n"; 
    FunctionUtility::xsWrite(oss.str(),15);

    // Loop over the two temperatures to use in interpolation

    for (size_t itabt=0; itabt<2; itabt++) {

      // Find whether this temperature index is already in interpTindex.

      int ifound = -1;
      for (size_t i=0; i<interpTindex.size(); i++) {
	if ( interpTindex[i] == itemp[itabt] ) ifound = i;
      }

      // if it is already there then add in the appropriate fraction of IonFrac
      // if not then push a new entry.

      if ( ifound != -1 ) {
	vector<RealArray>& rIonFrac = interpIonFrac[ifound];
	for (size_t i=0; i<rIonFrac.size(); i++) {
	  for (size_t j=0; j<rIonFrac[i].size(); j++) {
	    rIonFrac[i][j] += tcoeff[itabt] * IonFrac[it][i][j];
	  }
	}
	coeffSum[ifound] += tcoeff[itabt];
      } else {
	interpTindex.push_back(itemp[itabt]);
	vector<RealArray> tIonFrac(IonFrac[it].size());
	for (size_t i=0; i<tIonFrac.size(); i++) {
	  tIonFrac[i].resize(IonFrac[it][i].size());
	  tIonFrac[i] = tcoeff[itabt] * IonFrac[it][i];
	}
	interpIonFrac.push_back(tIonFrac);
	coeffSum.push_back(tcoeff[itabt]);
      }
    
    }

  }

  // Loop over input (ion) temperatures

  for (size_t it=0; it<Tbinput.size(); it++) {

    int nmtval = tval.size();

    int itemp[2];

    // Find tabulated temperatures to be used in the interpolation and
    // calculate interpolation coefficients

    if ( Tbinput[it] < tval[0] ) {
      std::ostringstream oss;
      oss << " Warning: temperature " << Tbinput[it] << "< minimum in model file " << tval[0] << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
      itemp[0] = 0;
      itemp[1] = 0;
    } else if ( Tbinput[it] > tval[nmtval-1] ) {
      std::ostringstream oss;
      oss << " Warning: temperature " << Tbinput[it] << "> maximum in model file " << tval[nmtval-1] << "\n";
      FunctionUtility::xsWrite(oss.str(),10);
      itemp[0] = nmtval-1;
      itemp[1] = nmtval-1;
    } else {
      itemp[1] = 0;
      while ( Tbinput[it] >= tval[itemp[1]] ) itemp[1]++;
      itemp[0] = itemp[1] - 1;
      if ( Tbinput[it] == tval[itemp[0]] ) itemp[1] = itemp[0];
    }

    Real tdelta;
    if ( logTempInterpolation ) {
      tdelta = log(tval[itemp[1]]) - log(tval[itemp[0]]);
    } else {
      tdelta = tval[itemp[1]] - tval[itemp[0]];
    }
    Real tcoeff[2];
    if ( tdelta >  0.0 ) {
      if ( logTempInterpolation ) {
	tcoeff[0] = (log(tval[itemp[1]]) - log(Tbinput[it]))/tdelta;
	tcoeff[1] = (log(Tbinput[it]) - log(tval[itemp[0]]))/tdelta;
      } else {
	tcoeff[0] = (tval[itemp[1]] - Tbinput[it])/tdelta;
	tcoeff[1] = (Tbinput[it] - tval[itemp[0]])/tdelta;
      }
    } else {
      tcoeff[0] = 1.0;
      tcoeff[1] = 0.0;
    }

    std::ostringstream oss;
    oss  << "For lines interpolating between " << tval[itemp[0]]
	 << " and " << tval[itemp[1]] << " for " << Tbinput[it] << "\n"
	 << "with coefficients " << tcoeff[0] << " and " << tcoeff[1] << "\n"; 
    FunctionUtility::xsWrite(oss.str(),15);

    // Loop over the two temperatures to use in interpolation

    for (size_t itabt=0; itabt<2; itabt++) {

      // Find whether this temperature index is already in interpTbindex.

      int ifound = -1;
      for (size_t i=0; i<interpTbindex.size(); i++) {
	if ( interpTbindex[i] == itemp[itabt] ) ifound = i;
      }

      // if it is already there then add in the appropriate fraction of IonFrac
      // if not then push a new entry.

      if ( ifound != -1 ) {
	vector<RealArray>& rIonFrac = interpbIonFrac[ifound];
	for (size_t i=0; i<rIonFrac.size(); i++) {
	  for (size_t j=0; j<rIonFrac[i].size(); j++) {
	    rIonFrac[i][j] += tcoeff[itabt] * IonFrac[it][i][j];
	  }
	}
	coeffbSum[ifound] += tcoeff[itabt];
      } else {
	interpTbindex.push_back(itemp[itabt]);
	vector<RealArray> tIonFrac(IonFrac[it].size());
	for (size_t i=0; i<tIonFrac.size(); i++) {
	  tIonFrac[i].resize(IonFrac[it][i].size());
	  tIonFrac[i] = tcoeff[itabt] * IonFrac[it][i];
	}
	interpbIonFrac.push_back(tIonFrac);
	coeffbSum.push_back(tcoeff[itabt]);
      }
    
    }

  }

  // interpTindex and interpTbindex now contain the list of temperature indices 
  // we will need. First find out if we need to read the data for any of these 
  // from the files

  vector<int> needTindex;
  for (size_t it=0; it<interpTindex.size(); it++) {
    if (!IsTemperatureLoaded(interpTindex[it])) needTindex.push_back(interpTindex[it]);
  }
  for (size_t it=0; it<interpTbindex.size(); it++) {
    if (!IsTemperatureLoaded(interpTbindex[it])) needTindex.push_back(interpTbindex[it]);
  }
  if (!needTindex.empty()) ReadTemperature(needTindex);

  // First do the continuum by looping over the tabulated temperatures required 
  // using the indices stored in interpTindex and the ion fractions in interpIonFrac.

  for (size_t it=0; it<interpTindex.size(); it++) {

    ApedTemperatureRecord& thisT = TemperatureRecord[interpTindex[it]];

    // Loop over elements listed in Zinput.

    for (size_t ielt=0; ielt<Zinput.size(); ielt++) {

      // temporary arrays to accumulate flux for continuum and pseudo-continuum
      RealArray contFluxArray(0.0,energyArray.size()-1);
      RealArray pseudoFluxArray(0.0,energyArray.size()-1);

      // Find the matching element in thisT.

      int iindex = -1;
      for (size_t i=0; i<thisT.ElementRecord.size(); i++) {
	if ( thisT.ElementRecord[i].AtomicNumber == Zinput[ielt] ) {
	  iindex = i;
	}
      }

      // If there is no tabulated information for this element at this
      // temperature then move on to next element

      if ( iindex == -1 ) continue;

      ApedElementRecord& thisElt = thisT.ElementRecord[iindex];
      size_t eltZ = thisElt.AtomicNumber;

      // Set the array over ions for this temperature and element.
      
      const RealArray& IonFracElt = interpIonFrac[it][ielt];

      // Loop over ions for this temperature and element in IonFrac
      // setting up the IonIndex list

      size_t nions = IonFracElt.size();

      vector<int> IonIndex;

      for (size_t iion=0; iion<nions+1; iion++) {

	// Find the matching ion in thisElt. Note the special case
	// of Ion=0 which is for the CEI continuum which is labelled as Ion=0.

	for (size_t i=0; i<thisElt.IonRecord.size(); i++) {
	  if ( thisElt.IonRecord[i].Ion == (int)iion ) {
	    IonIndex.push_back(i);
	  }
	}

      }

      // Now loop over the ions required calculating the continuum. If this
      // is CEI then just use the first element of IonIndex

      Real acoeff = abunZ[ielt];
      size_t nindex = IonIndex.size();
      if ( isCEI ) nindex = 1;

      for (size_t jindex=0; jindex<nindex; jindex++) {

	ApedIonRecord& thisIon = thisElt.IonRecord[IonIndex[jindex]];
	size_t iion = thisIon.Ion;

	Real coeff = acoeff;
	if ( iion == 0 ) {
	  if ( isCEI ) {
	    coeff *= IonFracElt[0];
	  } else {
	    coeff = 0.0;
	  }
	} else {
	  coeff *= IonFracElt[iion-1];
	}

	if ( coeff > 0.0 ) {

	  // Interpolate the continuum onto the output energy array. Note that 
	  // the input is flux/keV values at the energies tabulated so need to 
	  // multiply by the bin width to convert to the units that XSPEC 
	  // requires. Also note that this algorithm assumes that the response 
	  // energies are monotonically increasing.

	  apedInterpFlux(thisIon.ContinuumEnergy, thisIon.ContinuumFlux, 
			 coeff, sourceFrameEnergy, contFluxArray);

	  // Repeat for the pseudo-continuum.

	  if ( !noLines ) {
	    apedInterpFlux(thisIon.PseudoContinuumEnergy, 
			   thisIon.PseudoContinuumFlux, coeff,
			   sourceFrameEnergy, pseudoFluxArray);
	  }

	}

	// End loop over ions

      }

      // sum continuum flux for this temperature and element into total flux array

      fluxArray += contFluxArray;

      // Repeat for the pseudo-continuum. If the pseudo-continuum is to be
      // broadened then we do that by treating each energy bin as its own line.

      if ( broadenPseudoContinuum && (thermalBroadening || velocityBroadening>0.0) ) {

	// Set line widths and related quantities for any lines from this element

	Real widthCoeff;
	if ( thermalBroadening ) {
	  widthCoeff = sqrt( thisT.Temperature * KEVTOERG /
			     getAtomicMass(eltZ) / AMU + 
			     velocityBroadening*velocityBroadening*km_cm*km_cm )
	    / (LIGHTSPEED*km_cm);
	} else {
	  widthCoeff = velocityBroadening / LIGHTSPEED;
	}

	RealArray midEnergy(sourceFrameEnergy.size()-1);
	for (size_t i=0; i<midEnergy.size(); i++) {
	  midEnergy[i] = 0.5*(sourceFrameEnergy[i]+sourceFrameEnergy[i+1]);
	}

	std::vector<RealArray> width(midEnergy.size());
	for (size_t i=0; i<midEnergy.size(); i++) width[i].resize(1);
	for (size_t i=0; i<midEnergy.size(); i++) width[i][0] = midEnergy[i]*widthCoeff;

	RealArray pseudoLineNorms(pseudoFluxArray);
	pseudoFluxArray = 0.0;
	calcManyLines(sourceFrameEnergy, midEnergy, width, pseudoLineNorms,
		      (Real)1.0e-6, 0, true, pseudoFluxArray);

      }

      fluxArray += pseudoFluxArray;

      // End loop over elements

    }

    // End loop over temperatures.

  }

  // Now do the lines by looping over the tabulated temperatures required using
  // the indices stored in interpTbindex and the ion fractions in interpbIonFrac.

  Real maxLineflux = 0.0;
  
  for (size_t it=0; it<interpTbindex.size(); it++) {

    ApedTemperatureRecord& thisT = TemperatureRecord[interpTbindex[it]];

    // Loop over elements listed in Zinput.

    for (size_t ielt=0; ielt<Zinput.size(); ielt++) {

      // Find the matching element in thisT.

      int iindex = -1;
      for (size_t i=0; i<thisT.ElementRecord.size(); i++) {
	if ( thisT.ElementRecord[i].AtomicNumber == Zinput[ielt] ) {
	  iindex = i;
	}
      }

      // If there is no tabulated information for this element at this
      // temperature then move on to next element

      if ( iindex == -1 ) continue;

      ApedElementRecord& thisElt = thisT.ElementRecord[iindex];
      size_t eltZ = thisElt.AtomicNumber;

      // Set line widths and related quantities for any lines from this
      // element

      Real crtsig = 6.0;
      Real widthCoeff;
      if ( thermalBroadening ) {
	widthCoeff = sqrt( thisT.Temperature * KEVTOERG /
			   getAtomicMass(eltZ) / AMU + 
			   velocityBroadening*velocityBroadening*km_cm*km_cm ) / (LIGHTSPEED*km_cm);
      } else {
	widthCoeff = velocityBroadening / LIGHTSPEED;
      }

      Real crit = widthCoeff * crtsig;
      Real Emin = sourceFrameEnergy[0] - crit;
      Real Emax = sourceFrameEnergy[sourceFrameEnergy.size()-1] + crit;

      // Set the array over ions for this temperature and element.
      
      const RealArray& IonFracElt = interpbIonFrac[it][ielt];

      // Loop over ions for this temperature and element in IonFrac
      // setting up the IonIndex list

      size_t nions = IonFracElt.size();

      vector<int> IonIndex;

      for (size_t iion=0; iion<nions+1; iion++) {

	// Find the matching ion in thisElt. Note the special case
	// of Ion=0 which is for the CEI continuum which is labelled as Ion=0.

	for (size_t i=0; i<thisElt.IonRecord.size(); i++) {
	  if ( thisElt.IonRecord[i].Ion == (int)iion ) {
	    IonIndex.push_back(i);
	  }
	}

      }

      Real acoeff = abunZ[ielt];

      // Now do the lines. We do not include IonIndex[0] because that is
      // the special case for CEI continuum.

      for (size_t jindex=1; jindex<IonIndex.size(); jindex++) {

	ApedIonRecord& thisIon = thisElt.IonRecord[IonIndex[jindex]];

	// Now add in the lines if there are any

	bool haveLines = false;
	size_t nLines = thisIon.LineEnergy.size();
	if ( nLines  > 0 ) {
	  if ( thisIon.LineEnergy[0] <= Emax && 
	       thisIon.LineEnergy[nLines-1] >= Emin ) {
	    haveLines = true;
	  }
	}

	if ( haveLines && !noLines && acoeff > 0.0 ) {

	  // Find the first and last lines for this ion within the
	  // required energy range

	  int firstLine = Numerics::BinarySearch(thisIon.LineEnergy, Emin) + 1;
	  int lastLine = Numerics::BinarySearch(thisIon.LineEnergy, Emax);

	  if ( lastLine < 0 ) lastLine = nLines-1;

	  // note that that the width array is set up as the vector<RealArray>
	  // input required for calcManyLines even though the RealArray has
	  // size 1. this is to allow generalization to line shapes with multiple
	  // shape parameters such as Voigt.

	  size_t nLinesNeeded = lastLine - firstLine + 1;
	  RealArray energy(nLinesNeeded);
	  std::vector<RealArray> width(nLinesNeeded);
	  for (size_t i=0; i<nLinesNeeded; i++) width[i].resize(1);
	  RealArray lineflux(nLinesNeeded);

	  for (size_t i=0; i<nLinesNeeded; i++) {

	    size_t il = i + firstLine;
	    energy[i] = thisIon.LineEnergy[il];

	    // Multiply by the ion fraction for the driver ion
	    lineflux[i] = acoeff * IonFracElt[thisIon.IonDriver[il]-1] * 
	      thisIon.LineEmissivity[il];

	    // If the no resonance line flag is set then zero out the line
	    // flux if appropriate - this doesn't do anything at the moment
	    // since we have no way of identifying resonance lines

	    if ( noResonanceLines ) {
	    }

	    if ( lineflux[i] > maxLineflux ) maxLineflux = lineflux[i];

	    // check whether the line flux is above the minimum specified by
	    // the user for broadening

	    if ( lineflux[i] >= minLineflux ) {

	      // Find the width of the line using the equation
	      //  W = E/c  sqrt [ kT / (A mp ) + velocity^2 ] if thermal broadening
	      // is on else
	      //  W = (E/c) velocity

	      width[i][0] = energy[i] * widthCoeff;

	    } else {

	      width[i][0] = 0.0;

	    }

	    // End set up loop over lines.

	  }

	  // calculate all the lines for this element and ion in one go for efficiency
	  // The calcManyLines routine can be used for different line shapes. 0 
	  // indicates Gaussian and 1 Lorentzian.

	  calcManyLines(sourceFrameEnergy, energy, width, lineflux, 
			(Real)1.0e-6, 0, true, fluxArray);

	  // End if statement on #lines > 0

	}

	//End loop over ions

      }

      // End loop over elements.

    }

    // End loop over temperatures.

  }

  // Fix the normalization by multiplying by 1e14 and correcting for
  // time dilation and stored the max line flux

  FunctionUtility::loadDbValue("ApecMaxLineflux", maxLineflux*normfactor);
  
  fluxArray *= normfactor;

  return;
}

// methods for ApedTemperatureRecord

ApedTemperatureRecord::ApedTemperatureRecord()
{
}

ApedTemperatureRecord::~ApedTemperatureRecord()
{
}

void ApedTemperatureRecord::LoadElementRecord(ApedElementRecord input)
{
  ElementRecord.push_back(input);
  return;
}

// deep copy

ApedTemperatureRecord& ApedTemperatureRecord::operator=(const ApedTemperatureRecord& beta)
{
  Temperature = beta.Temperature;
  ElementRecord.resize(beta.ElementRecord.size());
  for (size_t i=0; i<ElementRecord.size(); i++) {
    ElementRecord[i] = beta.ElementRecord[i];
  }
  return *this;
}

// methods for ApedElementRecord

ApedElementRecord::ApedElementRecord()
{
}

ApedElementRecord::~ApedElementRecord()
{
}

void ApedElementRecord::LoadIonRecord(ApedIonRecord input)
{
  IonRecord.push_back(input);
  return;
}

// deep copy

ApedElementRecord& ApedElementRecord::operator=(const ApedElementRecord& beta)
{
  AtomicNumber = beta.AtomicNumber;
  IonRecord.resize(beta.IonRecord.size());
  for (size_t i=0; i<IonRecord.size(); i++) {
    IonRecord[i] = beta.IonRecord[i];
  }
  return *this;
}

// methods for ApedIonRecord

ApedIonRecord::ApedIonRecord()
{
}

ApedIonRecord::~ApedIonRecord()
{
}

// deep copy

ApedIonRecord& ApedIonRecord::operator=(const ApedIonRecord& beta)
{
  Ion = beta.Ion;
  ContinuumEnergy.resize(beta.ContinuumEnergy.size());
  ContinuumEnergy = beta.ContinuumEnergy;
  ContinuumFlux.resize(beta.ContinuumFlux.size());
  ContinuumFlux = beta.ContinuumFlux;
  ContinuumFluxError.resize(beta.ContinuumFluxError.size());
  ContinuumFluxError = beta.ContinuumFluxError;
  PseudoContinuumEnergy.resize(beta.PseudoContinuumEnergy.size());
  PseudoContinuumEnergy = beta.PseudoContinuumEnergy;
  PseudoContinuumFlux.resize(beta.PseudoContinuumFlux.size());
  PseudoContinuumFlux = beta.PseudoContinuumFlux;
  PseudoContinuumFluxError.resize(beta.PseudoContinuumFluxError.size());
  PseudoContinuumFluxError = beta.PseudoContinuumFluxError;
  LineEnergy.resize(beta.LineEnergy.size());
  LineEnergy = beta.LineEnergy;
  LineEnergyError.resize(beta.LineEnergyError.size());
  LineEnergyError = beta.LineEnergyError;
  LineEmissivity.resize(beta.LineEmissivity.size());
  LineEmissivity = beta.LineEmissivity;
  LineEmissivityError.resize(beta.LineEmissivityError.size());
  LineEmissivityError = beta.LineEmissivityError;
  ElementDriver.resize(beta.ElementDriver.size());
  ElementDriver = beta.ElementDriver;
  IonDriver.resize(beta.IonDriver.size());
  IonDriver = beta.IonDriver;
  UpperLevel.resize(beta.UpperLevel.size());
  UpperLevel = beta.UpperLevel;
  LowerLevel.resize(beta.LowerLevel.size());
  LowerLevel = beta.LowerLevel;

  return *this;
}

// Handy routines

Real getAtomicMass(const int& AtomicNumber)
{
  const int NUMZ = 30;

  const Real atomM[NUMZ] = {1.00794, 4.002602, 6.941, 9.012182, 10.811,
			    12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
			    22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
			    32.065, 35.4527, 39.948, 39.0983, 40.078,
			    44.955910, 47.867, 50.9415, 51.9961, 54.938049,
			    55.845, 58.933200, 58.6934, 63.456, 65.39};

  return atomM[AtomicNumber-1];

}

// return the filenames if the version set by APECROOT is the same as that
// entered return true

bool getApedFileNames(string& version, string& continuumFile, string& lineFile)
{
  string pname = "APECROOT";
  string newVersion = FunctionUtility::getModelString(pname);
  if ( !newVersion.length() || newVersion == FunctionUtility::NOT_A_KEY() )
    newVersion = version;

  bool sameVersion = (newVersion == version);
  version = newVersion;

  // set the names of the input files

  // first test whether only a version number has been set. Then assume 
  // that only a root filename has been given. If neither these work then 
  // assume that an entire directory path was included.

  const string& datadir = FunctionUtility::modelDataPath();
  string rootfil = datadir + "apec_v" + version;
  continuumFile = rootfil + "_coco.fits";
  lineFile = rootfil + "_line.fits";
  std::ifstream test1(continuumFile.c_str());
  if (!test1) {
    rootfil = datadir + version;
    continuumFile = rootfil + "_coco.fits";
    lineFile = rootfil + "_line.fits";
    std::ifstream test2(continuumFile.c_str());
    if (!test2) {
      rootfil = version;
      continuumFile = rootfil + "_coco.fits";
      lineFile = rootfil + "_line.fits";            
    }
  }
  return sameVersion;
}


bool getNEIApedFileNames(string& version, string& eigenVersion, string& continuumFile, string& lineFile)
{

  // check for NEIVERS xset
  string pname = "NEIVERS";
  string newEigenVersion = FunctionUtility::getModelString(pname);
  if ( !newEigenVersion.length() || newEigenVersion == FunctionUtility::NOT_A_KEY() )
    newEigenVersion = eigenVersion;

  // Check for a valid eigenVersion string
  if ( newEigenVersion != "1.0" && newEigenVersion != "1.1" && 
       newEigenVersion != "2.0" && newEigenVersion.substr(0,2) != "3." ) {
    std::ostringstream oss;
    oss << newEigenVersion << " is not a valid NEIVERS, ignoring." << "\n";
    FunctionUtility::xsWrite(oss.str(),10);
    newEigenVersion = eigenVersion;
  }

  bool sameEigenVersion = (newEigenVersion == eigenVersion);
  eigenVersion = newEigenVersion;

  // If eigenVersion is 1.* then this is the old code which doesn't need filenames
  // so return immediately
  if ( eigenVersion.substr(0,1) == "1" ) return sameEigenVersion;

  // Check for NEIAPECROOT xset
  pname = "NEIAPECROOT";
  string newVersion = FunctionUtility::getModelString(pname);
  if ( !newVersion.length() || newVersion == FunctionUtility::NOT_A_KEY() )
    newVersion = version;

  bool sameVersion = (newVersion == version);
  string oldVersion(version);
  version = newVersion;

  // set the names of the input files. if eigenVersion is 2.0 then special case
  // otherwise use the version number to construct the filenames

  const string& datadir = FunctionUtility::modelDataPath();

  if ( eigenVersion == "2.0" ) {

    continuumFile = datadir + "APEC_nei_v11_comp.fits";
    lineFile = datadir + "APEC_nei_v11_line.fits";

  } else {

    // first test whether only a version number has been set. Then assume 
    // that only a root filename has been given. If neither these work then 
    // assume that an entire directory path was included. If that doesn't work
    // try assuming the user included the _nei in the name. If even that doesn't
    // work then just stay with the old version.

    string rootfil = datadir + "apec_v" + version;
    continuumFile = rootfil + "_nei_comp.fits";
    lineFile = rootfil + "_nei_line.fits";
    std::ifstream test1(continuumFile.c_str());
    if (!test1) {
      rootfil = datadir + version;
      continuumFile = rootfil + "_nei_comp.fits";
      lineFile = rootfil + "_nei_line.fits";
      std::ifstream test2(continuumFile.c_str());
      if (!test2) {
	rootfil = version;
	continuumFile = rootfil + "_nei_comp.fits";
	lineFile = rootfil + "_nei_line.fits";
	std::ifstream test3(continuumFile.c_str());
	if (!test3) {
	  continuumFile = rootfil + "_comp.fits";
	  lineFile = rootfil + "_line.fits";
	  std::ifstream test4(continuumFile.c_str());
	  if (!test4) {
	    std::ostringstream oss;
	    oss << "Cannot find any NEI APEC files corresponding to the root\n"
		<< rootfil << "\n"
		<< "continuing to use those from " << oldVersion << "\n";
	    FunctionUtility::xsWrite(oss.str(),10);
	    sameVersion = true;
	  }
	}
      }
    }
  }
  return (sameVersion && sameEigenVersion);
}

// Interpolate from inputEnergy, inputFlux onto energyArray, fluxArray.
// Used to interpolate both continuum and pseudo-continuum.

void apedInterpFlux(const RealArray& inputEnergy, const RealArray& inputFlux,
		    const Real& coeff, const RealArray& energyArray, 
		    RealArray& fluxArray)
{

  if ( inputEnergy.size() <= 0 ) return;
  if ( energyArray.size() <= 0 ) return;

  size_t inputLastE = inputEnergy.size() - 1;
  size_t lastE = energyArray.size() - 1;

  if ( energyArray[0] > inputEnergy[inputLastE] ) return;
  if ( energyArray[lastE] < inputEnergy[0] ) return;

  size_t ilow = 0;
  for (size_t ien=0; ien<energyArray.size()-1; ien++) {

    Real emin = energyArray[ien];
    Real emax = energyArray[ien+1];

    // if there are no input points in this bin then go onto the next bin

    if ( emin >= inputEnergy[inputLastE] || emax <= inputEnergy[0] ) continue;

    // find last input point below the start of the bin

    while ( inputEnergy[ilow] < emin ) ilow++;
    if ( ilow > 0 ) ilow--;
    
    Real eInLow = inputEnergy[ilow];
    Real eInHigh = inputEnergy[ilow+1];
    Real eInRange = eInHigh - eInLow;

    // Case of first input point being in this bin - assume flux density
    // is equal to that at first point below first point
    if ( eInLow > emin ) {
      fluxArray[ien] += coeff * inputFlux[ilow] * ( eInLow - emin );
    } else {
      // ilow is below start of energy bin so now do case when ilow+1 is above
      // the end of the energy bin
      if ( eInHigh > emax ) {
	Real emiddle = (emin+emax)*0.5;
	fluxArray[ien] += coeff * (emax-emin) *
	  ( inputFlux[ilow] * (eInHigh-emiddle) + inputFlux[ilow+1] * (emiddle-eInLow) )
	  / eInRange;
      } else {
	// case of ilow+1 within this bin
	Real emiddle = (emin+eInHigh)*0.5;
	fluxArray[ien] += coeff * (eInHigh-emin) *
	  ( inputFlux[ilow] * (eInHigh-emiddle) + inputFlux[ilow+1] * (emiddle-eInLow) )
	  / eInRange;
	// now run through input points until ilow is above the end of
	// energy bin
	ilow++;
        if (ilow == inputLastE)
        {
	   fluxArray[ien] += coeff * inputFlux[ilow] * ( emax - eInHigh);
           continue;
        }
	eInLow = eInHigh;
	eInHigh = inputEnergy[ilow+1];
	eInRange = eInHigh - eInLow;
	while ( eInLow <= emax && ilow < inputLastE ) {
	  if ( eInHigh <= emax ) {
	    emiddle = (eInLow+eInHigh)*0.5;
	    fluxArray[ien] += coeff * (eInHigh-eInLow) *
	      ( inputFlux[ilow] * (eInHigh-emiddle) + inputFlux[ilow+1] * (emiddle-eInLow) )
	      / eInRange;
	  } else {
	    emiddle = (eInLow+emax)*0.5;
	    fluxArray[ien] += coeff * (emax-eInLow) *
	      ( inputFlux[ilow] * (eInHigh-emiddle) + inputFlux[ilow+1] * (emiddle-eInLow) )
	      / eInRange;
	  }	    
	  ilow++;
          if (ilow == inputLastE)
          {
	     fluxArray[ien] += coeff * inputFlux[ilow] * ( emax - eInHigh);
          }
          else
          {
	     eInLow = eInHigh;
	     eInHigh = inputEnergy[ilow+1];
	     eInRange = eInHigh - eInLow;
          }
	}
	if (ilow < inputLastE) {
	  // ilow is now above the end of the current bin so
	  // decrease it by 1 so that it will be below the start of the next
	  // bin
	  ilow--;
	  eInLow = inputEnergy[ilow];
	  eInHigh = inputEnergy[ilow+1];
	  eInRange = eInHigh - eInLow;
	}
      }
    }

    // end loop over output energy bins
  }

  return;
}

// Wrap-ups to read the file and calculate the spectrum

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  DemArr[0] = Dem;

  return calcCEISpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
  			 DemArr, qtherm, velocity, fluxArray, fluxErrArray);
}

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray Tbinput(Tinput.size());
  Tbinput = Tinput;

  return calcCEISpectrum(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput,
  			 Dem, qtherm, velocity, fluxArray, fluxErrArray);
}

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput, 
		    const Real& Tbinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray TbinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  TbinputArr[0] = Tbinput;
  DemArr[0] = Dem;

  return calcCEISpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
  			 TbinputArr, DemArr, qtherm, velocity, fluxArray, fluxErrArray);
}


int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput, 
		    const RealArray& Tbinput,
		    const RealArray& Dem, const bool qtherm, 
		    const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  static bool noLines(false);
  return calcCEISpectrum(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput,
			 Dem, qtherm, velocity, noLines, fluxArray, fluxErrArray);
}

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    const bool noLines, RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  DemArr[0] = Dem;

  return calcCEISpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
  			 DemArr, qtherm, velocity, noLines, fluxArray, fluxErrArray);
}

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Dem, const bool qtherm, const Real velocity,
		    const bool noLines, RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray Tbinput(Tinput.size());
  Tbinput = Tinput;

  return calcCEISpectrum(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput,
  			 Dem, qtherm, velocity, noLines, fluxArray, fluxErrArray);
}

int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput, 
		    const Real& Tbinput,
		    const Real& Dem, const bool qtherm, const Real velocity,
		    const bool noLines, RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray TbinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  TbinputArr[0] = Tbinput;
  DemArr[0] = Dem;

  return calcCEISpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
  			 TbinputArr, DemArr, qtherm, velocity, noLines, fluxArray, fluxErrArray);
}


int calcCEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput, 
		    const RealArray& Tbinput,
		    const RealArray& Dem, const bool qtherm, 
		    const Real velocity, const bool noLines,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  static bool isFirst(true);
  static Aped AtomicData;
  static string version = FunctionUtility::atomdbVersion();

  // check the setting of APECUSENEI to see whether or not we should route this
  // calculation through the NEI code to avoid potential problems interpolating on
  // the CEI tables

  string pname = "APECUSENEI";
  string pvalue = FunctionUtility::getModelString(pname);
  bool useNEI = false;
  if (pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY()) {
    pvalue = XSutility::lowerCase(pvalue);
    if (pvalue.substr(0,1) == "t" || pvalue.substr(0,1) == "y" || 
	pvalue.substr(0,2) == "on") {
      useNEI = true;
    } else {
      useNEI = false;
    }
  }

  if ( useNEI ) {

    // routing this calculation through calcNEISpectrum

    vector<vector<RealArray> > IonFrac(Tinput.size());
    for (size_t iT=0; iT<Tinput.size(); iT++) {
      calcCEIfractions(Tinput[iT], Zinput, IonFrac[iT]);
      for (size_t iZ=0; iZ<IonFrac[iT].size(); iZ++) IonFrac[iT][iZ] *= Dem[iT];
    }

    int status = calcNEISpectrum(energyArray, Zinput, abundance, Redshift, Tinput, 
				 Tbinput, IonFrac, qtherm, velocity, fluxArray, 
				 fluxErrArray);

    return status;

  } else {

    string contfil, linefil;
    bool noChange = getApedFileNames(version, contfil, linefil);

    if ( isFirst || !noChange ) {

      std::ostringstream oss;
      oss << "Reading APEC data from " << version << "\n";
      FunctionUtility::xsWrite(oss.str(),10);

      std::ostringstream oss2;
      oss2 << "Reading continuum data from " << contfil << "\n"
	   << "Reading line      data from " << linefil << "\n";
      FunctionUtility::xsWrite(oss2.str(),25);

      int status = AtomicData.Read(contfil, linefil);
      if ( status != 0 ) {
	if ( status == 1 ) {
	  string msg = "Failed to read " + contfil;
	  FunctionUtility::xsWrite(msg, 5);
	} else if ( status == 2 ) {
	  string msg = "Failed to read " + linefil;
	  FunctionUtility::xsWrite(msg, 5);
	}
	return status;
      }
      isFirst = false;
    }

    // Check for values of APECTHERMAL, APECVELOCITY, APECMINFLUX,
    // APECNOLINES and APECNORES, APECLOGTINTERP xset parameters

    AtomicData.SetLineSpecs(noLines, false, qtherm, velocity, 0.0, false, false);

    // calculate the spectrum

    AtomicData.SumEqSpectra(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput, Dem,
			    fluxArray, fluxErrArray);
  }

  return 0;
}

int calcRSSpectrum(const RealArray& energyArray, 
		   const IntegerVector& Zinput, const RealArray& abundance,
		   const Real Redshift, const Real& Tinput,
		   const Real& Dem, const bool qtherm, const Real velocity,
		   RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  DemArr[0] = Dem;

  return calcRSSpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
			DemArr, qtherm, velocity, fluxArray, fluxErrArray);
}

int calcRSSpectrum(const RealArray& energyArray, 
		   const IntegerVector& Zinput, const RealArray& abundance,
		   const Real Redshift, const RealArray& Tinput,
		   const RealArray& Dem, const bool qtherm, const Real velocity,
		   RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray Tbinput(Tinput.size());
  Tbinput = Tinput;

  return calcRSSpectrum(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput,
			Dem, qtherm, velocity, fluxArray, fluxErrArray);
}

int calcRSSpectrum(const RealArray& energyArray, 
		   const IntegerVector& Zinput, const RealArray& abundance,
		   const Real Redshift, const Real& Tinput, 
		   const Real& Tbinput,
		   const Real& Dem, const bool qtherm, const Real velocity,
		   RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray TbinputArr(1);
  RealArray DemArr(1);
  TinputArr[0] = Tinput;
  TbinputArr[0] = Tbinput;
  DemArr[0] = Dem;

  return calcRSSpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
			TbinputArr, DemArr, qtherm, velocity, fluxArray, fluxErrArray);
}


int calcRSSpectrum(const RealArray& energyArray, 
		   const IntegerVector& Zinput, const RealArray& abundance,
		   const Real Redshift, const RealArray& Tinput, 
		   const RealArray& Tbinput,
		   const RealArray& Dem, const bool qtherm, 
		   const Real velocity,
		   RealArray& fluxArray, RealArray& fluxErrArray)
{
  static bool isFirst(true);
  static Aped AtomicData;

  if ( isFirst ) {
    const string& datadir = FunctionUtility::modelDataPath();
    string contfil = datadir + "RS93_coco.fits";
    string linefil = datadir + "RS93_line.fits";

    int status = AtomicData.Read(contfil, linefil);
    if ( status != 0 ) return status;
    isFirst = false;
  }

  AtomicData.SetLineSpecs(false, false, qtherm, velocity, 0.0, false, false);

  // calculate the spectrum

  AtomicData.SumEqSpectra(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput, Dem,
			  fluxArray, fluxErrArray);

  return 0;
}

int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const vector<RealArray>& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  vector<vector<RealArray> > IonFracArr(1);
  TinputArr[0] = Tinput;
  IonFracArr[0].resize(IonFrac.size());
  for (size_t i=0; i<IonFrac.size(); i++) {
    IonFracArr[0][i].resize(IonFrac[i].size());
    for (size_t j=0; j<IonFrac[i].size(); j++) {
      IonFracArr[0][i][j] = IonFrac[i][j];
    }
  }

  return calcNEISpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
			 IonFracArr, qtherm, velocity, fluxArray, fluxErrArray);
}

int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const vector<vector<RealArray> >& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray Tbinput(Tinput.size());
  Tbinput = Tinput;

  return calcNEISpectrum(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput,
			 IonFrac, qtherm, velocity, fluxArray, fluxErrArray);
}

int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const Real& Tinput,
		    const Real& Tbinput,
		    const vector<RealArray>& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  RealArray TinputArr(1);
  RealArray TbinputArr(1);
  vector<vector<RealArray> > IonFracArr(1);
  TinputArr[0] = Tinput;
  TbinputArr[0] = Tbinput;
  IonFracArr[0].resize(IonFrac.size());
  for (size_t i=0; i<IonFrac.size(); i++) {
    IonFracArr[0][i].resize(IonFrac[i].size());
    for (size_t j=0; j<IonFrac[i].size(); j++) {
      IonFracArr[0][i][j] = IonFrac[i][j];
    }
  }

  return calcNEISpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
			 TbinputArr, IonFracArr, qtherm, velocity, fluxArray, 
			 fluxErrArray);
}

// prototype for old Kazik code from KBcalcNEISpectrum.cxx
int KBcalcNEISpectrum(const RealArray& energyArray, const IntegerVector& Zinput,
		      const RealArray& abundance, const Real Redshift,
		      const RealArray& Tinput, 
		      const vector<vector<RealArray> >& IonFrac,
		      const bool qtherm, const Real velocity,
		      RealArray& fluxArray, RealArray& fluxErrArray);


int calcNEISpectrum(const RealArray& energyArray, 
		    const IntegerVector& Zinput, const RealArray& abundance,
		    const Real Redshift, const RealArray& Tinput,
		    const RealArray& Tbinput,
		    const vector<vector<RealArray> >& IonFrac, 
		    const bool qtherm, const Real velocity,
		    RealArray& fluxArray, RealArray& fluxErrArray)
{
  static Aped AtomicData;
  static string contfil, linefil;
  static bool isFirst(true);


  // check for NEIVERS and NEIAPECROOT having been xset and return the
  // filenames to get the spectral data.

  static string version = FunctionUtility::atomdbVersion();
  static string eigenVersion = FunctionUtility::neiVersion();
  bool sameVersion = getNEIApedFileNames(version, eigenVersion, contfil, linefil);

  // If the eigenfunction version number is 1.x then it is not Aped so call the old
  // code

  if ( eigenVersion.substr(0,1) == "1" ) {
         int status = KBcalcNEISpectrum(energyArray, Zinput, abundance, Redshift, 
					Tinput, IonFrac, qtherm, velocity, 
					fluxArray, fluxErrArray);
     if ( status != 0 ) {
       std::stringstream oss;
       oss << "KBcalcNEISpectrum failed: status = " << status << "\n";
       FunctionUtility::xsWrite(oss.str(),10);       
     }
     return status;
  }

  if ( !sameVersion || isFirst ) {

    std::ostringstream oss;
    oss << "Reading NEI APEC spectral data from " << version << "\n";
    oss << "and eigenfunction data from " << eigenVersion << "\n";
    FunctionUtility::xsWrite(oss.str(),10);

    std::ostringstream oss2;
    oss2 << "Reading continuum data from " << contfil << "\n"
	 << "Reading line      data from " << linefil << "\n";
    FunctionUtility::xsWrite(oss2.str(),25);

    int status = AtomicData.Read(contfil, linefil);
    if ( status != 0 ) return status;
    isFirst = false;
  }

  // Check for values of APECTHERMAL, APECVELOCITY, APECMINFLUX, 
  // APECBROADPSEUDO, APECNOLINES, APECNORES, and APECLOGTINTERP xset 
  // parameters

  AtomicData.SetLineSpecs(false, false, qtherm, velocity, 0.0, false, false);

  std::ostringstream oss;
  oss << "Ion Fractions used:" << "\n";
  for (size_t it=0; it<Tinput.size(); it++) {
    oss << "Temperature = " << Tinput[it] << " keV" << "\n";
    for (size_t i=0; i<Zinput.size(); i++) {
      oss << Zinput[i] << ": ";
      for (size_t j=0; j<(size_t)(Zinput[i]+1); j++) {
	oss << IonFrac[it][i][j] << " ";
      }
      oss << "\n";
    }
  }
  FunctionUtility::xsWrite(oss.str(), 25);

  // calculate the spectrum

  AtomicData.SumNeqSpectra(energyArray, Zinput, abundance, Redshift, Tinput, Tbinput,
			   IonFrac, fluxArray, fluxErrArray);

  return 0;
}

