// Code for IonBalNei, IonBalTemperatureRecord, and IonBalElementRecord classes


#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSstreams.h>
#include <sstream>
#include <cmath>
#include <CCfits/CCfits>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <memory>

#include "IonBalNei.h"

using namespace CCfits;

// Methods for IonBalNei class.

// default constructor

IonBalNei::IonBalNei()
{
  TemperatureRecord.resize(0);
  Version = FunctionUtility::neiVersion();
}

// destructor

IonBalNei::~IonBalNei()
{
  TemperatureRecord.resize(0);
}

//  Set and get the version string. Note that if the version is changed
//  we call Clear to reset the object and setVersion returns true.
//  The setVersion method with no input checks NEIVERS and uses that if
//  set.

string IonBalNei::getVersion()
{
  return Version;
}

bool IonBalNei::setVersion()
{
  bool change = false;

  const string pname("NEIVERS");
  string version(FunctionUtility::getModelString(pname));
  if ( version.length() && version != FunctionUtility::NOT_A_KEY() ) 
    change = setVersion(version);

  return change;
}

bool IonBalNei::setVersion(string inVersion)
{

  string version = inVersion;

  // check for a valid version string

  if ( version != "1.0" && version != "1.1" && 
       version != "2.0" && version.substr(0,2) != "3." ) {
    std::ostringstream oss;
    oss << version << " is not a valid NEIVERS, ignoring." << "\n";
    xs_write(const_cast<char*>(oss.str().c_str()),5);    
    return false;
  }

  bool change = false;
  if ( version != Version ) change = true;
  if ( change ) {
    std::ostringstream oss;
    oss << "Resetting internal eigenvector data since version has changed from " 
	  << Version << " to " << inVersion << "\n";
    xs_write(const_cast<char*>(oss.str().c_str()),10);
    Version = inVersion;
    Clear();
  }
  return change;
}

int IonBalNei::ReadElements(int Z, string dirname)
{
  vector<int> Zarray(1);
  Zarray[0] = Z;
  int status = ReadElements(Zarray, dirname);
  return status;
}

// read from file for the requested elements

int IonBalNei::ReadElements(vector<int> Zarr, string dirname) 
{

  if ( Zarr.size() == 0 ) return(0);

  // construct the version name used in the files based on the version in use
  string versionname;
  if ( Version.substr(0,2) == "3." ) {
    versionname = "_v"+Version;
  } else if ( Version.substr(0,3) == "1.1" || Version.substr(0,3) == "2.0") {
    versionname = "0502";
  } else if ( Version.substr(0,3) == "1.0" ) {
    versionname = "";
  } else {
    std::ostringstream oss;
    oss << Version << " is not a valid NEIVERS version number." << "\n";
    xs_write(const_cast<char*>(oss.str().c_str()),15);
    return -1;
  }

  string filename = dirname + "eigen" + versionname + ".fits";
  vector<string> hduName(Zarr.size());
  for (size_t i=0; i<Zarr.size(); i++) hduName[i] = atomName[Zarr[i]-1];

  // open file
  std::unique_ptr<FITS> pIonBalfile((FITS*)0);

  try {
    pIonBalfile.reset(new FITS(filename,CCfits::Read));
  } catch(...) {
    FunctionUtility::xsWrite("Failed to read eigenvector data from "+filename,5);
    return(1);
  }

  FunctionUtility::xsWrite("Reading eigenvector data from "+filename,15);

  // loop round elements

  for (size_t iZ=0; iZ<Zarr.size(); iZ++) {

    int Z = Zarr[iZ];

    try {
      ExtHDU& ionbalExt = pIonBalfile->extension(hduName[iZ]);

      // get the number of temperatures from the number of rows in the 
      // PARAMETERS extension and check it is the number we are expecting

      size_t NumberTemperatures = (size_t)ionbalExt.rows();
      if ( NumberTemperatures != NUMBTEMP ) {
	FunctionUtility::xsWrite("Unexpected number of temperatures in "+filename+"["+hduName[iZ]+"]",5);
	return(2);
      }

      // if the temperature records do not exist then set them up, if they do
      // exist then check for consistency

      if ( TemperatureRecord.size() == 0 ) {
	TemperatureRecord.resize(NumberTemperatures);
	Temperature.resize(NumberTemperatures);
      } else {
	if ( TemperatureRecord.size() != NumberTemperatures ) {
	  FunctionUtility::xsWrite("Unexpected number of temperatures in "+filename+"["+hduName[iZ]+"]",5);
	  return(2);
	}
      }

      // read all the eigenvector data from the file. More efficient to do this
      // in one go then copy into objects

      vector<RealArray> FEQB(NumberTemperatures);
      vector<RealArray> EIG(NumberTemperatures);
      vector<RealArray> VR(NumberTemperatures);
      vector<RealArray> VL(NumberTemperatures);
      ionbalExt.column("FEQB").readArrays(FEQB, (long)1, (long)NumberTemperatures);
      ionbalExt.column("EIG").readArrays(EIG, (long)1, (long)NumberTemperatures);
      ionbalExt.column("VR").readArrays(VR, (long)1, (long)NumberTemperatures);
      ionbalExt.column("VL").readArrays(VL, (long)1, (long)NumberTemperatures);

      // loop round temperatures

      for (size_t iTemp=0; iTemp<NumberTemperatures; iTemp++) {

	// set up temperature record and set temperature. Note that for convenience
	// we store the temperature both in the record and in a separate array in
	// the IonBalNei object.

	IonBalTemperatureRecord& Trecord = TemperatureRecord[iTemp];

	Trecord.Temperature = pow(10.0,MINLOGT+iTemp*DELTALOGT)/KEVTOK;
	Temperature[iTemp] = Trecord.Temperature;

	// check whether the element record exists for this temperature

	int EltIndex = -1;
	for (size_t i=0; i<Trecord.ElementRecord.size(); i++) {
	  if (Trecord.ElementRecord[i].AtomicNumber == Z ) {
	    EltIndex = i;
	  }
	}

	if ( EltIndex == -1 ) {
	  IonBalElementRecord Erecord;
	  Erecord.AtomicNumber = Z;
	  Trecord.LoadElementRecord(Erecord);
	}
	EltIndex = Trecord.ElementRecord.size()-1;

	IonBalElementRecord& Erecord = Trecord.ElementRecord[EltIndex];

	// Load the eigenvector data into the object arrays

	Erecord.EquilibriumPopulation.resize(Z+1);
	Erecord.EquilibriumPopulation = FEQB[iTemp];
	Erecord.Eigenvalues.resize(Z);
	Erecord.Eigenvalues = EIG[iTemp];
	Erecord.RightEigenvectors.resize(Z);
	for (size_t i=0; i<(size_t)Z; i++) {
	  Erecord.RightEigenvectors[i].resize(Z);
	  for (size_t j=0; j<(size_t)Z; j++) {
	    Erecord.RightEigenvectors[i][j] = VR[iTemp][i*Z+j];
	  }
	}
	Erecord.LeftEigenvectors.resize(Z);
	for (size_t i=0; i<(size_t)Z; i++) {
	  Erecord.LeftEigenvectors[i].resize(Z);
	  for (size_t j=0; j<(size_t)Z; j++) {
	    Erecord.LeftEigenvectors[i][j] = VL[iTemp][i*Z+j];
	  }
	}
	
	// end loop over temperatures

      }

    } catch(...) {
      FunctionUtility::xsWrite("Failed to open "+filename+"["+hduName[iZ]+"]",5);
      return(1);
    }

    // end loop over elements

  }

  return(0);

}

// load a temperature record

void IonBalNei::LoadTemperatureRecord(IonBalTemperatureRecord input)
{
  TemperatureRecord.push_back(input);
  return;
}

// return the temperatures stored

RealArray IonBalNei::Temperatures()
{
  RealArray Temperature(TemperatureRecord.size());
  for (size_t i=0; i<TemperatureRecord.size(); i++) {
    Temperature[i] = TemperatureRecord[i].Temperature;
  }
  return Temperature;
}

// return the number of temperatures

int IonBalNei::NumberTemperatures()
{
  return (int)TemperatureRecord.size();
}

// return the number of elements - assumes that all elements are included
// for all temperatures

int IonBalNei::NumberElements()
{
  return (int)TemperatureRecord[0].ElementRecord.size();
}

// return true if the requested element is included - assumes that all elements
// are included for all temperatures

bool IonBalNei::ContainsElement(const int& Z)
{
  bool found(false);
  if ( TemperatureRecord.size() > 0 ) {
    for (size_t i=0; i<TemperatureRecord[0].ElementRecord.size(); i++) {
      if ( TemperatureRecord[0].ElementRecord[i].AtomicNumber == Z ) found = true;
    }
  }
  return found;
}

// deep copy

IonBalNei& IonBalNei::operator=(const IonBalNei& beta)
{
  TemperatureRecord.resize(beta.TemperatureRecord.size());
  for (size_t i=0; i<TemperatureRecord.size(); i++) {
    TemperatureRecord[i] = beta.TemperatureRecord[i];
  }
  return *this;
}

// clear out the contents of the object

void IonBalNei::Clear()
{
  for (size_t i=0; i<TemperatureRecord.size(); i++) {
    TemperatureRecord[i].Clear();
  }
  TemperatureRecord.resize(0);
  return;
}

// return the collisional equilibrium ionization fractions for electron
// temperature Te and element Z

RealArray IonBalNei::CEI(const Real& Te, const int& Z)
{
  RealArray IonFrac(Z+1);
  IonFrac = 0.0;

  // Find the the tabulated temperatures which bracket Te using a binary
  // search

  int iTlow = locateIndex(Temperature, Te);
  int iThigh;
  if ( iTlow == (int)Temperature.size()-1 ) {
    iThigh = iTlow;
  } else {
    iThigh = iTlow + 1;
  }

  IonBalTemperatureRecord& TrecordLow = TemperatureRecord[iTlow];
  IonBalTemperatureRecord& TrecordHigh = TemperatureRecord[iThigh];

  // Now find the ElementRecord which corresponds to Z

  int EltIndex = -1;
  for (size_t i=0; i<TrecordLow.ElementRecord.size(); i++) {
    if ( TrecordLow.ElementRecord[i].AtomicNumber == Z ) {
      EltIndex = i;
    }
  }

  if ( EltIndex == -1 ) return IonFrac;

  IonBalElementRecord& ErecordLow = TrecordLow.ElementRecord[EltIndex];  
  IonBalElementRecord& ErecordHigh = TrecordHigh.ElementRecord[EltIndex];  

  Real Tdiff = TrecordHigh.Temperature - TrecordLow.Temperature;
  if ( Tdiff > 0.0 ) {
    Real factorLow = (TrecordHigh.Temperature-Te)/Tdiff;
    Real factorHigh = (Te-TrecordLow.Temperature)/Tdiff;
    IonFrac = factorLow * ErecordLow.EquilibriumPopulation +
      factorHigh * ErecordHigh.EquilibriumPopulation;
  } else {
    IonFrac = ErecordLow.EquilibriumPopulation;
  }

  return IonFrac;

}



// calculate the NEI ion fractions for electron temperature Te, ionization
// nt tau and element Z.

RealArray IonBalNei::Calc(const Real& Te, const Real& tau, const int& Z)
{
  RealArray initIonFrac(Z+1);
  initIonFrac = 0.0;
  initIonFrac[0] = 1.0;

  return this->Calc(Te, tau, Z, initIonFrac);
}

RealArray IonBalNei::Calc(const Real& Te, const Real& tau, const int& Z, 
			  const RealArray& initIonFrac)
{
  RealArray TeArr(1), tauArr(1), weight(1);
  TeArr[0] = Te;
  tauArr[0] = tau;
  weight[0] = 1.0;

  return this->Calc(TeArr, tauArr, weight, Z, initIonFrac);
}

// Calculates ionization fractions at electron temperatures
// Te and a set of ionization parameters tau(i), i=1,..,n,
// where the contribution of each tau is weight(i). Electron 
// temperature is assumed to be linear function of tau.

RealArray IonBalNei::Calc(const RealArray& Te, const RealArray& tau, 
			  const RealArray& weight, const int& Z)
{
  RealArray initIonFrac(Z+1);
  initIonFrac = 0.0;
  initIonFrac[0] = 1.0;

  return this->Calc(Te, tau, weight, Z, initIonFrac);
}

RealArray IonBalNei::Calc(const RealArray& Te, const RealArray& tau, 
			  const RealArray& weight, const int& Z, const RealArray& initIonFrac)
{
  RealArray IonFrac(Z+1);
  IonFrac = 0.0;

  // special case if Z=1 or 2 which we assume is always completely ionized
  if ( Z == 1 ) {
    IonFrac[0] = 0.0;
    IonFrac[1] = 1.0;
    return IonFrac;
  } else if ( Z == 2 ) {
    IonFrac[0] = 0.0;
    IonFrac[1] = 0.0;
    IonFrac[2] = 1.0;
    return IonFrac;
  }
    
  // If tau and weight have different sizes then give up

  if ( tau.size() != weight.size() ) return IonFrac;

  // As a first approximation find the tabulated temperature closest to Te
  // for each Te and store information in Tindex array

  IntegerVector Tindex(Te.size());
  for (size_t it=0; it<Tindex.size(); it++) {
    Real Tdiff = 1e10;
    for (size_t i=0; i<TemperatureRecord.size(); i++) {
      if ( fabs(Te[it]-TemperatureRecord[i].Temperature) < Tdiff ) {
	Tindex[it] = i;
	Tdiff = fabs(Te[it]-TemperatureRecord[i].Temperature);
      }
    }
  }

  Real del = 1.0/((Real)(Te.size()-1));

  // loop over tau values

  for (size_t itau=0; itau<tau.size(); itau++) {

    // Starting ionization

    RealArray fs(Z+1);
    fs = initIonFrac;

    // loop over temperatures

    for (size_t it=0; it<Te.size()-1; it++) {

      IonBalTemperatureRecord& Trecord = TemperatureRecord[Tindex[it]];

      // Now find the ElementRecord which corresponds to Z

      int EltIndex = -1;
      for (size_t i=0; i<Trecord.ElementRecord.size(); i++) {
	if ( Trecord.ElementRecord[i].AtomicNumber == Z ) {
	  EltIndex = i;
	}
      }

      if ( EltIndex == -1 ) return IonFrac;

      IonBalElementRecord& Erecord = Trecord.ElementRecord[EltIndex];

      // Decompose into eigenvectors

      RealArray work(Z);
      for (size_t i=0; i<(size_t)Z; i++) {
	work[i] = fs[i+1] - Erecord.EquilibriumPopulation[i+1];
      }
      RealArray fspec(Z);
      fspec = 0.0;
      for (size_t i=0; i<(size_t)Z; i++) {
	for (size_t j=0; j<(size_t)Z; j++) {
	  fspec[i] += Erecord.LeftEigenvectors[i][j] * work[j];
	}
      }

      // propagate forward

      for (size_t i=0; i<(size_t)Z; i++) {
	work[i] = fspec[i]*exp(Erecord.Eigenvalues[i] * del * tau[itau]);
      }

      Real sum(0.0);
      for (size_t i=0; i<(size_t)Z; i++) {
	fs[i+1] = 0.0;
	for (size_t j=0; j<(size_t)Z; j++) {
	  fs[i+1] += work[j]*Erecord.RightEigenvectors[j][i];
	}
	fs[i+1] += Erecord.EquilibriumPopulation[i+1];
	if ( fs[i+1] < 0.0 ) fs[i+1] = 0.0;
	sum += fs[i+1];
      }
      fs[0] = 1.0 - sum;
      if ( fs[0] < 0.0 ) fs[0] = 0.0;

      // end loop over temperatures
    }

    Real fract = weight[itau];

    for (size_t i=0; i<IonFrac.size(); i++) {
      IonFrac[i] += fract * fs[i];
    }

    // end loop over tau

  }

  return IonFrac;
}

// Calculates ionization fractions at electron temperature
// Te(n) and ionization parameter tau(n), for electron 
// temperatures Te given in a tabular form as a function of 
// ionization parameter tau. 

RealArray IonBalNei::Calc(const RealArray& Te, const RealArray& tau, 
			  const int& Z)
{
  RealArray initIonFrac(Z+1);
  initIonFrac = 0.0;
  initIonFrac[0] = 1.0;

  return this->Calc(Te, tau, Z, initIonFrac);
}

RealArray IonBalNei::Calc(const RealArray& Te, const RealArray& tau, 
			  const int& Z, const RealArray& initIonFrac)
{
  // Starting ionization

  RealArray IonFrac(Z+1);
  IonFrac = initIonFrac;

  // special case if Z=1 or 2 which we assume is always completely ionized
  if ( Z == 1 ) {
    IonFrac[0] = 0.0;
    IonFrac[1] = 1.0;
    return IonFrac;
  } else if ( Z == 2 ) {
    IonFrac[0] = 0.0;
    IonFrac[1] = 0.0;
    IonFrac[2] = 1.0;
    return IonFrac;
  }
    
  // As a first approximation find the tabulated temperature closest to Te
  // for each Te and store information in Tindex array

  IntegerVector Tindex(Te.size());
  for (size_t it=0; it<Tindex.size(); it++) {
    Real Tdiff = 1e10;
    for (size_t i=0; i<TemperatureRecord.size(); i++) {
      if ( fabs(Te[it]-TemperatureRecord[i].Temperature) < Tdiff ) {
	Tindex[it] = i;
	Tdiff = fabs(Te[it]-TemperatureRecord[i].Temperature);
      }
    }
  }

  // loop over temperatures

  for (size_t it=0; it<Te.size()-1; it++) {

    IonBalTemperatureRecord& Trecord = TemperatureRecord[Tindex[it]];

    // Now find the ElementRecord which corresponds to Z

    int EltIndex = -1;
    for (size_t i=0; i<Trecord.ElementRecord.size(); i++) {
      if ( Trecord.ElementRecord[i].AtomicNumber == Z ) {
	EltIndex = i;
      }
    }
    if ( EltIndex == -1 ) return IonFrac;

    IonBalElementRecord& Erecord = Trecord.ElementRecord[EltIndex];

    // Decompose into eigenvectors

    RealArray work(Z);
    for (size_t i=0; i<(size_t)Z; i++) {
      work[i] = IonFrac[i+1] - Erecord.EquilibriumPopulation[i+1];
    }
    RealArray fspec(Z);
    fspec = 0.0;
    for (size_t i=0; i<(size_t)Z; i++) {
      for (size_t j=0; j<(size_t)Z; j++) {
	fspec[i] += Erecord.LeftEigenvectors[i][j] * work[j];
      }
    }
    
    // propagate forward

    for (size_t i=0; i<(size_t)Z; i++) {
      work[i] = fspec[i]*exp(Erecord.Eigenvalues[i] * (tau[it+1]-tau[it]));
    }

    Real sum(0.0);
    for (size_t i=0; i<(size_t)Z; i++) {
      IonFrac[i+1] = 0.0;
      for (size_t j=0; j<(size_t)Z; j++) {
	IonFrac[i+1] += work[j]*Erecord.RightEigenvectors[j][i];
      }
      IonFrac[i+1] += Erecord.EquilibriumPopulation[i+1];
      if ( IonFrac[i+1] < 0.0 ) IonFrac[i+1] = 0.0;
      sum += IonFrac[i+1];
    }
    IonFrac[0] = 1.0 - sum;
    if ( IonFrac[0] < 0.0 ) IonFrac[0] = 0.0;

    // end loop over temperatures
  }

  return IonFrac;
}

// methods for IonBalTemperatureRecord

IonBalTemperatureRecord::IonBalTemperatureRecord()
{
}

IonBalTemperatureRecord::~IonBalTemperatureRecord()
{
}

void IonBalTemperatureRecord::LoadElementRecord(IonBalElementRecord input)
{
  ElementRecord.push_back(input);
  return;
}

// clear out the contents of the object

void IonBalTemperatureRecord::Clear()
{
  for (size_t i=0; i<ElementRecord.size(); i++) {
    ElementRecord[i].Clear();
  }
  ElementRecord.resize(0);
  return;
}

// deep copy

IonBalTemperatureRecord& IonBalTemperatureRecord::operator=(const IonBalTemperatureRecord& beta)
{
  Temperature = beta.Temperature;
  ElementRecord.resize(beta.ElementRecord.size());
  for (size_t i=0; i<ElementRecord.size(); i++) {
    ElementRecord[i] = beta.ElementRecord[i];
  }
  return *this;
}

// methods for IonBalElementRecord

IonBalElementRecord::IonBalElementRecord()
{
}

IonBalElementRecord::~IonBalElementRecord()
{
}

// clear out the contents of the object

void IonBalElementRecord::Clear()
{
  AtomicNumber = 0;
  EquilibriumPopulation.resize(0);
  Eigenvalues.resize(0);
  for (size_t i=0; i<LeftEigenvectors.size(); i++) {
    LeftEigenvectors[i].resize(0);
    RightEigenvectors[i].resize(0);
  }
  LeftEigenvectors.resize(0);
  RightEigenvectors.resize(0);

  return;
}

// deep copy

IonBalElementRecord& IonBalElementRecord::operator=(const IonBalElementRecord& beta)
{
  AtomicNumber = beta.AtomicNumber;
  EquilibriumPopulation.resize(beta.EquilibriumPopulation.size());
  EquilibriumPopulation = beta.EquilibriumPopulation;
  Eigenvalues.resize(beta.Eigenvalues.size());
  Eigenvalues = beta.Eigenvalues;
  LeftEigenvectors.resize(beta.LeftEigenvectors.size());
  for (size_t i=0; i<LeftEigenvectors.size(); i++) {
    LeftEigenvectors[i].resize(beta.LeftEigenvectors[i].size());
    LeftEigenvectors[i] = beta.LeftEigenvectors[i];
  }
  RightEigenvectors.resize(beta.RightEigenvectors.size());
  for (size_t i=0; i<RightEigenvectors.size(); i++) {
    RightEigenvectors[i].resize(beta.RightEigenvectors[i].size());
    RightEigenvectors[i] = beta.RightEigenvectors[i];
  }
  return *this;
}

// Functions outside the class but wrapping up loading and calc methods. Try to
// be clever here and funnel everything through the most general case

void calcNEIfractions(const Real& Te, const Real& tau, const int& Z, RealArray& IonFrac)
{
  RealArray initIonFrac(0);
  calcNEIfractions(Te, tau, Z, initIonFrac, IonFrac);
  return;
}
void calcNEIfractions(const Real& Te, const Real& tau, const IntegerVector& Z, vector<RealArray>& IonFrac)
{
  vector<RealArray> initIonFrac(0);
  calcNEIfractions(Te, tau, Z, initIonFrac, IonFrac);
  return;
}

void calcNEIfractions(const Real& Te, const Real& tau, const int& Z, 
			   const RealArray& initIonFrac, RealArray& IonFrac)
{
  RealArray TeArray(2), tauArray(1), weight(1);
  TeArray[0] = Te;
  TeArray[1] = Te;
  tauArray[0] = tau;
  weight[0] = 1.0;

  calcNEIfractions(TeArray, tauArray, weight, Z, initIonFrac, IonFrac);

  return;
}

void calcNEIfractions(const Real& Te, const Real& tau, const IntegerVector& Z, 
		      const vector<RealArray>& initIonFrac, vector<RealArray>& IonFrac)
{
  RealArray TeArray(2), tauArray(1), weight(1);
  TeArray[0] = Te;
  TeArray[1] = Te;
  tauArray[0] = tau;
  weight[0] = 1.0;

  calcNEIfractions(TeArray, tauArray, weight, Z, initIonFrac, IonFrac);
  return;
}

void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const int& Z, RealArray& IonFrac)
{
  RealArray weight(0), initIonFrac(0);

  calcNEIfractions(Te, tau, weight, Z, initIonFrac, IonFrac);
  return;
}
void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const IntegerVector& Z, vector<RealArray>& IonFrac)
{
  vector<RealArray> initIonFrac(0);
  RealArray weight(0);

  calcNEIfractions(Te, tau, weight, Z, initIonFrac, IonFrac);
  return;
}

void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const int& Z, const RealArray& initIonFrac, RealArray& IonFrac)
{
  RealArray weight(0);

  calcNEIfractions(Te, tau, weight, Z, initIonFrac, IonFrac);
  return; 
}
void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const IntegerVector& Z, const vector<RealArray>& initIonFrac,
		      vector<RealArray>& IonFrac)
{
  RealArray weight(0);

  calcNEIfractions(Te, tau, weight, Z, initIonFrac, IonFrac);
  return; 
}

void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const RealArray& weight, const int& Z, RealArray& IonFrac)
{
  RealArray initIonFrac(0);
  calcNEIfractions(Te, tau, weight, Z, initIonFrac, IonFrac);
  return; 
}

void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const RealArray& weight, const IntegerVector& Z, vector<RealArray>& IonFrac)
{
  vector<RealArray> initIonFrac(0);
  calcNEIfractions(Te, tau, weight, Z, initIonFrac, IonFrac);
  return; 
}

void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const RealArray& weight, const int& Z, 
		      const RealArray& initIonFrac, RealArray& IonFrac)
{

  IntegerVector Zarray(1);
  Zarray[0] = Z;
  vector<RealArray> tinitIonFrac(1);
  tinitIonFrac[0].resize(initIonFrac.size());
  tinitIonFrac[0] = initIonFrac;

  vector<RealArray> tIonFrac;
  calcNEIfractions(Te, tau, weight, Zarray, tinitIonFrac, tIonFrac);

  IonFrac.resize(tIonFrac[0].size());
  IonFrac = tIonFrac[0];

  return;
}

// this is the general case. only recalculates IonFrac if necessary

void calcNEIfractions(const RealArray& Te, const RealArray& tau, 
		      const RealArray& weight, const IntegerVector& Zarray, 
		      const vector<RealArray>& initIonFrac, 
		      vector<RealArray>& IonFrac)
{

  static IonBalNei NeiData;
  static RealArray saveTe;
  static RealArray saveTau;
  static RealArray saveWeight;
  static IntegerVector saveZarray;
  static vector<RealArray> saveInitIonFrac;
  static vector<RealArray> saveIonFrac;

  // it will be useful to have a boolean on whether initIonFrac has been set
  bool noInit = true;
  if ( initIonFrac.size() > 0 ) noInit = false;

  // update the version if necessary

  bool changed = NeiData.setVersion();

  // check for a change in the input arrays

  if ( !identicalArrays(Te,saveTe) || !identicalArrays(tau,saveTau) ||
       !identicalArrays(weight, saveWeight) || 
       !identicalArrays(Zarray, saveZarray) ||
       !identicalArrays(initIonFrac, saveInitIonFrac) ) {
    changed = true;
    saveTe.resize(Te.size());
    saveTe = Te;
    saveTau.resize(tau.size());
    saveTau = tau;
    saveZarray.resize(Zarray.size());
    saveZarray = Zarray;
    saveWeight.resize(weight.size());
    saveWeight = weight;
    saveInitIonFrac.resize(initIonFrac.size());
    for (size_t iZ=0; iZ<initIonFrac.size(); iZ++) {
      saveInitIonFrac[iZ].resize(initIonFrac[iZ].size());
      saveInitIonFrac[iZ] = initIonFrac[iZ];
    }
  }

  // if nothing has changed then just return the saved values from the
  // previous run

  if ( !changed ) {
    IonFrac.resize(saveIonFrac.size());
    for (size_t iZ=0; iZ<IonFrac.size(); iZ++) {
      IonFrac[iZ].resize(saveIonFrac[iZ].size());
      IonFrac[iZ] = saveIonFrac[iZ];
    }
    return;
  }

  // need to calculate so loop over elements
  // first we need to find out which ones require data to be read and

  vector<int> ZreadList;
  for (size_t iZ=0; iZ<Zarray.size(); iZ++) {
    int Z = Zarray[iZ];
    if ( Z > 2 ) {
      if ( !NeiData.ContainsElement(Z) ) {
	bool fileShouldExist;
	if ( NeiData.Version.substr(0,1) == "1" || 
	     NeiData.Version.substr(0,1) == "2" ) {
	  fileShouldExist = false;
	  for (size_t i=0; i<NoldVersionZ; i++) {
	    if ( Z == oldVersionZ[i] ) fileShouldExist = true;
	  }
	} else {
	  fileShouldExist = true;
	}
	if ( fileShouldExist ) ZreadList.push_back(Z);
      }
    }
  }

  // read in the data requested.

  if ( ZreadList.size() > 0 ) {
    const string& datadir = FunctionUtility::modelDataPath();
    NeiData.ReadElements(ZreadList, datadir);
  }

  // calculate ion fractions

  IonFrac.resize(Zarray.size());
  saveIonFrac.resize(Zarray.size());
  for (size_t iZ=0; iZ<Zarray.size(); iZ++) {
    int Z = Zarray[iZ];

  // first do the special case for Z=1 or 2, which are assumed to be
  // fully ionized

    if ( Z == 1 ) {
      IonFrac[iZ].resize(2);
      IonFrac[iZ][0] = 0.0;
      IonFrac[iZ][1] = 1.0;
    } else if ( Z == 2 ) {
      IonFrac[iZ].resize(3);
      IonFrac[iZ][0] = 0.0;
      IonFrac[iZ][1] = 0.0;
      IonFrac[iZ][2] = 1.0;
    } else {

      // now can go onto the other possibilities

      // Resize and initialize the output array
      
      IonFrac[iZ].resize(Z+1);
      IonFrac[iZ] = 0.0;

      // if something has gone wrong and data for this element is not available
      // then jump over it

      if ( NeiData.ContainsElement(Z) ) {

	// first special case of negative tau which means want 
	// the equilibrium fractions

	if ( tau[0] < 0.0 ) {
	  IonFrac[iZ] = NeiData.CEI(Te[0], Z);
	} else {

	  // now apply the appropriate Calc method

	  if ( Te.size() == 1 && tau.size() == 1 ) {
	    if ( noInit ) {
	      IonFrac[iZ] = NeiData.Calc(Te[0], tau[0], Z);
	    } else {
	      IonFrac[iZ] = NeiData.Calc(Te[0], tau[0], Z, initIonFrac[iZ]);
	    }
	  } else if ( weight.size() == 0 ) {
	    if ( noInit ) {
	      IonFrac[iZ] = NeiData.Calc(Te, tau, Z);
	    } else {
	      IonFrac[iZ] = NeiData.Calc(Te, tau, Z, initIonFrac[iZ]);
	    }
	  } else {
	    if ( noInit ) {
	      IonFrac[iZ] = NeiData.Calc(Te, tau, weight, Z);
	    } else {
	      IonFrac[iZ] = NeiData.Calc(Te, tau, weight, Z, initIonFrac[iZ]);
	    }
	  }

	}
      }

    }

    // save these ion fractions in a static array so we don't have to
    // recalculate on next call if not necessary

    saveIonFrac[iZ].resize(IonFrac[iZ].size());
    saveIonFrac[iZ] = IonFrac[iZ];

  }

  return;

}


// routine to calculate the collisional equilibrium ion fractions. This calls
// the NEI routine with the special case of a negative tau.

void calcCEIfractions(const Real Te, const IntegerVector& Z, vector<RealArray>& IonFrac)
{
  Real tau(-1.0);
  calcNEIfractions(Te, tau, Z, IonFrac);
  return;
}

// helpful routine to return an index into the temperatures

int getNEItempIndex(const Real& tkeV)
{
  Real tK = tkeV * KEVTOK;
  int i = (int)((log10(tK)-MINLOGT)/DELTALOGT + 0.5);
  if ( i < 1 ) i = 1;
  if ( i > NUMBTEMP-1 ) i = NUMBTEMP-1;
  return i;
}

// and one to return the number of temperatures

int getNEInumbTemp()
{
  return NUMBTEMP;
}

// useful debugging routine

string writeIonFrac(const IntegerVector& Z, const vector<RealArray>& IonFrac)
{
  std::ostringstream oss;
  oss << "Ion Fractions:\n";
  for (size_t iZ=0; iZ<Z.size(); iZ++) {
    oss.width(2);
    oss << Z[iZ] << " : ";
    for (size_t i=0; i<IonFrac[iZ].size(); i++) {
      oss.precision(3);
      oss << IonFrac[iZ][i] << ", ";
    }
    oss << "\n";
  }
  return oss.str();
}

string writeIonFrac(const int& Z, const IntegerVector& Zarray, 
		    const vector<RealArray>& IonFrac)
{
  std::ostringstream oss;
  oss << "Ion Fractions for Z= ";
  oss.width(2);
  oss << Z << " : ";
  for (size_t iZ=0; iZ<Zarray.size(); iZ++) {
    if ( Zarray[iZ] == Z ) {
      for (size_t i=0; i<IonFrac[iZ].size(); i++) {
	oss.precision(3);
	oss << IonFrac[iZ][i] << ", ";
      }
      oss << "\n";
    }
  }
  return oss.str();
}


// helpful routine to do a binary search on a RealArray and return the index
// of the element of xx immediately less than x.

int locateIndex(const RealArray& xx, const Real x)
{
  size_t n = xx.size();
  
  size_t jl = -1;
  size_t ju = n;
  while ( ju - jl > 1 ) {
    size_t jm = ( ju + jl ) / 2;
    if ( (xx[n-1] > xx[0]) == (x > xx[jm]) ) {
      jl = jm;
    } else {
      ju = jm;
    }
  }
  
  if ( jl > 0 ) {
    return (int)jl;
  } else {
    return (int)0;
  }
  
}

// functions to check whether arrays are identical. Could do this with a
// template but not worth it at the moment. These routines should ideally
// be in utility routines. Better would be to set up RealArray and IntegerVector
// as classes with the comparison operators defined (note the C++11 standard
// does include these for arrays).

 bool identicalArrays(const vector<RealArray>& a, const vector<RealArray>& b)
{
  if ( a.size() != b.size() ) return false;
  for (size_t i=0; i<a.size(); i++) {
    if ( a[i].size() != b[i].size() ) return false;
  }
  for (size_t i=0; i<a.size(); i++) {
    for (size_t j=0; j<a[i].size(); j++) {
    // strictly should have this differ by greater than machine precision
      if ( a[i][j] != b[i][j] ) return false;
    }
  }
  return true;
}

bool identicalArrays(const RealArray& a, const RealArray& b)
{
  if ( a.size() != b.size() ) return false;
  for (size_t i=0; i<a.size(); i++) {
    // strictly should have this differ by greater than machine precision
    if ( a[i] != b[i] ) return false;
  }
  return true;
}

bool identicalArrays(const IntegerVector& a, const IntegerVector& b)
{
  if ( a.size() != b.size() ) return false;
  for (size_t i=0; i<a.size(); i++) {
    if ( a[i] != b[i] ) return false;
  }
  return true;
}
