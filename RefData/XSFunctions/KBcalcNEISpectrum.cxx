// A wrap-up for Kazik Borkowski's code to evaluate an NEI spectrum. This
// is an analogue for the calcNEISpectrum routine in Aped.h. Converts from
// the C++ arrays into C arrays to be passed to the Fortran neispec routine
// via cfortran.

#include <xsTypes.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSstreams.h>
#include <cfortran.h>

using namespace std;

PROTOCCALLSFSUB8(NEI1TSPEC,nei1tspec,FLOAT,INT,FLOATV,INT,FLOATV,
		 FLOAT,INT,FLOATV)
#define NEI1TSPEC(temp,nion,ionfrac,ne,ear,z,nelt,workphot) \
  CCALLSFSUB8(NEI1TSPEC,nei1tspec,FLOAT,INT,FLOATV,INT,FLOATV,FLOAT,	\
	      INT,FLOATV,temp,nion,ionfrac,ne,ear,z,nelt,workphot)

// prototypes

int KBcalcNEISpectrum(const RealArray& energyArray, const IntegerVector& Zinput,
		      const RealArray& abundance, const Real Redshift,
		      const Real& Tinput, 
		      const vector<RealArray>& IonFrac,
		      const bool qtherm, const Real velocity,
		      RealArray& fluxArray, RealArray& fluxErrArray);

int KBcalcNEISpectrum(const RealArray& energyArray, const IntegerVector& Zinput,
		      const RealArray& abundance, const Real Redshift,
		      const RealArray& Tinput, 
		      const vector<vector<RealArray> >& IonFrac,
		      const bool qtherm, const Real velocity,
		      RealArray& fluxArray, RealArray& fluxErrArray);

void KBconvertIonFrac(const vector<RealArray>& IonFrac, 
		      const IntegerVector& Zinput, float* ionf, int* ionel, 
		      int* ionstage);

// these are the elements which are included in the old NEI code.

const int nZincl = 12;
const int Zinclude[] = {1,2,6,7,8,10,12,14,16,20,26,28};
const int nIonp = 159;


// first the special case with a single temperature

int KBcalcNEISpectrum(const RealArray& energyArray, const IntegerVector& Zinput,
		      const RealArray& abundance, const Real Redshift,
		      const Real& Tinput, 
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

  return KBcalcNEISpectrum(energyArray, Zinput, abundance, Redshift, TinputArr,
			   IonFracArr, qtherm, velocity, fluxArray, fluxErrArray);
}

// now the general case with an array of temperatures

int KBcalcNEISpectrum(const RealArray& energyArray, const IntegerVector& Zinput,
		      const RealArray& abundance, const Real Redshift,
		      const RealArray& Tinput, 
		      const vector<vector<RealArray> >& IonFrac,
		      const bool qtherm, const Real velocity,
		      RealArray& fluxArray, RealArray& fluxErrArray)
{
  using namespace XSutility;

  static vector<RealArray> fluxSum(nZincl);

  int status = 0;


  // handy diagnostic info

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
  xs_write(const_cast<char*>(oss.str().c_str()),25);

  // set up the arrays required by the old-style code. note that we don't
  // actually need ionel and ionstage here but these were included in the
  // KBconvertIonFrac routine in case they were required elsewhere

  float* ionf = new float[nIonp];
  int* ionel = new int[nIonp];
  int* ionstage = new int[nIonp];

  // loop over temperatures

  int nE = (int)energyArray.size() - 1;
  float* workphot = new float[nE*nZincl];

  float* ear = new float[nE+1];
  for (size_t i=0; i<(size_t)(nE+1); i++) ear[i] = energyArray[i];

  float z = Redshift;

  for (size_t i=0; i<(size_t)nZincl; i++) {
    fluxSum[i].resize(nE);
    fluxSum[i] = 0.0;
  }

  for (size_t it=0; it<Tinput.size(); it++) {

    // convert IonFrac to ionf

    KBconvertIonFrac(IonFrac[it], Zinput, ionf, ionel, ionstage);

    NEI1TSPEC(Tinput[it], nIonp, ionf, nE, ear, z, nZincl, workphot);

    for (size_t i=0; i<(size_t)(nZincl); i++) {
      for (size_t j=0; j<(size_t)nE; j++) {
	fluxSum[i][j] += workphot[i*nE+j];
      }
    }

  }

  // Sum over elements, weighting by abundance

  fluxArray.resize(nE);
  for (size_t i=0; i<(size_t)(nZincl); i++) {
    Real abund = 0.0;
    for (size_t j=0; j<Zinput.size(); j++) {
      if ( Zinput[j] == Zinclude[i] ) abund = abundance[j];
    }
    Real weight = abund * FunctionUtility::getAbundance(Zinclude[i]);
    fluxArray += weight * fluxSum[i];
  }

  // Correct for redshift and normalize

  fluxArray *= 1.0e14 / (1.0 + z);
  fluxErrArray.resize(0);

  // tidy up memory

  delete[] ionf;
  delete[] ionel;
  delete[] ionstage;
  delete[] ear;
  delete[] workphot;

  return status;
}
//---------------------------------------------------------------------------------
// convert from the IonFrac structure into the old-style arrays used for NEI models
// this assumes that ionf, ionel, and ionstage have all been allocated enough
// memory.

void KBconvertIonFrac(const vector<RealArray>& IonFrac, 
		      const IntegerVector& Zinput, float* ionf, int* ionel, 
		      int* ionstage)
{
  // set up useful offsets into the ionf array

  IntegerVector Zoffset(nZincl);
  for (size_t i=0; i<3; i++) Zoffset[i] = i;
  for (size_t i=3; i<(size_t)nZincl; i++) {
    Zoffset[i] = Zoffset[i-1] + Zinclude[i-1] + 1;
  }

  // set ionel and ionstage arrays

  ionel[0] = 1;
  ionel[1] = 2;
  ionstage[0] = 2;
  ionstage[1] = 3;

  for (size_t i=2; i<(size_t)nZincl; i++) {
    for (size_t j=0; j<(size_t)(Zinclude[i]+1); j++) {
      ionel[Zoffset[i]+j] = Zinclude[i];
      ionstage[Zoffset[i]+j] = j+1;
    }
  }

  // loop round elements to include

  for (size_t i=0; i<(size_t)nZincl; i++) {

    // find match to input element

    int ifound = -1;
    for (size_t j=0; j<Zinput.size(); j++) {
      if ( Zinput[j] == Zinclude[i] ) ifound = j;
    }

    // update ionf - note special cases for H and He

    if ( ifound != -1 ) {
      if ( Zinclude[i] == 1 ) {
	ionf[Zoffset[i]] = IonFrac[0][1];
      } else if ( Zinclude[i] == 2 ) {
	ionf[Zoffset[i]] = IonFrac[1][2];
      } else {
	for (size_t k=0; k<IonFrac[ifound].size(); k++) {
	  ionf[Zoffset[i]+k] = IonFrac[ifound][k];
	}
      }
    }
  }

  return;

}
