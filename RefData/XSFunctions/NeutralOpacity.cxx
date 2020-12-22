// Methods for NeutralOpacity class.

#include "NeutralOpacity.h"
#include <cfortran.h>

// Fortran functions called

PROTOCCALLSFFUN4(FLOAT,GPHOTO,gphoto,FLOAT,FLOAT,INT,INT)
#define GPHOTO(E1,E2,Z,IS) CCALLSFFUN4(GPHOTO,gphoto,FLOAT,FLOAT,INT,INT,E1,E2,Z,IS)

PROTOCCALLSFFUN5(FLOAT,PHOTO,photo,FLOAT,FLOAT,INT,INT,INT)
#define PHOTO(E1,E2,Z,IOPT,IS) CCALLSFFUN5(PHOTO,photo,FLOAT,FLOAT,INT,INT,INT,E1,E2,Z,IOPT,IS)


// default constructor

NeutralOpacity::NeutralOpacity()
{
  AtomicNumber.resize(0);
}

// default destructor

 NeutralOpacity::~NeutralOpacity()
 {
 }


 void NeutralOpacity::Setup()
 {

   // set up arrays of atomic numbers and names

   int numb[] = {1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 17, 18, 20, 24, 26, 27, 28};
   AtomicNumber.resize(18);
   for (size_t i=0; i<18; i++) AtomicNumber[i] = numb[i];
   string name[] = {"H", "He", "C", "N", "O", "Ne", "Na", "Mg",
		  "Al", "Si", "S", "Cl", "Ar", "Ca", "Cr", "Fe",
		  "Co", "Ni"};
   ElementName.resize(18);
   for (size_t i=0; i<18; i++) ElementName[i] = name[i];

   // get the cross-sections in use and place in the string

   CrossSectionSource = FunctionUtility::XSECT();

   return;

 }

void NeutralOpacity::Get(RealArray inputEnergy, Real Abundance, Real IronAbundance, bool IncludeHHe, RealArray& Opacity)
{

  size_t Nabund(AtomicNumber.size());

  RealArray abund(Nabund);
  int status(0);

  Opacity.resize(inputEnergy.size());
  Opacity = 0.0;

  // set abundances

  for (size_t i=0; i<Nabund; i++) {

    abund[i] = FunctionUtility::getAbundance(ElementName[i]);
    if ( AtomicNumber[i] == 26 ) {
      abund[i] *= pow(10.0,IronAbundance);
    } else if ( AtomicNumber[i] > 2 ) {
      abund[i] *= pow(10.0,Abundance);
    } else {
      if ( !IncludeHHe ) abund[i] = 0.0;
    }

  }

  // loop over input energies

  for (size_t ie=0; ie<inputEnergy.size(); ie++) {

    // convert to keV from m_e c^2 units

    float EkeV = inputEnergy[ie] * 511.00;

    // if above 100 keV or below Lyman limit then leave opacity as zero

    if ( EkeV > LYLIMIT && EkeV < 100.0 ) {

      // loop over elements accumulating opacity

      for (size_t j=0; j<Nabund; j++) {

	if ( CrossSectionSource == "vern" ) {
	  Opacity[ie] += abund[j] * GPHOTO(EkeV, EkeV, AtomicNumber[j], status);
	} else if ( CrossSectionSource == "bcmc" ) {
	  Opacity[ie] += abund[j] * PHOTO(EkeV, EkeV, AtomicNumber[j], 3, status);
	} else if ( CrossSectionSource == "obcm" ) {
	  Opacity[ie] += abund[j] * PHOTO(EkeV, EkeV, AtomicNumber[j], 2, status);
	}

      }

    }

  }

  // scale opacity to match units used in IonizedOpacity

  Opacity *= 1.0e22;

  return;

}

void NeutralOpacity::GetValue(Real inputEnergy, Real Abundance, Real IronAbundance, bool IncludeHHe, Real& Opacity)
{

  RealArray E(1);
  RealArray op(1);

  E[0] = inputEnergy;

  Get(E, Abundance, IronAbundance, IncludeHHe, op);

  Opacity = op[0];

  return;

}
