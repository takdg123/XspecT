// Methods for IonizedOpacity class.

#include "IonizedOpacity.h"
#include <XSFunctions/Utilities/xsFortran.h>
#include <CCfits/CCfits>
#include <memory>

// default constructor

IonizedOpacity::IonizedOpacity()
{
  AtomicNumber.resize(0);
}

// default destructor

 IonizedOpacity::~IonizedOpacity()
 {
 }


void IonizedOpacity::LoadFiles()
{

  using namespace CCfits;

  // get memory for arrays - there are 10 different elements at present

  ion.resize(10);
  sigma.resize(10);
  AtomicNumber.resize(10);
  ElementName.resize(10);
  Energy.resize(721);

  //  read in data from mansig.fits

  string absname = FunctionUtility::modelDataPath() + "mansig.fits";
  std::unique_ptr<FITS> pInfile((FITS*)0);
  try {
    pInfile.reset(new FITS(absname, Read, false));
    ExtHDU& sData = pInfile->extension("SIGMAS");
    size_t nRows(sData.rows());
    vector<int> ZNumber(nRows), IonNumber(nRows);
    vector<RealArray> IonData(nRows), SigmaData(nRows);
    for (size_t i=0; i<nRows; i++) {
      IonData[i].resize(10);
      SigmaData[i].resize(721);
    }
    sData.column("Z").read(ZNumber, 1, nRows);
    sData.column("ION").read(IonNumber, 1, nRows);
    sData.column("SIGMA").readArrays(SigmaData, 1, nRows);
    sData.column("IONDATA").readArrays(IonData, 1, nRows);
    ExtHDU& eData = pInfile->extension("ENERGIES");
    eData.column("ENERGY").read(Energy, 1, eData.rows());
    
    // remap onto internal arrays

    size_t iZ(-1);
    int currentZ(-1);
    size_t iIon(-1);
    for (size_t i=0; i<ZNumber.size(); i++) {
      if ( ZNumber[i] != currentZ ) {
	iZ++;
	AtomicNumber[iZ] = ZNumber[i];
	currentZ = ZNumber[i];
	ion[iZ].resize(AtomicNumber[iZ]);
	sigma[iZ].resize(AtomicNumber[iZ]);
	iIon = -1;
      }
      iIon++;
      ion[iZ][iIon].resize(10);
      for (size_t k=0; k<10; k++) {
	ion[iZ][iIon][k] = IonData[i][k];
      }
      // change units of coefficients
      ion[iZ][iIon][1] *= 1.0E+10;
      ion[iZ][iIon][3] *= 1.0E+04;
      ion[iZ][iIon][4] *= 1.0E-04;
      ion[iZ][iIon][6] *= 1.0E-04;

      sigma[iZ][iIon].resize(721);
      for (size_t k=0; k<721; k++) {
	sigma[iZ][iIon][k] = SigmaData[i][k]/6.6e-27;
      }

    }
  } catch(...) {
    std::ostringstream msg;
    msg << "Failed to read data from " << absname;
    FunctionUtility::xsWrite(msg.str(), 5);
    return;
  }

  // set names of elements

  ElementName[0] = "H";
  ElementName[1] = "He";
  ElementName[2] = "C";
  ElementName[3] = "N";
  ElementName[4] = "O";
  ElementName[5] = "Ne";
  ElementName[6] = "Mg";
  ElementName[7] = "Si";
  ElementName[8] = "S";
  ElementName[9] = "Fe";

  return;
}

 void IonizedOpacity::Setup(Real Xi, Real Temp, RealArray inputEnergy, RealArray inputSpectrum)
 {

   // if the files have not been read do so

   if ( AtomicNumber.size() == 0 ) {

     this->LoadFiles();

     // set memory for population numbers array

     num.resize(AtomicNumber.size());

   }

   // calculate the population numbers for Xi, Temp and inputSpectrum.

   //  set up input spectrum on the cross-section energy grid
   //  inputSpectrum is E F(E) in m_e c^2 units. spec array is F(E) Delta E.


   size_t Nenergy = Energy.size();
   RealArray spec(0.0,Nenergy);
   RealArray evEnergy(inputEnergy.size());

   evEnergy = inputEnergy * 5.11e5;

   size_t i = 0;
   Real denom = 0.;
   for (size_t k=0; k<Nenergy; k++) {
     if ( Energy[k] >= evEnergy[0] && Energy[k] <= evEnergy[evEnergy.size()-1] ) {
       while ( evEnergy[i] < Energy[k] ) i++;
       spec[k] = ( inputSpectrum[i-1]*(evEnergy[i]-Energy[k]) +
		   inputSpectrum[i]*(Energy[k]-evEnergy[i-1]) ) 
	         / (evEnergy[i]-evEnergy[i-1]) / Energy[k];
       if ( k == 0 ) {
	 spec[k] *= (Energy[k+1]-Energy[k]);
       } else if ( k == Nenergy-1 ) {
	 spec[k] *= (Energy[k]-Energy[k-1]);
       } else {
	 spec[k] *= (Energy[k+1]-Energy[k-1])/2.0;
       }
       denom += spec[k];
     }
   }

   // normalise photon spectrum

   for (size_t k=0; k<Nenergy; k++) spec[k] /= denom*Energy[k];

   // some handy variables

   Real xil;
   if ( Xi <= 0.0 ) {
     xil = -100.0;
   } else {
     xil = log(Xi);
   }
   Real t4 = Temp*.0001;
   Real tfact = 1.033E-3/sqrt(t4);

   // loop over elements

   for (size_t i=0; i<AtomicNumber.size(); i++) {

     // calculate the photoionization integrals and population balance
     Real mult = 0.0;
     Real ratsum = 0.0;

     size_t Ne = (size_t)AtomicNumber[i];

     RealArray mul(Ne);
     RealArray ratio(Ne);

     for (size_t j=0; j<Ne; j++) {

       //  do integral in log10 (eV) from lowest energy
       //  to the maximum energy- use same energy grid
       //  as the cross-sections are defined on

       Real intgral = 0.0;
       for (size_t k=0; k<Nenergy; k++) {
	 intgral += sigma[i][j][k]*spec[k];
       }

       // calculate recombination coefficients - radiative and dielectric
       // from Aldrovandi+Pequignot 1973 Astro, astrophys 25 137

       Real arec;
       if ( j < Ne-1 ) {
	 Real e1 = exp(-ion[i][j][4]/t4);
	 Real e2 = exp(-ion[i][j][6]/t4);
	 arec = ion[i][j][1] * pow(t4,-ion[i][j][2])
	   + ion[i][j][3] * pow(t4,-1.5) * e1 * (1.0+ion[i][j][5]*e2);

       } else {

	 //  do radiative recomb to hydrogenic ion NB assumed gaunt factor=1
	 //  from Gould+Thakur 1970 ann phys 61 351.

	 Real z2 = AtomicNumber[i]*AtomicNumber[i];
	 Real y = 15.8*z2/t4;
	 arec = tfact*z2*(1.735+log(y)+1/(6.*y));

       }

       //  const=7.958d-2/2.43d+4

       ratio[j] = log(3.2749e-6*intgral/arec);
       ratsum += ratio[j];

       //  calculate population numbers
       //  set sum=1 as the sum=1+A21+A32A21+....
       //  complicated because floating range limit

       mul[j] = ratsum + (j+1)*xil;
       if ( mul[j] > mult ) mult = mul[j];

     }

     Real sum = 0.0;
     for (size_t j=0; j<Ne; j++) {
       mul[j] -=  mult;
       sum += exp(mul[j]);
     }

     sum += exp(-mult);

     num[i].resize(Ne+1);

     num[i][0] = -mult - log(sum);
     for (size_t j=1; j<=Ne; j++) {
       num[i][j] = num[i][j-1] + ratio[j-1] + xil;
     }
     for (size_t j=0; j<=Ne; j++) num[i][j] = exp(num[i][j]);

   }

   return;

 }

void IonizedOpacity::Get(RealArray inputEnergy, Real Abundance, Real IronAbundance, bool IncludeHHe, RealArray& Opacity)
{

  RealArray abund(AtomicNumber.size());


  Opacity.resize(inputEnergy.size());
  Opacity = 0.0;

  // set abundances

  for (size_t i=0; i<abund.size(); i++) {

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

  size_t k = 0;
  Real EnergyMax = Energy[Energy.size()-1];
  Real EnergyMin = Energy[0];

  for (size_t ie=0; ie<inputEnergy.size(); ie++) {

    Real xp = inputEnergy[ie] * 5.11e5;

    if ( xp >= EnergyMax ) {

      Real re = pow((xp/EnergyMax),-3.0);
      for (size_t i=0; i<AtomicNumber.size(); i++) {
	Real rp = abund[i]*re;
	for (size_t j=0; j<(size_t)AtomicNumber[i]; j++) {
	  Opacity[ie] += num[i][j]*sigma[i][j][720]*rp;
	}
      }

    } else if ( xp < EnergyMin ) {

      Opacity[ie] = 1.0e20;

    } else {

      // pick closest energy bin

      while ( Energy[k] < xp ) k++;

      Real re = (xp-Energy[k-1])/(Energy[k]-Energy[k-1]);
      for (size_t i=0; i<AtomicNumber.size(); i++) {
	Real rp1 = abund[i] * re;
	Real rp2 = abund[i] * (1.0 - re);
	for (size_t j=0; j<(size_t)AtomicNumber[i]; j++) {
	  Opacity[ie] += num[i][j]*(rp1*sigma[i][j][k] + rp2*sigma[i][j][k-1]);
	}
      }

    }

  }

  Opacity *= 6.6e-5;

  return;

}

void IonizedOpacity::GetValue(Real inputEnergy, Real Abundance, Real IronAbundance, bool IncludeHHe, Real& Opacity)
{

  RealArray E(1);
  RealArray op(1);

  E[0] = inputEnergy;

  Get(E, Abundance, IronAbundance, IncludeHHe, op);

  Opacity = op[0];

  return;

}
