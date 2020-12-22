#include <xsTypes.h>
#include <memory>
#include <cfortran.h>

PROTOCCALLSFSUB9(FMEKAL,fmekal,FLOATV,FLOATV,FLOATV,INT,FLOAT,FLOAT,FLOAT,FLOATV,PFLOAT)
#define FMEKAL(E1,E2,Flux,Nbin,Cem,H,T,Factor,Ed) \
           CCALLSFSUB9(FMEKAL,fmekal,FLOATV,FLOATV,FLOATV,INT,FLOAT,FLOAT,FLOAT,FLOATV,PFLOAT, \
           E1,E2,Flux,Nbin,Cem,H,T,Factor,Ed)

PROTOCCALLSFSUB17(FLINE,fline,FLOATV,FLOATV,FLOATV,INT,INT,INT,INT,FLOAT,FLOAT,FLOAT,FLOATV,FLOAT,FLOAT,FLOATV,FLOATV,FLOATV,PFLOAT)
#define FLINE(E1,E2,Flux,Nbin,Ic,Itype,Iem,Cem,H,T,Factor,Dilfac,Radt,Xzin,Sion,Alfa,Ed) \
            CCALLSFSUB17(FLINE,fline,FLOATV,FLOATV,FLOATV,INT,INT,INT,INT,FLOAT,FLOAT,FLOAT,FLOATV,FLOAT,FLOAT,FLOATV,FLOATV,FLOATV,PFLOAT, \
            E1,E2,Flux,Nbin,Ic,Itype,Iem,Cem,H,T,Factor,Dilfac,Radt,Xzin,Sion,Alfa,Ed)


void clcmek(int plasmaType, const RealArray& energyArray, const Real Temperature,
	    const RealArray& Abundance, const Real Density, const Real z,
            RealArray& fluxArray)
{
  // XSPEC subroutine to calculate the Mewe plasma emission spectrum.
  // Arguments
  //    plasmaType   I  i: Type - 3 = Mewe-Kaastra-Liedahl
  //                              5 = Mewe-Kaastra
  //    energyArray  R  i: energy ranges
  //    Temperature  R  i: plasma temperature (keV)
  //    Abundance    R  i: abundances (H,He,C,N,O,Ne,Na,Mg,Al,Si,S,Ar,Ca,Fe,Ni)
  //    Density      R  i: H density (cm^-3)
  //    z            R  i: redshift
  //    fluxArray    R  r: model spectrum

  // Translated from clcmek.f Feb 09 (C.Gordon)
  

  // Copy energyArray into a float* and apply redshift factor

  int ne = energyArray.size()-1;
  std::unique_ptr<float[]> apEz(new float[ne+1]);
  float* ez = apEz.get();

  // fill energy array in source frame.

  const float zfac = 1.0 + (float)z;
  for (size_t ie=0; ie<=(size_t)ne; ++ie) ez[ie] = energyArray[ie]*zfac;

  // Create a float* with abundances

  size_t NOEL = Abundance.size();
  std::unique_ptr<float[]> apAbun(new float[NOEL]);
  float* abun = apAbun.get();

  for (size_t i=0; i<NOEL; ++i) abun[i] = Abundance[i];

  // internal arrays which are required by the old FLINE routine

  const size_t NUMION = 207;
  float xz[NUMION], sion[NUMION], alfa[NUMION];

  // create float* for output flux and initialize it

  std::unique_ptr<float[]> apPhotar(new float[ne]);
  float* photar = apPhotar.get();

  for (size_t ie=0; ie<(size_t)ne; ++ie) photar[ie] = 0.0;

  const float t = (float)Temperature;
  const float h = (float)Density;

  // calculate the output emissivity

  float ed = 0.0;
  if (plasmaType == 3)
    FMEKAL(ez, ez+1, photar, ne, (float)1.0, h, t, abun, ed);
  else if (plasmaType == 5)
    FLINE(ez, ez+1, photar, ne, 1, 1, 1, (float)1.0, h, t, abun, (float)0.0, 
	  (float)0.0, xz, sion, alfa, ed);

  // convert 1e-4 ph/m**2/s/keV to ph/cm**2/s/bin for XSPEC
  // also correct for the effect that input normalisation is n_en_H V / d**2
  // (output of fline is ~ n_H n_H V / d**2)
  // 119.64 factor is to give meka model is the same normalisation as
  // the R-S model (factor = 4 * pi * (3.096)**2)
  // 1/zfac is time dilation.

  float dfact = 1.0;
  if (h*ed*zfac != 0.0) dfact = 119.64 / h / h / ed / zfac;
  for (size_t ie=0; ie<(size_t)ne; ++ie)
    fluxArray[ie] = photar[ie] * dfact*(ez[ie+1] - ez[ie]);

  return;
}
