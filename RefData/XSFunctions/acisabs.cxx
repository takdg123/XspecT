// c      FORTRAN version of xcont.pro, written by George Chartas
// c      Conversion performed by Konstantin Getman 07/17/02
// c
// c                  THE CONTENT OF xcont.pro      
// c       (1) Use observed decay with time in the ratio of observed 
// c       X-ray transmissions at 0.67keV and 5.895keV.
// c        Observed decay in S3 based on Catherine Grants' analysis:
// c        Fit parameters provided by Allyn Tennant
// c
// c        R(t) = (norm*exp(-tauinf*(1.0-exp(-t/tefold))))       eq(1)
// c          norm=0.00722+/-0.00007
// c          tefold=620.+/-66.
// c          tauinf=0.582+/-0.024
// c        
// c       
// c       (2) Model decay
// c
// c       Trans(E1,t) = TOBF(E1) * Tcont(E1,t)   eq(2)
// c
// c       where  TOBF(E1) is the transmission of the OBF at energy E1
// c        Tcont(E1,t) is the transmission of the contamination
// c        layer at energy E1 and time t
// c
// c       The modeled decay is:
// c
// c        R(t)=Trans(E1,t)/Trans(E2,t)=const*exp{-rho*d(t)*[mac(E1)-mac(E2)]} eq(3)
// c
// c
// c        where : rho is the density of the contaminant
// c        d(t) is the thickness of the contaminant
// c        mac(E) is the mass absorption coefficient of the 
// c        contaminant at energy E.
// c        mac(E) calculated from atomic scattering factor files
// c        provided at http://www-cxro.lbl.gov/optical_constants/asf.html
// c
// c      The plot of the thickness of the contaminant vs. time is 
// c       derived form eq(1) and eq(3):
// c
// c       d(t)*rho = ln(R(t)/const) / [mac(E2) - mac(E1)] 
// c
// c
// c
// c      (3) Suggested modification of ACIS arf files to account for contamination
// c
// c       Tcont(E,t) = exp(-mac(E)*rho*d(t))
// c
// c------------------------------------------------------------------------      
// c                  THE CONTENT OF acisabs xspec model
// c      Number of model parameters: 7
// c      1      Tdays      Days between launch and observation
// c       2       norm    Normalization factor
// c      3      tauinf      Optical depth at infinite time
// c      4      tefold      e-folding time for build-up of contaminant 
// c      5      nC      Number of carbon atoms in hydrocarbon molecule
// c      6      nH      Number of hydrogen atoms in hydrocarbon molecule
// c      7      nO      Number of oxygen atoms in hydrocarbon molecule
// c      8      nN      Number of nitrogen atoms in hydrocarbon molecule
// c      

// converted from Fortran to C++    kaa 1/2/18
// uses acisabsdata.fits file from modelData

#include <xsTypes.h>
#include <functionMap.h>
#include <Utilities/FunctionUtility.h>
#include <CCfits/CCfits>
#include <memory>

Real Mac(const RealArray& engrid, const RealArray& mu, const Real& enkeV);

void acisabs(const RealArray& energyArray, const RealArray& params,
             int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
             const string& initString)
{

  using namespace CCfits;

  const Real Ecal1(0.67);
  const Real Ecal2(5.895);
  const Real awC(12.01115);
  const Real awO(15.9994);
  const Real awH(1.00797);
  const Real awN(14.0067);

  static RealArray muO, muH, muC, muN, enO, enH, enC, enN;

  // Read in the data
  static bool first(true);
  if ( first ) {
    string filenm = FunctionUtility::modelDataPath() + "acisabsdata.fits";
    std::unique_ptr<FITS> pInfile((FITS*)0);
    try {
      pInfile.reset(new FITS(filenm, Read, false));
      // do the H data
      ExtHDU& HData = pInfile->extension("H");
      HData.column("ENERGY").read(enH,1,HData.rows());
      HData.column("MU").read(muH,1,HData.rows());
      muH *= 4.208e7 / enH / awH; 
      // do the C data
      ExtHDU& CData = pInfile->extension("C");
      CData.column("ENERGY").read(enC,1,CData.rows());
      CData.column("MU").read(muC,1,CData.rows());
      muC *= 4.208e7 / enC / awC; 
      // do the N data
      ExtHDU& NData = pInfile->extension("N");
      NData.column("ENERGY").read(enN,1,NData.rows());
      NData.column("MU").read(muN,1,NData.rows());
      muN *= 4.208e7 / enN / awN; 
      // do the O data
      ExtHDU& OData = pInfile->extension("O");
      OData.column("ENERGY").read(enO,1,OData.rows());
      OData.column("MU").read(muO,1,OData.rows());
      muO *= 4.208e7 / enO / awO; 
      first = false;
    } catch(...) {
      std::ostringstream msg;
      msg << "Failed to read data from " << filenm;
      FunctionUtility::xsWrite(msg.str(),5);
      return;
    }
  }

  size_t nE = energyArray.size();
  flux.resize(nE-1);
  fluxErr.resize(0);

  Real Tdays = params[0];
  Real norm = params[1];
  Real tauinf = params[2];
  Real tefold = params[3];
  Real nC = params[4];
  Real nH = params[5];
  Real nO = params[6];
  Real nN = params[7];
      
  // Calculate the relative abundance of each element

  Real xO = nO * awO;
  Real xH = nH * awH;
  Real xN = nN * awN;
  Real xC = nC * awC;
  Real elnorm = xO + xH + xN + xC;
  xO /= elnorm;
  xH /= elnorm;
  xN /= elnorm;
  xC /= elnorm;

  // Calculate mass absorption coefficient at Ecal1 and Ecal2
    
  Real MacE1 = xO*Mac(enO,muO,Ecal1)+xH*Mac(enH,muH,Ecal1)
    +xN*Mac(enN,muN,Ecal1)+xC*Mac(enC,muC,Ecal1);
  Real MacE2 = xO*Mac(enO,muO,Ecal2)+xH*Mac(enH,muH,Ecal2)
    +xN*Mac(enN,muN,Ecal2)+xC*Mac(enC,muC,Ecal2);

  // Estimate thickness of contaminate during observation 

  Real decay = norm*exp(-tauinf*(1.0 - exp(-Tdays/tefold)));
  Real thRho = log(decay/norm)/(MacE2 - MacE1);

  // Calculate transmission of contamination layer for input energies

  for (size_t i=0; i<energyArray.size()-1; i++) {
    Real eninpt = (energyArray[i] + energyArray[i+1])/2.0;
    flux[i] = exp (-(xO*Mac(enO,muO,eninpt)+xH*Mac(enH,muH,eninpt)
		     +xN*Mac(enN,muN,eninpt)+xC*Mac(enC,muC,eninpt))*thRho);
  }

  return;
}


// Function performs linear interpolation, calculates 
// mass absorption coefficient for specified energy enkeV      

Real Mac(const RealArray& engrid, const RealArray& mu, const Real& enkeV)
{
  Real eneV = 1000.0 * enkeV;

  if ( eneV <= engrid[0] ) return mu[0];
  if ( eneV >= engrid[engrid.size()-1] ) return mu[mu.size()-1];

  for (size_t i=0; i<engrid.size()-1; i++) {
    if ( eneV >= engrid[i] && eneV <= engrid[i+1] ) {
      return mu[i+1] - (engrid[i+1]-eneV)*(mu[i+1]-mu[i])/(engrid[i+1]-engrid[i]);
    }
  }

  return mu[mu.size()-1];
}
