//     Aug 99    Plane-parallel shock models with unequal ion and electron 
//               temperatures. May be used to model a shock section.
//		    additive model for XSPEC to use with general models.
//		    parameters (param):
//               (0):    postshock temperature t(keV)
//               (1):    postshock electron temperature t_e(keV)
//               (2):    hydrogen abundance (switches on and off free-free
//                       continuum
//               (3..31): abundances of all elements from Z=2 to Z=30
//                                  with respect to solar values
//               (32):   ionization parameter taul(cm^-3 s)
//               (33):   ionization parameter tauu(cm^-3 s)
//               (34):   redshift z
//
//               Parameter (2) should be usually set to one. Its only use
//               in the code is to set hydrogen abundance to zero; in this
//               case this parameter should be set to zero. It is not the
//               density parameter.

#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/XSstream.h>
#include <IonBalNei.h>
#include <Aped.h>

using namespace std;
using namespace XSutility;

// prototype for routine which calculates shock temperatures and timescales

void calcShock(const Real& tshock, const Real& taus, const Real& beta,
	       RealArray& tArr, RealArray& tauArr);

// prototype for function to calculate electron/mean temperature ratio

Real calcElRat(const Real& R, const Real& deltim, const Real& elRatI);



void vvnpshock(const RealArray& energyArray, const RealArray& params,
	       int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	       const string& initString)
{
  // for speed reasons set Zinput and abundance arrays only for elements
  // with non-zero abundance

  int nZ = 0;
  for (size_t i=0; i<TOTZ; i++) {
    if ( params[i+2] > 0.0 ) nZ++;
  }
  IntegerVector Zinput(nZ);
  RealArray abundance(nZ);
  int iZ = 0;
  for (size_t i=0; i<TOTZ; i++) {
    if ( params[i+2] > 0.0 ) {
      Zinput[iZ] = i+1;
      abundance[iZ] = params[i+2];
      iZ++;
    }
  }

  const size_t ntt = 5;

  static bool isFirst = true;
  static Real tkeVSave;
  static Real tekeVSave;
  static Real tauLowSave;
  static Real tauHighSave;
  static IntegerVector ZinputSave;

  static vector<vector<RealArray> > IonFrac;
  static vector<Real> MeanTemp;

  // if necessary set the ion fractions

  Real tkeV, tekeV, tauLow, tauHigh;
  if ( params[0] >= params[1] ) {
    tkeV = params[0];
    tekeV = params[1];
  } else {
    tkeV = params[1];
    tekeV = params[0];
  }
  if ( params[33] >= params[32] ) {
    tauLow = params[32];
    tauHigh = params[33];
  } else {
    tauLow = params[33];
    tauHigh = params[32];
  }
  Real beta = tekeV / tkeV;

  // if necessary recalculate the ion fractions

  if ( !identicalArrays(Zinput,ZinputSave) || tkeVSave != tkeV || 
       tekeVSave != tkeV || tauLowSave != tauLow || tauHighSave != tauHigh || 
       isFirst ) {

    tkeVSave = tkeV;
    tekeVSave = tekeV;
    tauLowSave = tauLow;
    tauHighSave = tauHigh;

    // set up zones

    size_t nzones;
    Real taumin, deltau, tauiml;

    if ( tauLow == 0.0 ) {

      nzones = 200;
      taumin = 1.0e8;
      deltau = log10(tauHigh/taumin)/(Real)(nzones-1);
      tauiml = 0.0;

    } else {

      nzones = 20;
      deltau = (tauHigh - tauLow)/(Real)nzones;
      tauiml = tauLow;

    }

    // loop over zones

    size_t nT = getNEInumbTemp();
    vector<vector<RealArray> > IonFracAcc(nT);
    for (size_t it=0; it<nT; it++) {
      IonFracAcc[it].resize(nZ);
      for (size_t i=0; i<(size_t)nZ; i++) {
	IonFracAcc[it][i].resize(Zinput[i]+1);
	IonFracAcc[it][i] = 0.0;
      }
    }
    RealArray delarr(0.0, nT);
    vector<RealArray> tIonFrac;

    Real taui, delem;
    for (size_t izone=0; izone<nzones; izone++) {

      if ( tauLow == 0.0 ) {
	taui = taumin * pow(10.0, izone*deltau);
	delem = (taui - tauiml) / tauHigh;
      } else {
	taui = tauLow + (izone+1)*deltau;
	delem = (taui - tauiml) / (tauHigh - tauLow);
      }
      Real tausiz = 0.5*(taui+tauiml);
      tauiml = taui;

      // calculate temperatures and timescales for the shock

      RealArray tArr, tauArr;
      calcShock(tkeV*KEVTOK, tausiz, beta, tArr, tauArr);
      tArr /= KEVTOK;
      calcNEIfractions(tArr, tauArr, Zinput, tIonFrac);

      // load into the appropriate temperature bin

      int it = getNEItempIndex(tArr[tArr.size()-1]);
      for (size_t i=0; i<(size_t)nZ; i++) {
	for (size_t j=0; j<(size_t)Zinput[i]+1; j++) {
	  IonFracAcc[it][i][j] += tIonFrac[i][j]*delem;
	}
      }
      delarr[it] += delem;

      // end loop over zones
    }

    // Rebin temperature zones for spectral calculations. This cuts
    // down the numbers of separate temperatures and ion fractions
    // by up to a factor of 5.

    MeanTemp.resize(0);
    IonFrac.resize(0);

    for (size_t i=0; i<nT-ntt; i+=ntt) {
      Real eitot(0.0);
      for (size_t j=i; j<i+ntt; j++) {
	if ( delarr[j] > 0.0 ) eitot = 1.0;
      }
      if ( eitot > 0.0 ) {
	Real sumem(0.0);
	Real sumt(0.0);
	tIonFrac.resize(nZ);
	for (size_t j=0; j<(size_t)nZ; j++) tIonFrac[j] = 0.0;
	for (size_t it=i; it<i+ntt; it++) {
	  Real tt = pow(10.0,MINLOGT+it*DELTALOGT)/KEVTOK;
	  sumem += delarr[it];
	  sumt += delarr[it]*tt;
	  for (size_t j=0; j<(size_t)nZ; j++) tIonFrac[j] += IonFracAcc[it][j];
	}
	IonFrac.push_back(tIonFrac);
	MeanTemp.push_back(sumt/sumem);
      }
    }

    ZinputSave.resize(Zinput.size());
    ZinputSave = Zinput;

    // end set up of ionization fractions
  }

  // calculate output spectrum

  Real Redshift(params[34]);
  RealArray Tinput(MeanTemp.size());
  for (size_t i=0; i<Tinput.size(); i++) Tinput[i] = MeanTemp[i];

  int status = calcNEISpectrum(energyArray, Zinput, abundance, 
			       Redshift, Tinput, IonFrac, false, 0.0, flux, fluxErr);
  if ( status != 0 ) {
    std::ostringstream oss;
    oss << "calcNEISpectrum failed in vvnpshock: status = " << status << "\n";
    oss << "set chatter to 25 and retry to get additional diagnostics.\n";
    xs_write(const_cast<char*>(oss.str().c_str()),10);
  }

  isFirst = false;

  return;

}

//------------------------------------------------------------------------------
// Calculate temperatures and ionization timescales from the shock. 
// Input is:
//   tshock  mean temperature at the shock
//   taus    ionization parameter at the shock
// Output is:
//   tArr    temperatures
//   tauArr  ionization timescales

void calcShock(const Real& tshock, const Real& taus, const Real& beta,
	       RealArray& tArr, RealArray& tauArr)
{

  Real tmin = 1.0;
  Real tmax = 1.0;
  Real taumax = 1.0;

  vector<Real> tVec, tauVec;

  if ( beta >= 1.0 ) {
    tVec.push_back(tmin);
    tVec.push_back(tmax);
    tauVec.push_back(0.0);
    tauVec.push_back(taumax);
  } else if ( beta <  1.0 ) {
    Real heabun = 1.0/pow(10.0,1.01);
    Real xmu = (1.+4.*heabun)/(2.+3*heabun);
    Real xmue = (1.+2.*heabun)/(1.+4*heabun);
    Real convn = 1./(xmu*xmue);
    Real r = convn*taus/pow(tshock,1.5);
    tauVec.push_back(0.0);
    tVec.push_back(beta*tmax);
    Real delf;
    Real delf1 = 0.01;
    Real tl1 = calcElRat(r,delf1,beta)*tmax;
    for (size_t l=2; l<=100; l++) {
      delf = (Real)(l)*0.01;
      Real tl = calcElRat(r,delf,beta)*tmax;
      if ( abs((tl-tVec[tVec.size()-1])/tl) > 0.02 ) {
	tauVec.push_back(delf1);
	tVec.push_back(tl1);
      }
      delf1 = delf;
      tl1 = tl;
    }
    delf = 1.0;
    tVec.push_back(calcElRat(r,delf,beta)*tmin);
    tauVec.push_back(taumax);
  }

  // Normalize temperature and ionization timescale to physical units.

  size_t npt = tVec.size();
  tArr.resize(npt);
  tauArr.resize(npt);
  for (size_t i=0; i<npt; i++) {
    tArr[i] = tshock*tVec[i];
    tauArr[i] = taus*tauVec[i];
  }

  return;

}

//------------------------------------------------------------------------------
//  Calculates electron/mean temperature ratio for adiabatic process with
//  Coulomb collisions between ions and electrons. 
//      Relative accuracy of the fit - 7x10**-5 (max. dev.).
//    INPUT: R       - total particle density/mean temperature^3/2 (cgs)
//           deltim  - elapsed time (s)
//           elRatI  - initial electron/mean temperature ratio

#define TOT 2.0/3.0
#define EOT 8.0/3.0
#define OOS 1.0/7.0
#define OON 1.0/9.0
#define COULOMB 30.0
#define CONSTANT 2.0 * COULOMB / 503.0

Real calcElRat(const Real& R, const Real& deltim, const Real& elRatI)
{

  if ( elRatI == 1.0 ) return elRatI;

  Real x0 = elRatI;
  Real sx0 = sqrt(x0);
  Real y0;
  if ( x0 > 0.01 ) {
    y0 = log( (1.+sx0)/(1.-sx0) ) - TOT*sx0*(x0+3.);
  } else {
    y0 = 2.*sx0*x0*x0*( 0.2 + x0*( OOS + x0*OON ) );
  }

  Real y = CONSTANT*deltim*R + y0;
  if ( y < 3.95 ) {
    Real z = pow(2.5*y,0.4);
    return z*(1.0000271+z*(-0.28649804+z*(0.58116559E-02+z*(
	   -0.22288787E-02+z*(0.48568998E-02-z*0.78976963E-03)))));
  } else {
    return pow(tanh(0.5*(y+EOT)),2.0);
  }
}

