//    Dec 98    Sedov models
//		   additive model for XSPEC to use with general models.
//		   parameters (param):
//              (0):    postshock temperature t(keV)
//              (1):    postshock electron temperature t_e(keV)
//              (2):    hydrogen abundance (switches on and off free-free
//                      continuum
//              (3..31): abundances of all elements from Z=2 to Z=30
//                        with respect to solar values
//              (32):   ionization timescale tau (cm^-3 s)
//              (33):   redshift z
//
//              Parameter (2) should be usually set to one. Its only use
//              in the code is to set hydrogen abundance to zero; in this
//              case this parameter should be set to zero. It is not the
//              density parameter.


#include <xsTypes.h>
#include <XSstreams.h>
#include <functionMap.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Utils/XSutility.h>
#include <XSUtil/Utils/XSstream.h>
#include <IonBalNei.h>
#include <Aped.h>
#include <CCfits/CCfits>
#include <memory>

using namespace std;
using namespace XSutility;
using namespace CCfits;

// prototype for routine which calculates shock temperatures and timescales

void calcSedov(const Real& emz, const Real& tshock, const Real& taus, const Real& beta,
	       RealArray& tArr, RealArray& tauArr);

// prototype for function to calculate electron/mean temperature ratio
// this is in vvnpshock.cxx

Real calcElRat(const Real& R, const Real& deltim, const Real& elRatI);


void vvsedov(const RealArray& energyArray, const RealArray& params,
	     int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
	     const string& initString)
{
  // for speed reasons set Zinput and abundance arrays only for elements
  // with non-zero abundance

  int nZ = 0;
  for (size_t i=0; i<TOTZ; i++) {
    if ( params[i+1] > 0.0 ) nZ++;
  }
  IntegerVector Zinput(nZ);
  RealArray abundance(nZ);
  int iZ = 0;
  for (size_t i=0; i<TOTZ; i++) {
    if ( params[i+1] > 0.0 ) {
      Zinput[iZ] = i+1;
      abundance[iZ] = params[i+1];
      iZ++;
    }
  }

  const size_t ntt = 5;
  const size_t nem = 400;

  static bool isFirst = true;
  static Real tkeVSave;
  static Real tekeVSave;
  static Real tauSave;
  static IntegerVector ZinputSave;

  static vector<vector<RealArray> > IonFrac;
  static vector<Real> MeanTemp;

  // if necessary set the ion fractions

  Real tkeV, tekeV, tau;
  if ( params[0] >= params[1] ) {
    tkeV = params[0];
    tekeV = params[1];
  } else {
    tkeV = params[1];
    tekeV = params[0];
  }
  tau = 0.25*params[32];
  Real beta = tekeV / tkeV;

  if ( !identicalArrays(Zinput,ZinputSave) || tkeVSave != tkeV || 
       tekeVSave != tkeV || tauSave != tau || isFirst ) {

    tkeVSave = tkeV;
    tekeVSave = tekeV;
    tauSave = tau;

    // loop over emission measures

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
    Real deleml, emb, eme;
     
    for (size_t iem=0; iem<nem-2; iem++) {
      if ( iem < nem/2-1 ) {
	deleml = 3.0/(Real)(nem/2-1);
	emb = pow(10.0,-4.0+(Real)iem*deleml);
	eme = pow(10.0,-4.0+(Real)(iem+1)*deleml);
      } else {
	size_t jem = iem - (nem/2-1);
	deleml = (4.0+log10(0.9))/(Real)(nem/2-1);
	emb = 1.0 - pow(10.0,-4.0+(Real)(jem+1)*deleml);
	eme = 1.0 - pow(10.0,-4.0+(Real)jem*deleml);
      }
      Real emz = 0.5*(emb+eme);
      Real delem = eme - emb;

      // calculate temperatures and timescales for this emission measure

      RealArray tArr, tauArr;
      try {
	calcSedov(emz, tkeV*KEVTOK, tau, beta, tArr, tauArr);
      } catch(...) {
	return;
      }
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

  Real Redshift(params[33]);
  RealArray Tinput(MeanTemp.size());
  for (size_t i=0; i<Tinput.size(); i++) Tinput[i] = MeanTemp[i];

  int status = calcNEISpectrum(energyArray, Zinput, abundance, Redshift, Tinput,
			       IonFrac, false, 0.0, flux, fluxErr);
  if ( status != 0 ) {
    std::stringstream oss;
    oss << "calcNEISpectrum failed in vsedov: status = " << status << "\n";
    oss << "set chatter to 25 and retry to get additional diagnostics.\n";
    xs_write(const_cast<char*>(oss.str().c_str()),10);
  }

  isFirst = false;
  return;

}

//------------------------------------------------------------------------------
// Calculate temperatures and ionization timescale for Sedov models
// Input is:
//   emz     fractional emission measure (instead of fractional radius)
//   tshock  mean temperature at the shock (K)
//   taus    ionization parameter at the shock (s/cm^3)
//   beta    electron/mean temperature ratio at the shock
// Output is:
//   tArr    temperatures
//   tauArr  ionization timescales

void calcSedov(const Real& emz, const Real& tshock, const Real& taus, 
	       const Real& beta, RealArray& tArr, RealArray& tauArr)
{
  static bool qfirst = true;
  static RealArray em, v, rksi, zeta, p, rho, g, t;

  // if first time through read in the Sedov solution data

  if ( qfirst ) {
    const string& datadir = FunctionUtility::modelDataPath();
    string filename = datadir + "sedov.fits";

    const string hduName("SEDOV");
    const vector<string> hduKeys;
    const vector<string> primaryKeys;
    std::unique_ptr<FITS> pSedovfile((FITS*)0);

    try {
      pSedovfile.reset(new FITS(filename,CCfits::Read,hduName,false,hduKeys,
				primaryKeys,(int)1));
    } catch(...) {
      std::stringstream oss;
      oss << "Failed to open " << filename << " for sedov model." << "\n";
      xs_write(const_cast<char*>(oss.str().c_str()),10);
      throw;
    }

    ExtHDU& sedovExt = pSedovfile->currentExtension();

    long nRows = (long)sedovExt.rows();
    sedovExt.column("EM").read(em, (long)1, nRows);
    sedovExt.column("V").read(v, (long)1, nRows);
    sedovExt.column("RKSI").read(rksi, (long)1, nRows);
    sedovExt.column("ZETA").read(zeta, (long)1, nRows);
    sedovExt.column("P").read(p, (long)1, nRows);
    sedovExt.column("RHO").read(rho, (long)1, nRows);
    sedovExt.column("G").read(g, (long)1, nRows);
    sedovExt.column("T").read(t, (long)1, nRows);

    qfirst = false;
  }

  int iz = locateIndex(em,emz);
  int iz1 = iz + 1;
  Real pp = (emz-em[iz])/(em[iz1]-em[iz]);
  Real tmin=(1.-pp)*p[iz]/rho[iz] + pp*p[iz1]/rho[iz1];
  Real zetaiz=(1.-pp)*zeta[iz] +pp*zeta[iz1];
  Real tmax=1./(zetaiz*zetaiz*zetaiz);
  Real rhop = pow(tmin/tmax,1.5);
  int j = locateIndex(rho,rhop);
  pp = (rhop-rho[j])/(rho[j+1]-rho[j]);
  Real zetap = zeta[j]*(1.-pp) + zeta[j+1]*pp;
  Real gp = g[j]*(1.-pp) + g[j+1]*pp;
  Real taumax = gp*pow(zetaiz,2.5);
  if ( taumax < 0.0 ) taumax = 0.0;
  
  //  Assumes cosmic H and He abundances (only trace of heavy elements).
  Real heabun = 1.0/pow(10.0,1.01);
  Real xmu = (1.+4.*heabun)/(2.+3*heabun);
  Real xmue = (1.+2.*heabun)/(1.+4*heabun);
  Real convn = 1./(xmu*xmue);
  Real r = 4.*convn*taus*pow(zetaiz,7.0)*(1./pow(zetaiz,2.5)-1.)/pow(tshock,1.5);
  Real delg=0.1*taumax/pow(zetaiz,2.5);

  vector<Real> tVec, tauVec;
  tauVec.push_back(0.0);
  tVec.push_back(beta*tmax);
  gp = (Real)1*delg;
  j = locateIndex(g,gp);
  pp = (gp-g[j])/(g[j+1]-g[j]);
  zetap = zeta[j]*(1.-pp) + zeta[j+1]*pp;
  rhop = rho[j]*(1.-pp) + rho[j+1]*pp;
  Real delf = (1./pow(zetap,2.5)-1.)/(1./pow(zetaiz,2.5)-1.);
  Real tl1 = calcElRat(r,delf,beta)*tmax*pow(rhop,2./3.);
  Real gp1 = gp;
  for (size_t l=2; l<=10; l++) {
    gp = (Real)l*delg;
    j = locateIndex(g,gp);
    pp = (gp-g[j])/(g[j+1]-g[j]);
    zetap = zeta[j]*(1.-pp) + zeta[j+1]*pp;
    rhop = rho[j]*(1.-pp) + rho[j+1]*pp;
    delf = (1./pow(zetap,2.5)-1.)/(1./pow(zetaiz,2.5)-1.);
    Real tl = calcElRat(r,delf,beta)*tmax*pow(rhop,2./3.);
    tVec.push_back(tl1);
    tauVec.push_back(pow(zetaiz,2.5)*gp1);
    tl1 = tl;
    gp1 = gp;
  }
  delf = 1.0;
  tVec.push_back(calcElRat(r,delf,beta)*tmin);
  tauVec.push_back(taumax);

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
