// Methods for the Magzdiarz and Zdziarski Compton reflection models class

#include "MZCompRefl.h"
#include <XSUtil/Numerics/AdaptiveIntegrate.h>

// default constructor

MZCompRefl::MZCompRefl()
{
  X.resize(0);
  Spinc.resize(0);

  pm1.resize(19);
  pmx.resize(25);
  ap.resize(3);

  ap[2] = 0.381;
}

// default destructor

MZCompRefl::~MZCompRefl()
{
}

// Calculate the reflection spectrum

void MZCompRefl::CalcReflection(string RootName, Real cosIncl, Real xnor, Real inXmax, 
				RealArray& InputX, RealArray& InputSpec, RealArray& Spref)
{
  using namespace Numerics;

  pm1y(cosIncl,pm1);
  pmxy(cosIncl,pmx);

  ap[0] = xnor/1.21;
  ap[1] = apf(cosIncl);

  pmymax = 1-sqrt(1-(cosIncl-.05)*(cosIncl-.05));

  // set the integration defaults

  Precision = 0.01;
  Xmax = inXmax;

  string pname = RootName + "_PRECISION";
  string pvalue(FunctionUtility::getModelString(pname));

  if ( pvalue.length() && pvalue != FunctionUtility::NOT_A_KEY() ) {
    std::istringstream iss(pvalue);
    if ( !(iss >> Precision) || !iss.eof() ) {
      std::ostringstream err;
      err << "Invalid " << RootName << "_PRECISION value: " << pvalue;
      xs_write(const_cast<char*>(err.str().c_str()), 10);
    }
  }

  // load the input spectrum

  X.resize(InputX.size());
  Spinc.resize(InputX.size());
  X = InputX;
  Spinc = InputSpec;

  // loop over energies

  int totalEval = 0;
  int totalEnergiesCalc = 0;

  Spref = 0.0;

  for (size_t j=0; j<X.size(); j++) {

      Real Xval(X[j]);

      // for energies between 0.01957 and Xmax

      if ( Xval >= XTRANSL && Xval <= Xmax ) {

	Real y = 1/Xval;

	IntegrandSetup(y);

	// evaluate the integral from y0 = y-max(dym,2) to y-pmymax

	Real dym = y - YMIN;
	Real dy = dym;
	if ( dy > 2 ) dy = 2;

	Real sr = 0;
	Real error = 0;

	// adaptive integration to given precision using Gauss-Kronrod

	//	int order = AdaptiveIntegral(y-dy, y-pmymax, sr, error);
	RealArray Parameters(0);
	int order = AdaptiveIntegrate<MZCompReflIntegrand>(y-dy, y-pmymax, this, 
							   Precision, sr, error);
	totalEval += (2*order+1) * 15;

	// evaluate the integral from y0 = y-dym to y-2

	Real srp = 0.;
	if ( dym > 2 ) {

	  // adaptive integration to given precision using Gauss-Kronrod

	  int order = AdaptiveIntegrate<MZCompReflIntegrand>(y-dym, y-2.0, this, 
							     Precision, srp, error);
	  totalEval += (2*order+1) * 15;

	  sr += srp;

	}

	Spref[j] = sr*Xval;

	totalEnergiesCalc++;

      }
  }

  return;
}

// private routines

// get the spectrum for input values of y = 1/X

void MZCompRefl::InterpolatedSpectrum(const RealArray& y0, RealArray& sout)
{

  sout = 0.0;

  // Performs binary search to find first X matching y0 then steps
  // through doing interpolations.
  // presumes monotonic increasing array X

  // Find the first target y0 which lies within the limits of the X array

  size_t first = 0;
  Real Xtarg = 1.0/y0[first];
  while ( Xtarg < X[0] && Xtarg > X[X.size()-1] ) {
    first++;
    Xtarg = 1.0/y0[first];
  }

  // Find the first interpolation point by binary search

  size_t ixlo, ixhi;

  ixlo = 0;
  ixhi = X.size()-1;
  while ( ixhi-ixlo > 1 ) {
    size_t ix = (ixhi + ixlo) / 2;
    if ( Xtarg > X[ix] ) {
      ixlo = ix;
    } else {
      ixhi = ix;
    }
  }

  // Now step through the interpolations

  for (size_t i=first; i<y0.size(); i++) {

    if ( i != first ) {
      Xtarg = 1.0/y0[i];

      if ( Xtarg >= 1.0/y0[i-1] ) {
	while ( X[ixhi] < Xtarg ) {
	  ixlo = ixhi;
	  ixhi++;
	}
      } else {
	while ( X[ixlo] > Xtarg ) {
	  ixhi = ixlo;
	  ixlo--;
	}
      }
    }

    if ( Xtarg < X[0] || Xtarg > X[X.size()-1] ) {
      sout[i] = 0.0;
    } else if ( ixlo != ixhi ) {
      sout[i] = ( Spinc[ixlo]*(X[ixhi]-Xtarg) + 
		  Spinc[ixhi]*(Xtarg-X[ixlo]) ) / (X[ixhi]-X[ixlo]);
    } else {
      sout[i] = Spinc[ixlo];
    }

  }

  return;
}

void MZCompRefl::InterpolatedValue(const Real y, Real& sout)
{
  RealArray yarray(1);
  RealArray sarray(1);

  yarray[0] = y;

  InterpolatedSpectrum(yarray, sarray);

  sout = sarray[0];

  return;
}

#include <iostream>
#include <fstream>

void MZCompRefl::WriteSpectrum(const string filename)
{

  ofstream debugfile;
  debugfile.open (filename.c_str());
  for (size_t i=0; i<X.size(); i++) debugfile << X[i] << " " << Spinc[i] << endl;
  debugfile.close();

  return;
}
void MZCompRefl::IntegrandSetup(const Real yIn)
{
  ysource = yIn;
}

RealArray MZCompRefl::CalculateIntegrand(const RealArray& y0)
{
  RealArray OutArray(y0.size());
  RealArray sout(y0.size());
  RealArray grout(y0.size());

  InterpolatedSpectrum(y0,sout);
  grxy(pm1,pmx,ysource,y0,ap,grout);

  OutArray.resize(y0.size());
  OutArray = sout * grout;

  return OutArray;
}

// other functions required by these classes

#include "IonizedOpacity.h"
#include "NeutralOpacity.h"

//****************************************************************************
// Routine to calculate the total flux. This is in common between ireflct, pexriv and bexriv

void calcCompReflTotalFlux(string ModelName, Real Scale, Real cosIncl, Real Abund, 
			   Real FeAbund, Real Xi, Real Temp, Real inXmax, 
			   RealArray& X, RealArray& Spinc, RealArray& Sptot)
{

  static IonizedOpacity warmabs;
  static NeutralOpacity coldabs;

  Real absScale = abs(Scale);

  size_t nX(X.size());

  Sptot = 0.0;

  // Add the direct component if required

  if ( Scale >= 0.0 ) Sptot += Spinc;


  // Calculate and add the reflected component if required

  if ( Scale != 0.0 ) {

    // set up opacities

    if ( Xi > 0.0 ) {
      warmabs.Setup(Xi, Temp, X, Spinc);
    } else {
      coldabs.Setup();
    }

    // and calculate them

    Real xas(0.01947);
    Real xnor(0.0); 
    bool qHHe = false;

    if ( Xi > 0.0 ) {
      warmabs.GetValue(xas, Abund, FeAbund, qHHe, xnor);
    } else {
      coldabs.GetValue(xas, Abund, FeAbund, qHHe, xnor);
    }
    xnor *= 1.503e2*xas*xas*xas;

    RealArray opac(0.0,nX); 
    if ( Xi > 0.0 ) {
      warmabs.Get(X, Abund, FeAbund, qHHe, opac);
    } else {
      coldabs.Get(X, Abund, FeAbund, qHHe, opac);
    }
    opac *= 1.503e2;

    // calculate the non-relativistic Compton reflection factors

    RealArray NonRel(0.0,nX);
    CalcNonRelComp(cosIncl, X, opac, NonRel);

    // set the maximum X value required in the reflected component

    Real Xmax(inXmax);
    Real xrefmax = 1/(YMIN+1-sqrt(1-(cosIncl-.05)*(cosIncl-.05)));
    if ( Xmax < 0.0 ) {
      Xmax = xrefmax;
    } else {
      if ( xrefmax < Xmax ) Xmax = xrefmax;
    }

    // do the Greens' function calculation for the relativistic Compton reflection

    MZCompRefl greenir;
    RealArray Spref(nX);
    greenir.CalcReflection(ModelName, cosIncl, xnor, Xmax, X, Spinc, Spref);

    // add the reflected component into the output spectrum

    for (size_t j=0; j<nX; j++) {

      Real Xval(X[j]);

      if ( Xval < XMIN ) {

      } else if ( Xval >= XMIN && Xval < XTRANSL ) {

	Sptot[j] += absScale*Spinc[j]*NonRel[j];

      } else if ( Xval < XTRANSH ) {

	Real fjc = sin(M_PI*((log(Xval)-XJL)/XJD-0.5));
	Sptot[j] += absScale*0.5*((1-fjc)*Spinc[j]*NonRel[j] + (1+fjc)*Spref[j]);

      } else if ( Xval <= Xmax ) {

	Sptot[j] += absScale*Spref[j];

      }

    }

  }

  return;

}

//****************************************************************************
// Routine to calculate the total flux. This is in common between ireflct, pexriv and bexriv
// This is a version included for comparison with older code which used a power-law
// approximation to calculate ionized opacities if Xi > 0.


void calcCompReflTotalFluxOld(string ModelName, Real Scale, Real cosIncl, Real Abund, 
			      Real FeAbund, Real Xi, Real Temp, Real inXmax, Real Gamma,
			      RealArray& X, RealArray& Spinc, RealArray& Sptot)
{

  static IonizedOpacity warmabs;
  static NeutralOpacity coldabs;

  Real absScale = abs(Scale);

  size_t nX(X.size());

  Sptot = 0.0;

  // Add the direct component if required

  if ( Scale >= 0.0 ) Sptot += Spinc;


  // Calculate and add the reflected component if required

  if ( Scale != 0.0 ) {

    // set up opacities using a spectrum based on the input power-law index Gamma
    RealArray Spgamma(Spinc.size());
    for (size_t i=0; i<X.size(); i++) Spgamma[i] = pow(X[i],-Gamma);

    if ( Xi > 0.0 ) {
      warmabs.Setup(Xi, Temp, X, Spgamma);
    } else {
      coldabs.Setup();
    }

    // and calculate them

    Real xas(0.01947);
    Real xnor(0.0); 
    bool qHHe = false;

    if ( Xi > 0.0 ) {
      warmabs.GetValue(xas, Abund, FeAbund, qHHe, xnor);
    } else {
      coldabs.GetValue(xas, Abund, FeAbund, qHHe, xnor);
    }
    xnor *= 1.503e2*xas*xas*xas;

    RealArray opac(0.0,nX); 
    if ( Xi > 0.0 ) {
      warmabs.Get(X, Abund, FeAbund, qHHe, opac);
    } else {
      coldabs.Get(X, Abund, FeAbund, qHHe, opac);
    }
    opac *= 1.503e2;

    // calculate the non-relativistic Compton reflection factors

    RealArray NonRel(0.0,nX);
    CalcNonRelComp(cosIncl, X, opac, NonRel);

    // set the maximum X value required in the reflected component

    Real Xmax(inXmax);
    Real xrefmax = 1/(YMIN+1-sqrt(1-(cosIncl-.05)*(cosIncl-.05)));
    if ( Xmax < 0.0 ) {
      Xmax = xrefmax;
    } else {
      if ( xrefmax < Xmax ) Xmax = xrefmax;
    }

    // do the Greens' function calculation for the relativistic Compton reflection

    MZCompRefl greenir;
    RealArray Spref(nX);
    greenir.CalcReflection(ModelName, cosIncl, xnor, Xmax, X, Spinc, Spref);

    // add the reflected component into the output spectrum

    for (size_t j=0; j<nX; j++) {

      Real Xval(X[j]);

      if ( Xval < XMIN ) {

      } else if ( Xval >= XMIN && Xval < XTRANSL ) {

	Sptot[j] += absScale*Spinc[j]*NonRel[j];

      } else if ( Xval < XTRANSH ) {

	Real fjc = sin(M_PI*((log(Xval)-XJL)/XJD-0.5));
	Sptot[j] += absScale*0.5*((1-fjc)*Spinc[j]*NonRel[j] + (1+fjc)*Spref[j]);

      } else if ( Xval <= Xmax ) {

	Sptot[j] += absScale*Spref[j];

      }

    }

  }

  return;

}

//****************************************************************************
void CalcNonRelComp(const Real cosIncl, const RealArray& X, const RealArray& opac, RealArray& NonRel)
{
  // calculate the reflection factors for nonrelativistic Compton reflection
  //  reflection orders >1 approximated with isotropic phase function,
  //  and approximation of H-function by Basko ApJ 223,268,(1978),
  //  valid up to 14keV;

  Real xl = log(1+1/cosIncl);
  Real x2 = cosIncl*cosIncl;
  Real x4 = x2*x2;

  for (size_t i=0; i<X.size(); i++) {
    Real al0 = 1/(1+opac[i]/1.21);
    Real al1 = sqrt(1-al0);
    Real al2 = SQ3*al1;
    NonRel[i] = .5*al0*cosIncl*(.375*((3-2*x2+3*x4)*xl+(3*x2-1)*(.5-cosIncl))
	       +((1+(1/al1-1)*(1-log(1+al2)/al2))*(1+cosIncl*SQ3)/(1+cosIncl*al2)-1)*xl);
  }

  return;

}

//****************************************************************************
void grxy(const RealArray& pm1, const RealArray& pmx, const Real y, const RealArray& y0array, const RealArray& ap, RealArray& grout)
{

//***********************************************************************
//  Green's functions for Compton reflection with asymtotic
//  absorption included; output is y*G(mu,y,y0)*W(mu,y,y0),
//  for .0316<y0<100.; in two parts according to y - y0 is < or > 2.
//  pm1   - vector of angle-dependent coefficients,
//          set-up by pm1y subroutine
//  pmx   - vector of angle-dependent coefficients,
//          set-up by pmxy subroutine
//  y     - final wavelength y=1/x
//  y0    - incident photon wavelength y0=1/x0
//  ap    - absorption parameters

  RealArray pyx(10);
  Real gr;

  for (size_t i=0; i<y0array.size(); i++) {

    Real y0 = y0array[i];
    Real dy = y - y0;

    if ( dy <= 2.0 ) {

      Real p1 = y0;
      if ( p1 > 10.0 ) p1 = 10.0;
      p1 = log(p1);
      Real p2 = p1*p1;
      Real p3 = p2*p1;
      pyx[0] = p1*(pmx[0]+pmx[1]*p1+pmx[2]*p2+pmx[3]*p3);
      pyx[1] = p1*(pmx[4]+pmx[5]*p1+pmx[6]*p2+pmx[7]*p3);
      pyx[2] = p1*(pmx[8]+pmx[9]*p1+pmx[10]*p2+pmx[11]*p3);
      pyx[3] = p1*(pmx[12]+pmx[13]*p1+pmx[14]*p2+pmx[15]*p3+pmx[16]*p3*p1);
      pyx[8] = 1 + pmx[17]*(.326*(pow(y0*.1,1.692)-1)+p1-2.303);
      if ( y0 < 10. ) {
	p1 = y0 + 2;
	p1 = log(p1/y)/log(p1/y0);
	p2 = p1*p1;
	p3 = p2*p1;
	gr = exp(pyx[0]+pyx[1]*p1+pyx[2]*p2+pyx[3]*p3);
      } else {
	Real p4 = 1/(10+dy);
	p1 = log(12*p4)*5.4848;
	p2 = p1*p1;
	p3 = p2*p1;
	p4 = pyx[8]*(y0+dy)*p4;
	gr = p4*exp(pyx[0]+pyx[1]*p1+pyx[2]*p2+pyx[3]*p3);
      }
      p1 = pm1[0] + pm1[1]*dy;
      p3 = pm1[4]-dy;
      if ( p3 < 0 ) p3 = 0.0;
      p2 = pow(pm1[2]+pm1[3]*pow(p3,pm1[5])+pm1[6]*pow(pm1[7]-dy,
						       pm1[8]),pm1[9])/(1+exp(pm1[10]*(pm1[11]-pow(dy,pm1[12]))));
      if ( p1 > p2 ) {
	gr *= p2;
      } else {
	gr *= p1;
      }
      Real xla0 = 1/(1+ap[0]*y0*y0*y0);
      Real xla = 1/(1+ap[0]*y*y*y);
      gr *= xla0*(ap[1]+(1-ap[1])*(1+(exp(ap[2]*dy)-2)*sqrt(1-xla))
		  *exp(aic(y0,dy,ap[0])));

    } else {

      Real p1 = log(y0);
      Real p4 = 1/y0;
      pyx[4] = .230 + .033*p1 + .763*pow(y0,-0.507) + .007*pow(y0,-1.553);
      pyx[5] = 75.180*p4;
      pyx[6] = -.100 + .115*p1 + .100*pow(y0,1.623);
      if ( y0 < 1. ) {
	pyx[7] = pmx[18] + pmx[19]*p1;
      } else {
	pyx[7] = pmx[20] + pmx[21]*p1 + pmx[22]*pow((pmx[23]+p4),pmx[24]);
      }
      p1 = 1/y;
      gr = pyx[4]*pow((1+pyx[5]*p1),(pyx[6]*(1+pyx[7]*p1)));
      p1 = 1 + dy;
      gr *= pm1[17]*pow(dy,-1.5)*p1*pow((1+pm1[13]/p1),(pm1[14]+pm1[15]*pow(p1,pm1[16])));
      gr *= exp(aic(y0,(dy-.15),ap[0]));

    }

    grout[i] = gr;

  }

  return;
}

//****************************************************************************
Real aic(const Real y0, const Real dy, const Real apc)
{
  Real yp = y0 + dy;
  Real a0 = apc;
  Real a1 = pow(apc,0.3333);
  Real c0 = .5/a1;
  Real c1 = -3.4641*c0;
  Real c2 = 0.5774;
  Real d0l = a0*y0*y0*y0;
  Real d0h = a0*yp*yp*yp;
  Real d1l = a1*y0;
  Real d1h = a1*yp;
  Real ail = y0*(3-log(1+d0l)) + c0*log((1-d1l+d1l*d1l)/(1+d1l)/(1+d1l))
    + c1*atan(c2*(2*d1l-1));
  Real aih = yp*(3-log(1+d0h)) + c0*log((1-d1h+d1h*d1h)/(1+d1h)/(1+d1h))
    + c1*atan(c2*(2*d1h-1));

  return aih - ail;
}

//****************************************************************************
void pmxy(const Real Xm, RealArray& Pmx)
{
//  Angle dependent coefficients for Green's function for Compton
//  reflection;
//  xm  - cosine of the angle between normal to the slab and the
//        observation direction
//  pmx - output coefficients vector

  Real p2 = Xm*Xm;
  Real p3 = p2*Xm;
  Real p4 = p3*Xm;
  Pmx[0] = .8827 - .5885*Xm + .7988*p2 - .4117*p3;
  Pmx[1] = .0517 + .1076*Xm - .1691*p2 + .0678*p3;
  Pmx[2] = .0014 + .0043*Xm;
  Pmx[3] = -.0003 - .0027*Xm + .0021*p2;
  Pmx[4] = -1.3259 + 1.6214*Xm - 3.7007*p2 + 2.1722*p3 + exp(20.3*(Xm-1.063));
  Pmx[5] = .0790 - .4029*Xm + .6619*p2 + .2210*p3 - exp(17.7*(Xm-1.041));
  Pmx[6] = .0751 - .1848*Xm + .4068*p2 - .4126*p3 + exp(9.8*(Xm-1.185));
  Pmx[7] = -.0020 - .0394*Xm + .1004*p2 - .0597*p3;
  Pmx[8] = 3.2953 - 3.6996*Xm + 7.9837*p2 - 4.6855*p3;
  Pmx[9] = -.5278 + .9857*Xm - 1.7454*p2 - .3464*p3 + exp(27.6*(Xm-0.959));
  Pmx[10] = -.1919 + .5798*Xm - 1.2879*p2 + 1.4885*p3 - exp(18.3*(Xm-0.986));
  Pmx[11] = .0200 + .0832*Xm - .0333*p2 - .2370*p3 + exp(16.6*(Xm-1.086));
  Pmx[12] = -2.2779 + 2.4031*Xm - 4.0733*p2 + 1.9499*p3 - exp(19.6*(Xm-0.968));
  Pmx[13] = .4790 - 1.0166*Xm + 3.1727*p2 - 4.0108*p3 + 3.0545*p4 - exp(30.4*(Xm-0.957));
  Pmx[14] = .1122 - .3580*Xm + .4985*p2 - .3750*p3 - .5349*p4 + exp(31.2*(Xm-0.972));
  Pmx[15] = -.0160 - .0471*Xm - .0536*p2 + .2006*p3 + .2929*p4 - exp(30.6*(Xm-1.009));
  Pmx[16] = .0005 + .0021*Xm - .0447*p2 + .1749*p3 - .2303*p4 + exp(16.9*(Xm-1.130));
  Pmx[17] = .006 + .089*Xm - .102*p2 + .056*p3;
  Pmx[18] = .618 - .830*Xm;
  Pmx[19] = .128 - .132*Xm;
  Pmx[20] = .632 - .875*Xm;
  Pmx[21] = -.672 - .286*Xm + .717*p2 - .481*p3;
  Pmx[22] = .0126 - .0160*Xm + .0077*p2;
  Pmx[23] = .0111 + .0030*Xm - .0014*p2;
  Pmx[24] = -2.437 - .328*Xm - .260*p2 + .279*p3;

  return;

}

//****************************************************************************
void pm1y(const Real Xm, RealArray& Pm1)
{
  //  Angle dependent coefficients for Green's function for Compton
  //  reflection of incident photon at energy x0=1;
  //  xm  - cosine of the angle between normal to the slab and the
  //        observation direction
  //  pm1 - output coefficients vector

  Real p2 = .911 - .549*pow(1-Xm,1.471);
  Real p3 = .254 - .041*pow(Xm,-3.798);
  Pm1[0] = p2 - 2*p3;
  Pm1[1] = p3;
  Pm1[2] = .161 + .439*pow(Xm,.189) - .791*pow(Xm,7.789);
  Pm1[3] = -.871+1.740*pow(Xm,.176)-1.088*pow(Xm,.907);
  if ( Pm1[3] < 0.0 ) Pm1[3] = 0.0;
  Pm1[4] = .934 + .054*pow(Xm,-.666);
  Real p4 = Xm;
  if ( p4 > .524 ) p4 = .524;
  Pm1[5] = 4.647 - 11.574*p4 + 11.046*p4*p4;
  Pm1[6] = .012 + .199*pow(Xm,.939) + .861*pow(Xm,8.999);
  Pm1[7] = 2.150;
  p2 = Xm*Xm;
  p3 = p2*Xm;
  Pm1[8] = -1.281 + 2.374*Xm - 4.332*p2 + 2.630*p3;
  Pm1[9] = 2.882 + .035*pow(Xm,-1.065);
  Pm1[10] = 30.28 + 69.29*Xm;
  Pm1[11] = 1.037*pow(1-sqrt(1-.874*p2),.123);
  Pm1[12] = .123;
  p4 = p3*Xm;
  Pm1[13] = 56.50 - 27.54*Xm + 671.8*p2 - 1245.*p3 + 708.4*p4;
  Pm1[14] = -.897 - .826*Xm + 2.273*p2 - .984*p3;
  Pm1[15] = 1.240 - 1.297*Xm;
  Pm1[16] = -.490 + 1.181*Xm - 1.038*p2;
  Pm1[17] = 1.858*(pow(Xm,2.033)+0.745*pow(Xm,1.065));
  Pm1[18] = 1 - sqrt(1-(Xm-.05)*(Xm-.05));

  return;
}

//****************************************************************************
Real apf (const Real x)
{

#define c0 0.802
#define c1 -1.019
#define c2 2.528
#define c3 -3.198
#define c4 1.457
#define c4inv 0.030

  Real xinv = (c4inv/x);
  Real xinv2 = xinv*xinv;
  Real xinv4 = xinv2*xinv2;

  return c0 + x*(c1 + x*(c2 + x*(c3 + x*c4))) + xinv4;
}

//****************************************************************************
RealArray MZCompReflIntegrand(const RealArray& y0, void *p)
{
  return static_cast<MZCompRefl*>(p)->CalculateIntegrand(y0);
}
