/*
 *   LIBRARIES 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <XSFunctions/Utilities/funcType.h>
#include <xsTypes.h>
//#include <cmath>
//#include <iomanip>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSUtil/Numerics/CosmologyFunction.h>
#include <XSUtil/Utils/XSstream.h>
#include <sstream>
#include <memory>
#include <XSstreams.h>

extern "C" {
 /*
  *                          NEW STRUCTURES AND FUNCTIONS
  */
struct ABUN {
  double   cc, xx, rr;
     };   
struct HYDR {
  double cc, ff, cx, cr, gx, gr;
};

struct  TEMP {
  double  cc, dt, ix, ir, cx, cr, tx, tr;
};
/*
 *                               DECLARATIONS
 */
double aa(double r, struct ABUN Ab);
double dei(double r, double innsq, double outsq, struct HYDR nH);
double eint(double r1,double r2,double outsq, double innsq,struct HYDR nH);
double hh(double r, struct HYDR nH);
double tt(double r, struct TEMP kT);
double DSQR(double a);

//      SUBROUTINE sumdem(itype, switch, ear, ne, abun, dens, z,
//     &                  ninputt, inputt, dem, ifl, 
//                      qtherm, velocity, photar, status)
//
//      INTEGER itype, switch, ne, ninputt, ifl, status
//      REAL ear(0:*), abun(*), inputt(*), dem(*), photar(*)
//      REAL dens, z

void sumdem_(int& itype, int& flag, float* ear, int& ne,
                    float* abun, float& dens, float& z, int& ninputt,
                    float* inputt, float* dem, int& ifl, int& qtherm,
                    float& velocity, float* photar, int& status);
}

/*---------------------------------------------------------------------------*/

double aa(double r, struct ABUN Ab)
     /*
      *    PARAMETERS
      * 
      *    r   i :   radius [Mpc]
      *    Ab  i :   structure of the metal distributon parameters
      *    aa  r :   metal abundance [solar units]
      */
{
  return (  Ab.cc / pow( (1.0 + DSQR(r/Ab.rr) ), Ab.xx ) );
}

/*---------------------------------------------------------------------------*/

double dei(double r, double innsq, double outsq, struct HYDR nH)   
     /*
      *  PARAMETERS 
      *  input 
      *      r        i :  radius [Mpc]
      *      innsq    i :  squared inner rim of the projected annulus  [Mpc**2]
      *      outsq    i :  squared outer rim of the projected annulus  [Mpc**2]
      *      nH       i :  structure of the hydrogen distribution
      *      dei      r :  non-normalised differential emission integral  
      *                    of the region bounded   by the projected ring  [inner, outer]
      *                    the actual differential emission integral is obtained by 
      *                    multiplying dei  by the correction factor for electron density 
      *                    <elden>~1.21 and by the  geometrical angular factor <angfac>
      */     
{
  double k, r_sq;
  r_sq = DSQR(r); 
  k    = sqrt(r_sq - innsq);
  if ( outsq < r_sq )
    k -= sqrt(r_sq - outsq);
  return ( r * k * DSQR( hh(r,nH) ) );
}

/*---------------------------------------------------------------------------*/

double eint(double r1,double r2,double outsq,double innsq, struct HYDR nH)
     /*
      * PARAMETERS
      *      r1     i :  inner radius of the spherical shell
      *      r2     i :  outer radius of the spherical shell
      *      innsq  i :  squared inner rim of the projected annulus
      *      outsq  i :  squared outer rim of the projected annulus
      *      nH     i :  structure of the hydrogen distribution
      *      eint   r :  non-normalised  emission integral of the region
      *                  bounded by the projected ring [inner, outer]
      *                  and comprised within the spherical shell [r1, r2] 
      *                  the actual  emission integral is obtained by 
      *                  multiplying <eint>   by the correction factor for electron density 
      *                  <elden>~1.21 and by the  geometrical angular factor <angfac>
      *
      * NOTE          
      *                  <eint> is calculated  with a 10-points Gauss-Legendre formula
      * 
      */     
{
  double d_m, d_p, xr,xm,dx,s;
  static double x[]={0.0,0.1488743389,0.4333953941,
		     0.6794095682,0.8650633666,0.9739065285};
  static double w[]={0.0,0.2955242247,0.2692667193,
		     0.2190863625,0.1494513491,0.0666713443};
  xm=0.5*(r2+r1);
  xr=0.5*(r2-r1);
  s=0.0;
  for (size_t j=1;j<=5;j++) {
    dx=xr*x[j];
    d_m = dei(xm-dx, innsq, outsq, nH);
    d_p = dei(xm+dx, innsq, outsq, nH);
    s += w[j]*(d_m + d_p);
  }
  return s *= xr;
}

/*---------------------------------------------------------------------------*/

double hh(double r, struct HYDR nH)
     /*
      * PARAMETERS
      *
      *    r    i :  radius [Mpc]
      *    nH   i :  structure with the hydrogen distribution parameters
      *    hh   r :  hydrogen density [cm**-3]
      *
      * LIST AND MEANING OF THE MAIN VARIABLES
      *
      *  h_c   pure ICM  component
      *  h_g   dominant  galaxy component (if any)
      */    
{ 
  double   h_c, h_g;
  h_c =  nH.ff          /  pow( 1.0 + DSQR(r/nH.cr) , nH.cx);
  h_g =  (1.0 - nH.ff)  /  pow( 1.0 + DSQR(r/nH.gr) , nH.gx );
  return ( nH.cc * (h_c + h_g) );
}

/*---------------------------------------------------------------------------*/

double tt(double r, struct TEMP kT)
     /*
      *   PARAMETERS
      *
      *    r   i :  radius [Mpc]
      *    kT  i :  structure of  the temperature parameters
      *    tt  r :  electronic   temperature [keV]
      */
{  
  double ans, k, l, rcx,rix;
  if (r==0)
    return kT.cc;
  else
    {
      k = pow( 1.0 +  DSQR(r/kT.tr), -kT.tx );
      if (kT.dt==0.0)
	ans  = k * kT.cc;
      else    
	{ 
	  rcx  = pow(kT.cr / r, kT.cx);
	  rix  = pow(r / kT.ir, kT.ix);
	  l    = M_2_PI * atan(rix) / (1 + rcx);
	  ans  = k * (  kT.cc + l * kT.dt );
	}
      return ans;
    }
}

/*---------------------------------------------------------------------------*/
double DSQR(double a)
{
  return (a = 0.0 ? 0.0 : a*a);
}

/*---------------------------------------------------------------------------*/

/*	
 * XSPEC subroutine for the analytical deprojection of an extended 
 * and optically thin  source. This  model is fed  with a set of spectra 
 * extracted in annular sectors, concentric about the X-ray emission peak.  
 * Spherical symmetry is assumed. The mode parameters define the 3D 
 * distributions of hydrogen density, temperature and metal abundance, 
 * as well as the redshift and other options (see the parameter list below).
 * For each annular sector, the inner boundary (in arcmin), the outer boundary
 * (also in arcmin) and the sector width (in degrees) are  specified 
 * (respectively) by three XFLTnnnn keywords identified by the keys inner,
 * outer and width, respectively. For backwards compatibility, if these are not
 * found then it will try XFLT0001, XFLT0002, and  XFLT0003. These must be 
 * added to the spectrum extension  in each input file  (e.g. with the FKEYPAR 
 * ftool).  
 * If the interactive chattiness level in XSPEC is set to a value > 10, 
 * smaug also prints:
 *
 *  - H0   =  Hubble constant [km/s/Mpc]
 *  - q0   = deceleration parameter
 *  - L0   = cosmological constant
 *  - DA   = source angular distance [Mpc]
 *  - DSET = dataset no. to which the quantities listed below are referred
 *  - IN   = inner rim of the projected  annular sector [Mpc]
 *  - OUT  = outer rim of the projected  annular sector [Mpc]
 *  - WID  = width of the  projected  annular sector  [deg]
 *  - EVOL = emitting volume  within the integration radius cutoff  [Mpc**3]
 *  - EINT = emission integral  within the integration 
 *           radius cutoff [Mpc**3 cm**-6].If nH.cc is frozen to 1, 
 *           the actual EI is obtaned   by multiplying this figure 
 *           by the square root of the model normalisation
 *           
 *   
 *  AUTHOR       - Fabio Pizzolato
 *
 *  LAST REVISED - 30th July  2003
 *
 *  ADDITIONAL NOTES
 *
 *        - the cosmological parameters H0 (Hubble's constant), q0 (deceleration parameter) and
 *           L0 (cosmological constant) must be set  with the XSPEC < cosmo > command
 *           before running smaug
 *        - He fixed to solar, metals frozen together
 *        - the output model spectrum is multiplied by the constant <norm>  
 *           defined below, so that the normalization  imposed by XSPEC is 
 *           the central hydrogen density squared, in [cm**-6],  provided  
 *           that  the fit paramter  nH.cc  is frozen to 1.0.	  
 *
 *  ARGUMENTS
 *              
 *
 *    ear        i:   energy ranges  [keV]
 *    ne         i:   number of elements in photar array
 *    param      i:   array of model parameters
 *                           1  kT.cc   central temperature [keV]
 *                           2  kT.dt   max difference of temperature [keV]
 *                           3  kT.ix   exponent of the inner temperature
 *                           4  kT.ir   radius of the inner temperature [Mpc] 
 *                           5  kT.cx   exponent of the middle temperature
 *                           6  kT.cr   radius of the middle temperature [Mpc]  
 *                           7  kT.tx   exponent of the outer temperature
 *                           8  kT.tr   radius of the outer temperature [Mpc] 
 *                           9  nH.cc   central hydrogen density [cm**-3]
 *                          10  nH.ff   fraction of nH.cc  relative 
 *                                       to the 1st beta component
 *                          11  nH.cx   exponent of the first beta component
 *                          12  nH.cr   radius of the 1st beta  component [Mpc] 
 *                          13  nH.gx   exponent of the 2nd beta component
 *                          14  nH.gr   radius of the 2nd beta  component [Mpc] 
 *                          15  Ab.cc   central metallicity [solar units]
 *                          16  Ab.xx   exponent of the metal distribution 
 *                          17  Ab.rr   radius of the  metal distribution [Mpc]
 *                          18  z       redshift of the source        
 *                          19  meshpts integration mesh-points
 *                          20  rcutoff integration cutoff radius [Mpc]	  
 *                          21  mode    mode of spectral evaluation:
 *                                              0 = calculate 
 *	                                        1 = interpolate 
 *                                              2 = APEC interpolate
 *                          22  itype     type of plasma emission code
 *                                              1 = Raymond-Smith
 *                                              2 = Mekal
 *                                              3 = Meka
 *                                              4 = APEC 
 *                          23  norm      model normalisation 
 *                                        (squared central H density [cm**-6])
 *	
 *    ifl        i:  dataset number
 *    photar[ ]  r:  calculated model spectrum [photons/cm**2/s/bin]
 *    photer[ ]  r:  model spectrum uncertainties  (unused)
 *
 *
 * LIST AND MEANING OF THE MAIN VARIABLES 
 *
 *     Ab            : structure defining the heavy elements profile
 *     abund[ ]      : array of metallicities (one slot per element) 
 *                      passed to the routine SUMDEM
 *     angfac        : angular factor of the spectrum extraction region 
 *     da            : angular distance of the source  [Mpc]
 *     elden         : conversion factor hydrogen -> electron density 
 *                      (Anders & Grevesse abundances)
 *     ei            : emission integral
 *     ei_shell      : emission integral of each integration  shell,
 *                   : corrected for the electron density <elden> 
 *                      but not for the geometrical angular factor <angfac> 
 *     H0            : Hubble constant [km/s/Mpc]
 *     h1 (h2)       : integration step for r < outer  (resp. outer < r)
 *     inner         : inner rim of the spectrum extraction region
 *     kT            : structure of the (electronic) temperature profile
 *     kT_shell[ ]   : one-element array with the average temperature
 *                   : of an integration shell
 *     L0            : cosmological constant 
 *     nH            : structure of  the H profile
 *     outer         : outer rim of the spectrum extraction region
 *     phoshell[ ]   : non normalised photon array  from one integration shell
 *     q0            : deceleration paramter
 */

/*
 *                                             MAIN ROUTINE                 
 */

extern "C" 
void xsmaug(Real *ear, int ne, Real *param, int ifl, 
                Real *photar, Real *photer, const char* initstring)
{  

  struct ABUN  Ab;   
  struct HYDR  nH;
  struct TEMP  kT;
  const double MIN2RAD =  3437.75;
  int    i, ie, iel; 
  int    mesh2, meshpt;
  int    itype, mode, no_el, status;
  double angfac, da, ei, evol, inner, H0, L0, outer, q0, z, zfac;
  double a1, a2, Ab_shell, ei_shell,  h1, h2,  r1, r2, rcutoff, t1, t2,  w1, w2, w12;
  double elden, norm;
  float *abund,  *phoshell;            /* dynamic arrays */
  float kT_shell[1], dem[ ] = {1.0};   /* static arrays  */

// forcecalc not implemented in xspec12
// check that the force-calc logical is set to ensure that this routine is
//   called for every dataset 
//
 // if ( !getForceCalc() ) {
 //   errStrg = (char*)malloc(256*sizeof(char));
 //   strcpy(errStrg, "Doing an xset forcecalc...");
 //   XWRITE(errStrg, 10);
 //   setForceCalc(!getForceCalc());
 //   free(errStrg);
 // }

  /*   INPUT MODEL  PARAMETERS   */

  kT.cc     = param[0]; 
  kT.dt     = param[1]; 
  kT.ix     = param[2];
  kT.ir     = param[3];
  kT.cx     = param[4]; 
  kT.cr     = param[5];
  kT.tx     = param[6]; 
  kT.tr     = param[7];
  nH.cc     = param[8];
  nH.ff     = param[9];
  nH.cx     = param[10];
  nH.cr     = param[11];
  nH.gx     = param[12];
  nH.gr     = param[13];
  Ab.cc     = param[14];
  Ab.xx     = param[15]; 
  Ab.rr     = param[16];
  z         = param[17];
  meshpt    = (int) param[18];
  rcutoff   = param[19];
  mode      = (int) param[20];
  itype     = (int) param[21];


  /*   GET THE COSMOLOGICAL PARAMETERS   */

  q0 = (double)FunctionUtility::getq0();
  H0 = (double)FunctionUtility::getH0();
  L0 = (double)FunctionUtility::getlambda0();

  /*    SET THE GEOMETRICAL FACTORS  */

  Numerics::FZSQ fzsq;
  zfac     = 1.0 + z;
  da = (2.9979E+05 / H0) * sqrt(  fzsq(z ,q0 , L0) ) / DSQR(zfac);
  mesh2    = 2 * meshpt; 

  if ( FunctionUtility::inXFLT(ifl, "inner") ) {
    inner = (Real)FunctionUtility::getXFLT(ifl, "inner");
  } else if ( FunctionUtility::inXFLT(ifl, 1) ) {
    inner = (Real)FunctionUtility::getXFLT(ifl, 1);
  } else {
    std::ostringstream errStrg;
    errStrg <<  "XSmaug: cannot find XFLTnnnn keyword for inner annulus for spectrum "  << ifl << "\n";
    throw FunctionUtility::FunctionException(errStrg.str());
  }

  if ( FunctionUtility::inXFLT(ifl, "outer") ) {
    outer = (Real)FunctionUtility::getXFLT(ifl, "outer");
  } else if ( FunctionUtility::inXFLT(ifl, 2) ) {
    outer = (Real)FunctionUtility::getXFLT(ifl, 2);
  } else {
    std::ostringstream errStrg;
    errStrg <<  "XSmaug: cannot find XFLTnnnn keyword for outer annulus for spectrum "  << ifl << "\n";
    throw FunctionUtility::FunctionException(errStrg.str());
  }

  if ( FunctionUtility::inXFLT(ifl, "width") ) {
    angfac = (Real)FunctionUtility::getXFLT(ifl, "width");
  } else if ( FunctionUtility::inXFLT(ifl, 3) ) {
    angfac = (Real)FunctionUtility::getXFLT(ifl, 3);
  } else {
    std::ostringstream errStrg;
    errStrg <<  "XSmaug: cannot find XFLTnnnn keyword for width for spectrum "  << ifl << "\n";
    throw FunctionUtility::FunctionException(errStrg.str());
  }

  inner  *= da / MIN2RAD;     
  outer  *= da / MIN2RAD;     
  angfac *= 4.0 * M_PI/360.0;   
  if ( (rcutoff <= outer) || (outer <= inner) || (angfac ==0.0 )) {
    std::ostringstream errStrg;
    errStrg << "XSmaug: for of dataset " << ifl 
	    << " either the outer ring exceeds the cutoff radius, the outer ring"
	    << " is less than or equal to the inner, or the sector width is zero\n";
    throw FunctionUtility::FunctionException(errStrg.str());
  }

  /*
   * ALLOCATE MEMORY FOR THE ARRAY OF THE METAL ABUNDANCES  
   * AND FOR THE SPECTRUM EMITTED BY EACH INTEGRATION SHELL
   */

  switch (itype)
    {
    case 1:            
      no_el   = 12;     break;  /* 12 elements for  Raymond-Smith  */
    case 2: case 3:            
      no_el   = 14;     break;  /* 14 elements for  Meka and mekal */
    case 4:                
      no_el   = 13;     break;  /*  13 elements for APEC  code     */
    default:
      {
                throw FunctionUtility::FunctionException
                   ("XSMaug: unknown  number of elements");
      }
    }  

  abund     =  (float*) malloc( no_el * sizeof(float));
  phoshell  =  (float*)  malloc( ne   * sizeof(float)); 


  for(ie=0; ie < ne; photar[ie] = 0.0, ++ie)
    ;  /* reset the output photon array */

    h1 = (outer - inner)/DSQR(meshpt); 
    h2 = (rcutoff - outer)/DSQR(meshpt);      

    r1=inner;
    t1=tt(r1,kT);

    a1=aa(r1,Ab);
    w1 = 0.0;

 /*   INTEGRATION LOOP    */       

  for(i=1, ei = 0.0; i<=mesh2; ++i)          
    {  
      r2     = ((i<= meshpt) ?  \
		(inner+h1*DSQR(i))  : (outer+h2*DSQR(i-meshpt)));
      w2   = dei(r2, DSQR(inner), DSQR(outer), nH);
      w12 = w1 + w2;      

      t2  = tt(r2,kT);
      kT_shell[0]  = (w1 * t1 + w2 *t2 )/w12;  

      a2  = aa(r2,Ab);
      Ab_shell      = (w1 * a1 + w2 * a2 )/w12;
      for(iel=1, abund[0]=1.0; iel < no_el; ++iel)
	abund[iel] = (float) Ab_shell;   

      /*   CALCULATE THE SPECTRUM OF EACH INTEGRATION  SHELL  */       

      elden    = 1.195 + 1.363E-02 * Ab_shell;
      ei_shell = elden * eint(r1, r2, DSQR(outer), DSQR(inner), nH); 
      std::unique_ptr<float[]> pFear(new float[ne+1]);
      float* fear(pFear.get());
      for (size_t l = 0; l <= (size_t)ne; ++l) fear[l] = ear[l];
      static float ONE(1.);
      static int one(1);
      float fz (z);
      static int NOT(0);
      static float ZERO(0.0);
      sumdem_(itype, mode, fear, ne, abund, ONE , fz, one, kT_shell, dem, ifl, 
                NOT, ZERO, phoshell, status);

      /* ADD THE SHELL CONTRIBUTION TO THE OUTPUT SPECTRUM, 
       * CORRECTING IT FOR <angfac>, DISTANCE AND REDSHIFT, 
       * AND FOR A FINAL FACTOR GIVING THE CORRECT XSPEC NORMALISATION
       */

      norm = 2.456E+09 * ei_shell * angfac / DSQR(zfac*da);     

      for(ie=0; ie<ne; ++ie)
	photar[ie] +=  norm * phoshell[ie];

      /*
       * ADD THE SHELL CONTRIBUTION TO THE THE DATASET EMISSION INTEGRAL
       * AND UPDATE FOR THE NEXT INTEGRATION LOOP
       */

      ei  += ei_shell;     
      r1  = r2;
      w1  = w2;
      t1  = t2;
      a1 =  a2; 
    }  

  /*   PRINT SOME USEFUL QUANTITIES IF THE CHATTINESS LEVEL IS > 15  */

  using namespace std;
  ostringstream oss;
  evol = angfac * 
    ( pow(DSQR(rcutoff)-DSQR(inner), 1.5) -  
      pow(DSQR(rcutoff)-DSQR(outer), 1.5 ))/3.0;

  if ( ifl == 1 ) {
    oss << "\n H0 = " << fixed << H0 << "  [km/s/Mpc]"
	<<"\nq0   = " << q0
	<< "\nL0   = " << L0 ;
    oss.precision(3);
    oss << "\nDA   =  " <<  scientific << da << " [Mpc]" << "\n";
  }

  oss << "\nDSET " << ifl
      <<  "\nIN   =  " << inner << " [Mpc] "
      <<  "\nOUT  =  " << outer << " [Mpc] ";
  oss.precision(1);
  oss << "\nWID  = " << fixed 
      <<  angfac << " [rad]" ;
  oss << "\nEVOL (< " << rcutoff << " [Mpc] = ";
  oss.precision(3);
  oss << scientific << evol << "[Mpc**3]         ";
  oss.precision(1);
  oss << "\nEINT (<" << fixed << rcutoff << " [Mpc] = ";
  oss.precision(3);
  xs_write(const_cast<char*>(oss.str().c_str()),15);


  /*  
   *        DISPOSE MEMORY 
   */
  free(abund);
  free(phoshell);
}
