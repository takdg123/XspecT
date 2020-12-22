#include<stdlib.h>
#include<stdio.h>
#include <string.h>
#include<math.h>
#include "numeric.h"

#define N 4000

double gammln(double xx);
double zeta(double xx);
double interp_funct(double *array_x, double *array_y, int dimen, double x);

double abs_val(double a);
double* boost(double *array_x, double *cflux);
double mod_bbnorm (double x, void* params);
double find_csi (double theta, double tau, double delta, double fb, double alpha1);

double h, a_norm, b_pow,  k_exp, q_fold, alpha_boost;

double * xsgrbcomp(double* energy, int Nflux, double* parameter, int spectrum, double*  flux, double* fluxError, char* init) {
  
  /* ############################################################################################# */

  double  x[Nflux+1], flux_bound[Nflux+1], unnorm_flux[Nflux], bb[Nflux+1], flux_bb[Nflux];

  double k1_norm, k2_norm,  k1_pow, k2_pow,  k1_exp,k1_fold, k2_fold, delta, delta_cor;
  double temp_ratio, nph_bb, params[2];

  double fb, fbflag, gamma, gamma2, kte, ktbb, bb_par,  tau, tau_cor, delta_ene,  z, cn, tmin, tmax;
  double A, beta,  theta, lambda_square, lambda_square2, csi, alpha1, norm_flux, nph_spec=0;

  int ie, ii, n;

  static double a[N+1], b[N+1], c[N+1], s[N+1], p[N+1], g[N+1], r[N+1];
  static double t[N+1], L[N+1], K[N+1], spec[N+1], xnum[N+1];

  /* ################################################################################ */
  /* Parameters of the model */
  /* ################################################################################ */

  ktbb=parameter[0];
  bb_par=parameter[1];
  kte=parameter[2];
  tau=parameter[3];
  beta=parameter[4];
  fbflag=parameter[5];
  A=pow(10,parameter[6]);
  z=parameter[7];
  alpha_boost=parameter[8];

  /* ################################################################################ */

  tmin=-12;
  tmax=12;
  h=(tmax-tmin)/(N);
  
  temp_ratio=kte/ktbb;
  fb=1+pow(beta,2)/(3*kte/511)*fbflag;

  /* ################################################################################ */

  theta=kte/511.;
  delta=-2./3*beta/(theta*tau);

  /* ########################################################################################### */
  /* Use the function for fitting the eigenvalue behaviour as a function */
  /* of the optical depth */
  /* lambda_square=a_norm(beta)*tau**[-b_exp(beta)]*Exp[-tau**[k_exp(beta)/q_fold(beta)] */
  /* and the parameters k1,k2,k3 and k4 in principle depends on  the albedo A */
  /* For GRBCOMP it is adopted the boundary condition A=0 */
  /* ########################################################################################### */

  k1_norm=4.05;
  k2_norm=-2.77;
  a_norm=k1_norm+k2_norm*beta;
  
  k1_pow=1.826;
  k2_pow=0.623;
  b_pow=k1_pow+k2_pow*beta;

  k1_exp=1.567;
  k_exp=k1_exp;

  k1_fold=3.19;
  k2_fold=-1.44;

  q_fold=k1_fold*pow(beta,k2_fold);

  lambda_square=a_norm*pow(tau, -b_pow)*exp(-pow(tau,k_exp)/q_fold);
  gamma=lambda_square/theta;


  /* Value of the spectral index for the case in which the optical depth is not corrected */
  /* for the outflow velocity, a factor of order of (1-beta) */

  alpha1=-(delta+3)/2. +sqrt(pow(delta-3,2)/4.+ gamma); 
  
  /* ########################################################################################################################## */
  /* In the presence of outflow bulk motion, the effective optical depth is less than that of the static case */
  /* by a factor of order (1-beta) */
  /* When considering the second-order bulk Comptonization effect term fb, the spectral index alpha becomes lower, */
  /* but this effect should be partially compensated by decreased value of the optical depth, which leads */
  /* instead to higher alpha values */
  /* The estimation is done by finding a numerical value of the bulk term such that value of the spectral index with */
  /* the fb-term is similar to the case with no fb-term and tau equal to the static value */
  /* ########################################################################################################################## */

  if (fbflag==1) {

    csi=find_csi (kte/511, tau, delta, fb, alpha1); 

    tau_cor=tau*csi;
    delta_cor=delta/csi;
    lambda_square2=a_norm*pow(tau_cor, -b_pow)*exp(-pow(tau_cor,k_exp)/q_fold);
    gamma2=lambda_square2/theta;


    tau=tau_cor;
    delta=delta_cor;
    gamma=gamma2;

  }

  /* ########################################################################## */
  /* Coefficients for the sweep method */
  /* ########################################################################## */
 
  for (n=0; n<=N; n++) {

    t[n]=tmin+n*h;
    xnum[n]=exp(t[n]);

    p[n]=fb;
    g[n]=(exp(t[n])-3*fb-delta);
    r[n]=(exp(t[n])-gamma+3*delta);



    s[n]=-exp(bb_par*t[n])/(exp(temp_ratio*exp(t[n]))-1);

    if (s[n] > -pow(10,-15)) s[n]=0; 

    a[n]=p[n]/pow(h,2);
    b[n]=-(2.*p[n]/pow(h,2)+g[n]/h-r[n]); 
    c[n]=(p[n]/pow(h,2)+g[n]/h);

  }


  /* ########################################################################## */
  /* Boundary conditions */
  /* ########################################################################## */


  a[0]=0; b[0]=1; c[0]=0; s[0]=0;
  a[N]=0; b[N]=1; c[N]=0; s[N]=0;


  L[0]=-c[0]/b[0];
  K[0]=0;					       

  /* ########################################################################## */

  for (n=1; n<N; n++) {
    L[n]=-c[n]/(a[n]* L[n-1] + b[n]);
    K[n]=(s[n]- a[n]*K[n-1])/(a[n]* L[n-1] + b[n]);
  }

  /* ########################################################################## */
  /* Boundary conditions for the intensity */
  /* ########################################################################## */

  spec[N]=0;
  spec[0]=0;


  for (n=N; n>0; n--) {
    spec[n-1]=L[n-1]*spec[n]+K[n-1]; 
    nph_spec=nph_spec+spec[n-1]*h;

    if (spec[n-1] < 0) {
      printf("Warning: with these parameters flux is negative\n\n"); 
    }

  }

 
  /* ############################################################### */

  double* boost_flux = (double*)malloc((N+1)*sizeof(double)); 
  boost_flux=boost(t, spec);

  params[0]=bb_par+1;
  params[1]=kte/ktbb;

  for (ie = 0; ie <=Nflux; ie++) {

    ii=ie-1;
    x[ie]=(1+z)*(energy[ie])/kte;


    flux_bound[ie]=interp_funct (xnum, boost_flux, N+1, x[ie]);
    flux_bound[ie]=flux_bound[ie]/energy[ie];

    /* Soft component thermal */

    bb[ie]=mod_bbnorm(x[ie], params)/energy[ie];

    if (ie > 0) {

      delta_ene=(energy[ie]-energy[ie-1]);
      unnorm_flux[ii]=1./(1+z)*0.5*(flux_bound[ie]+flux_bound[ie-1])*delta_ene;
      flux_bb[ii]=1./(1+z)*0.5*(bb[ie] +  bb[ie-1])*delta_ene;

    }

  }


  nph_bb=pow(kte/ktbb,-bb_par)*exp(gammln(bb_par))*zeta(bb_par);


  /* ################################################################################ */
  /* Normalization constant of the BB with radius in units of 10^9 cm  */
  /* and distance in units of Mpc    */
  /* ################################################################################ */
  
  cn=10.4;

  /* ########################################################################## */
  /* Numerical constant to allow conservation of the photon number */
  /* ########################################################################## */

  norm_flux=nph_bb/nph_spec; 

  for (ii=0; ii < Nflux; ii++) {
    flux[ii]=cn/(A+1)*pow(kte,bb_par)*(flux_bb[ii] + A*norm_flux*unnorm_flux[ii]);

  }



  return(flux);


}

/* ####################################################################*/
/* Second part of the model: solution of the radiative */
/* transfer equation takeing into account only the up-scattering term  */
/* ####################################################################*/

double*  boost  (double *t, double *cflux) {

  static double a[N+1], b[N+1], c[N+1], s[N+1], p[N+1], g[N+1], r[N+1];
  static double K[N+1], L[N+1], bspec[N+1];

  double gamma_boost=alpha_boost*(alpha_boost+3);

  int n;

  for (n=0; n<=N; n++) {


    p[n]=1;
    g[n]=-3;
    r[n]=-gamma_boost;

    a[n]=p[n]/pow(h,2);
    b[n]=-(2.*p[n]/pow(h,2)+g[n]/h-r[n]); 
    c[n]=(p[n]/pow(h,2)+g[n]/h);

  }

  /* ########################################################################## */
  /* Boundary conditions */
  /* ########################################################################## */

  a[0]=0; b[0]=1; c[0]=0; s[0]=0;
  a[N]=0; b[N]=1; c[N]=0; s[N]=0;


  K[0]=0;
  L[0]=-c[0]/b[0];

  /* ########################################################################## */

  for (n=1; n<N; n++) {

    s[n]=-cflux[n];

    L[n]=-c[n]/(a[n]* L[n-1] + b[n]);
    K[n]=(s[n]- a[n]*K[n-1])/(a[n]* L[n-1] + b[n]);

  }

  bspec[N]=0;

  /* ########################################################################## */

  for (n=N; n>0; n--) {
    bspec[n-1]=((L[n-1]*bspec[n]+K[n-1])); 

    if (bspec[n-1] < 0) {
      printf("Warning: with these parameters flux is negative\n\n");  
    }


  }

  for (n=0; n <= N; n++) {
    bspec[n]=gamma_boost*bspec[n];
  }

  return(bspec);

}


/* ######################################################################################## */

double mod_bbnorm (double x, void *params) {

  double bb_par, ratio_temp, value;

  bb_par=((double *)params)[0];
  ratio_temp=((double *)params)[1];

  value=pow(x,bb_par-1)/(exp(ratio_temp*x)-1);

  return(value);

}


/* ######################################################################################## */

double find_csi (double theta, double tau, double delta, double fb, double alpha1) {

  double tau_cor, delta_cor, gamma2, alpha2, csi, lambda_square2, funct;
  double csi_bound[2], csilow, csiup, csinew, csinew_previous, f_csilow, f_csiup;
  int i,j;

  csi_bound[0]=0.1;
  csi_bound[1]=1;

  csinew=csi_bound[0];
  csinew_previous=csi_bound[1];


  for (i=0; i<=100; i++) {

    if (abs_val(csinew-csinew_previous) < 0.001) break; 

    for (j=0; j<=1; j++)  {

      csi=csi_bound[j]; 

      tau_cor=tau*csi;
      delta_cor=delta/csi;
      lambda_square2=a_norm*pow(tau_cor, -b_pow)*exp(-pow(tau_cor,k_exp)/q_fold);
      gamma2=lambda_square2/theta;

      alpha2=(-3*fb-delta_cor+sqrt(pow(3*fb-delta_cor,2)+4*fb*gamma2))/(2*fb);

      funct=alpha2-alpha1;

      /* ################################################################################### */

      if (j==0) {
	f_csilow=funct;
	csilow=csi_bound[0];
 
      } else {

	f_csiup=funct;
	csiup=csi_bound[1];
 
      }

    }

    /* ################################################################################### */
    /* Determine the intersection point with the X-axis */
    /* ################################################################################### */

    if (f_csilow/f_csiup < 0) {

      csinew_previous=csinew;
      csinew=csiup-(csiup-csilow)/(f_csiup-f_csilow)*f_csiup;

    }


    /* End of loop over i */

    csi_bound[1]=csinew;

  }


  return(csinew);

}
		 

/* ######################################################################################## */

double abs_val(double a) {

  if (a > 0) {
    return(a);
  } else {
    return(-a);
  }

}



