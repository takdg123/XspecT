/* #################################################################################### */
/* Numerical routine for the solution of the Fokker-Planck approximation */
/* of the radiative transfer equation with cylindrical simmetry */
/* using an iteration method for partial differential equations */
/* Both the energy and space terms of the equation are considered */
/* #################################################################################### */

#include "numeric.h"

#define v_light 3*1e10
#define e_kin_fact 9.1094*1e-28* pow(v_light,2)*6.24*1e8
#define M 1000
#define bulkflag 2

/* ################################################################### */

double bb2d (double R, double x);
double slope(double energies[], double** spectrum, int qmin, double h, double ktbb, double kte, int NQ);
double check_conv(double previous_index,double index);

/* ################################################################### */

double * xscompmag(double* energy, int Nflux, double* parameter, int spectrum, double*  flux, double* fluxError, char* init) {

  double h, htau,  delta_ene,cn, qmin, qmax, taumin, taumax, r0, ht, norm_low, norm_int, norm_up, beta_max;
  double ktbb, kte, mdot, ratio_temp, theta, delta_hat, A, albedo,  eps, deriv_vel;
  double W, Z, P,Q, R_hat, G, H, csi, eta, zmax, aa;
  double previous_index, index=1.0, ztau;
  
/* Neutron star radius in Schwarzschild units*/
  double  z0=2.42;  

 /* Constant for the z(tau) formula */
  double k_zt=0.00217; 

  int i, j, ie, ii, m, count_conv=0, iter_val=0, betaflag, maxiter=0;
  int NTAU, NQ=200, ntaumin;

  /* matrices */
  double **a, **b, **c, **d, **e, **f, **s, **S_hat, **L, **K;
  double **L_tilde, **K_tilde, **S_tilde;
  double **spec_low, **spec_int, **spec_up, **flux_bulk;

  /* vectors */
  double *flux_comp, *x, *flux_bb, *q, *xnum, *flux_tot;
  double *vel, *tau, *phtot, *bbnorm, *integ_int, *integ_low, *integ_up, *fact;

 
/* ##########################################################*/
/* Model parameters */
/* ##########################################################*/

  ktbb=parameter[0];
  kte=parameter[1];
  taumax=parameter[2];
  eta=parameter[3];
  beta_max=parameter[4];
  r0=parameter[5];
  albedo=parameter[6];
  betaflag=parameter[7];

/* #############################################################*/

  H=(kte/511.)*(0.1/1e-3);
  taumin=0;
  zmax=2*z0;
  A=-beta_max*pow(z0,eta); 

/* #############################################################*/
/*Determine mdot from the continuity equation for both*/
/*velocity profiles*/
/*#############################################################*/

  if (betaflag==1) {

    mdot=-taumax*A*pow(r0,2)*(1+eta)/(0.002174*(pow(zmax,1+eta)-pow(z0,1+eta)));

  } else {

    mdot=50*pow(r0,3./2)/(pow(z0,1./2) *pow(zmax-z0, 1./2))*taumax;

  }

/* ############################################################*/

  qmin=-10;
  qmax=10;
  h=(qmax-qmin)/(NQ);


 /* ######################################################## */

  ratio_temp=kte/ktbb;
  theta=kte/511.;
 
/* ########################################################################## */

  G=(3/2.)*(1-albedo)/(1+albedo); 

/* #################################################################################### */
/* When using the definition J(x,tau)=R(tau) x^{-alpha} the inner boundary condition*/
/* becomes u[i][0]=L_tilde[i][0]* u[i][1] where  */
/* L_tilde[i][0]=1/(1+htau(-beta_max*(alpha+3)+G)); */
/* The necessary condition is L_tilde[i][0] > 0 */
/* #################################################################################### */
 
  if (-beta_max*(index+3)+G > 0) {
    htau=h/sqrt(H); 
    htau=0.5/(-beta_max*(index+3)+G);  

  } else {

    htau=-0.5/(-beta_max*(index+3)+G);

  }

  NTAU=taumax/htau; 


/* #################################################################### */

  ntaumin=10;

  if (NTAU < ntaumin) {

    double phi;
    phi=(taumax*(-beta_max*(index+3)+G))/ntaumin;
    htau=phi/(-beta_max*(index+3)+G);

    NTAU=taumax/htau; 
    
  }

/* #################################################################### */

  a=dmatrix(0,NQ,0,NTAU);
  b=dmatrix(0,NQ,0,NTAU);
  c=dmatrix(0,NQ,0,NTAU);
  d=dmatrix(0,NQ,0,NTAU);
  e=dmatrix(0,NQ,0,NTAU);
  f=dmatrix(0,NQ,0,NTAU);
  s=dmatrix(0,NQ,0,NTAU);
  S_hat=dmatrix(0,NQ,0,NTAU);
  L=dmatrix(0,NQ,0,NTAU);
  K=dmatrix(0,NQ,0,NTAU);
  L_tilde=dmatrix(0,NQ,0,NTAU);
  K_tilde=dmatrix(0,NQ,0,NTAU);
  S_tilde=dmatrix(0,NQ,0,NTAU);
  spec_low=dmatrix(0,NQ,0,NTAU);
  spec_int=dmatrix(0,NQ,0,NTAU);
  spec_up=dmatrix(0,NQ,0,NTAU);
  flux_bulk=dmatrix(0,NQ,0,NTAU);

/* #################################################################### */
 
  flux_comp=dvector(0,Nflux); 
  x=dvector(0,Nflux); 
  flux_bb=dvector(0,Nflux);
  q=dvector(0,NQ);
  xnum=dvector(0,NQ); 
  flux_tot=dvector(0,NQ);

  vel=dvector(0,NTAU);
  tau=dvector(0,NTAU);
  phtot=dvector(0,NTAU);
  bbnorm=dvector(0,NTAU);
  integ_int=dvector(0,NTAU);
  integ_low=dvector(0,NTAU);
  integ_up=dvector(0,NTAU);
  fact=dvector(0,NTAU);


/* ############################################################# */
/* Define the step-size for the time variable */
/* ############################################################# */

  ht=0.1*3*H*htau*htau;
 
/* ########################################################################## */
/* Velocity profile*/
/* ########################################################################## */

  csi=(15.77*r0)/mdot;
  aa=0.67/z0*csi;

  for (j=0; j<=NTAU; j++) {

    tau[j]=taumin+j*htau;
    ztau=pow((mdot*k_zt*pow(z0,1+eta)-tau[j]*A*pow(r0,2)*(1+eta))/(mdot*k_zt),1/(1.+eta));

/* ###################################################################### */
/* First velocity profile */
/* v(z)=-A*z^{-eta} */
/* ###################################################################### */

    if (betaflag==1) {

      vel[j]=A* pow(ztau,-eta);
      deriv_vel=pow(A,2)*pow(r0,2)*eta*pow(ztau,-1+1./(1+eta))*pow(ztau,-1)/(mdot*k_zt);

/* ###################################################################### */
/* Second velocity profile */
/* v(z)=-a*tau */
/* ###################################################################### */

    } else {

      aa=16.*sqrt(3.)/(49.*log(7./3.))*(1/z0)*csi;
      vel[j]= -aa*tau[j];
      deriv_vel=-aa;

    }
  }

 /* ########################################################################## */
/* Seed photon space and energy distrbution*/
/* ########################################################################### */
 
  for (j=0; j<=NTAU; j++) {
 
/* ###########################################################################*/
/* Exponential space distribution of the seed photons */
/* ###########################################################################*/

    cn=1.0344*1E-3;

    bbnorm[j]=cn*exp(-tau[j])/(H);
 
    for (i=0; i<=NQ; i++) {
      q[i]=qmin+i*h; 
      s[i][j]=bbnorm[j]*pow(kte,3)*bb2d(ratio_temp, exp(q[i])); 
    }
  }
  
/* ############################################################################### */
/* Build the initial guess of the spectrum equal to the seed blackbody spectrum */
/* ############################################################################### */

  for (j=0; j<=NTAU; j++) {
    for (i=0; i<=NQ; i++) {

      if (j==NTAU) s[i][j]=0;

      spec_low[i][j]=s[i][j];
      spec_int[i][j]=s[i][j];
      spec_up[i][j]=s[i][j];

    }
  }

/* ############################################################################ */
/* Building the solution of the first equation of 2D progonka method */
/* For any y=j, we calculate the solution of the equation for the Lx-operator */
/* ############################################################################ */

  for (j=0; j<=NTAU; j++) {

    delta_hat=1./(3*H)*deriv_vel;

/* ########################################################################## */
/* Functions of the space operator operator */
/* ########################################################################## */

    W=1./(3*H);
    Z=-vel[j]/H;

    for (i=0; i<= NQ; i++) {

/* ########################################################################## */
/* Functions of the energy operator */
/* ######################################################################### */

      flux_tot[i]=0;
      xnum[i]=exp(q[i]);

      P= (3*kte+ (bulkflag-1)*e_kin_fact*pow(vel[j],2))/(3*kte);
      Q= (3*kte*(exp(q[i])-3+delta_hat)- (bulkflag-1)*e_kin_fact *pow(vel[j],2))/(3*kte);
      R_hat= (exp(q[i])-3*delta_hat)-pow(vel[j],2)*pow(csi,2)/H - 1./ht;

      a[i][j]=P/pow(h,2.);
      b[i][j]=-2*P/pow(h,2)-Q/h+R_hat;
      c[i][j]=P/pow(h,2.)+ Q/h;

      if (j==0) {
  
	d[i][j]=0;
	e[i][j]=-2*W/pow(htau,2)- Z/htau;
	f[i][j]=W/pow(htau,2.)+ Z/htau;
  
      } else if (j==NTAU) {

	d[i][j]=W/pow(htau,2.);
	e[i][j]=-2*W/pow(htau,2)- Z/htau;
	f[i][j]=0;
  
      } else {


	d[i][j]=W/pow(htau,2.);
	e[i][j]=-2*W/pow(htau,2)- Z/htau;
	f[i][j]=W/pow(htau,2.)+ Z/htau;
  
      }

    }
  }


/* ########################################################################### */
/* Boundary conditions on energy (i=0,i=N) for coefficients and L,K                     */
/* ########################################################################### */
  
  for (j=0; j<= NTAU; j++) {
  
    a[0][j]=0; b[0][j]=1; c[0][j]=0; s[0][j]=0;
    a[NQ][j]=0; b[NQ][j]=1; c[NQ][j]=0; s[NQ][j]=0;
   
    L[0][j]=-c[0][j]/(b[0][j]);
    K[0][j]=0;

  }

/* ########################################################################### */
/* Starting the loop over m: at any cicle, we calculate the spectrum at three */
/* layers respect to the time (=m) */
/* ########################################################################### */

  m=0;  
  /* set norm_up here to suppress a compiler warning */
  norm_up = 1.0;
  while (m <= M){  

/* ########################################################################## */
/* Recursively building of the coefficients L, K for any i,j                  */
/* ########################################################################## */

    for (j=0; j<= NTAU; j++) {

/* ########################################################################## */
/* Set to zero at the outer energy-boundary the right-hand side */
/* source term when solving for the energy operator  */
/* ########################################################################## */

      S_hat[0][j]=0;
      S_hat[NQ][j]=0;
      
/* ########################################################################## */
/* Build the source function on the right-hand side of the energy operator */
/* ########################################################################## */

      for (i=1; i< NQ; i++) {
  
	if (j==0) {

	  S_hat[i][j]=-((e[i][j]+1./ht)*spec_low[i][j]+f[i][j]*spec_low[i][j+1]+s[i][j]);
      
	} else if (j==NTAU) {

	  S_hat[i][j]=-(d[i][j]*spec_low[i][j-1]+(e[i][j]+1./ht)*spec_low[i][j]+s[i][j]);
      
	} else {
  
	  S_hat[i][j]=-(d[i][j]*spec_low[i][j-1]+(e[i][j]+1./ht)*spec_low[i][j]+f[i][j]*spec_low[i][j+1]+s[i][j]);
	  
	}

/* ########################################################################## */

	L[i][j]=-c[i][j]/(a[i][j]* L[i-1][j] + b[i][j]);
	K[i][j]=(S_hat[i][j]- a[i][j]*K[i-1][j])/(a[i][j]* L[i-1][j] + b[i][j]); 

      }
   
    }

/* ########################################################################## */
/* Find the solution of the progonka first equation and normalization      */ 
/* ########################################################################## */

    for (j=0; j<=NTAU; j++) {

      spec_int[NQ][j]=0; 
 
      integ_int[j]=0;
      integ_low[j]=0;
      integ_up[j]=0;

/* ########################################################################## */

      for (i=NQ; i>0; i--) {
	spec_int[i-1][j]=L[i-1][j]*spec_int[i][j]+K[i-1][j]; 
	
      }

    }

/* ########################################################################## */
/* Boundary conditions on the optical depth for the coefficients and L,K */
/* of the second equation */
/* ########################################################################## */

    for (i=0; i<= NQ; i++) {

      if (betaflag==1) {
	L_tilde[i][0]=1./(1+htau*(-beta_max*(index+3)+G));
      } else {
	L_tilde[i][0]=1./(1+htau*G);
      }

      K_tilde[i][0]=0;
    }

/* ########################################################################## */
/* Recursively building of the coefficients L_tilde, K_tilde  */
/* for the second progonka over tau */
/* ########################################################################## */
 
    for (i=1; i< NQ; i++) {
      for (j=1; j< NTAU; j++) {
 
  
	L_tilde[i][j]=-f[i][j]/(d[i][j]* L_tilde[i][j-1] + e[i][j]-1./ht);

	if (j==NTAU) {

	  S_tilde[i][j]=(d[i][j]*spec_low[i][j-1]+e[i][j]*spec_low[i][j]-(1./ht)*spec_int[i][j]);

	} else {

	  S_tilde[i][j]=(d[i][j]*spec_low[i][j-1]+e[i][j]*spec_low[i][j]+f[i][j]*spec_low[i][j+1]-(1./ht)*spec_int[i][j]);

	}

	K_tilde[i][j]=(S_tilde[i][j]- d[i][j]*K_tilde[i][j-1])/(d[i][j]* L_tilde[i][j-1] + e[i][j]-1./ht);

      }
   
    }  

/* ########################################################################## */
/* Finding the solution of the progonka first equation and normalization      */ 
/* ########################################################################## */
  
    for (j=0; j<=NTAU; j++) {
      spec_up[0][j]=0;  
      spec_up[NQ][j]=0; 
      integ_up[j]=0;
    }  

    for (i=0; i<=NQ; i++) {
      spec_up[i][NTAU]=0; 
    }


/* ################################################################################ */


    for (j=NTAU; j>0; j--) {
      for (i=1; i< NQ; i++) {

	spec_up[i][j-1]=L_tilde[i][j-1]*spec_up[i][j]+K_tilde[i][j-1]; 
	if (spec_up[i][j-1] > 0 && spec_up[i][j-1] < 1e-30) spec_up[i][j-1]=0;

	integ_low[j]=integ_low[j]+spec_low[i][j-1]*h;
	integ_int[j]=integ_int[j]+spec_int[i][j-1]*h;
	integ_up[j]=integ_up[j]+spec_up[i][j-1]*h;
      }
      
    }
 
/* ############################################################################*/
/* Check the spectral slope */
/* ############################################################################*/

    if (m > 0) previous_index=index;

    index=slope(x, spec_up, qmin, h, ktbb, kte, NQ);


    if (m > 500) {
      double flag=check_conv(previous_index, index);

      if (flag==1 && index > 0.01)  {
	break; 
      }
    }

/* ############################################################################ */

    for (j=0; j<=NTAU; j++) {

      double  vel_tau=vel[j];

      if ((m==maxiter) && ((j==NTAU/2) ||  (j==NTAU-1))){

      }

      for (i=0; i<=NQ; i++) {
 
	if (m < M)  spec_low[i][j]=spec_up[i][j];

  /* ####################################################################################### */
/* Now write the expression for the flux */
/* ####################################################################################### */
  
	if ((j==NTAU) && (i==NQ)) {
	  flux_bulk[i][j]=-(spec_up[i][j]-spec_up[i][j-1])/(3*htau) + vel_tau*(spec_up[i][j]-(1./3)*(spec_up[i][j]-spec_up[i-1][j])/h);
 
	} else if ((j < NTAU) && (i==NQ)) {
	  flux_bulk[i][j]=-(spec_up[i][j+1]-spec_up[i][j])/(3*htau) + vel_tau*(spec_up[i][j]-(1./3)*(spec_up[i][j]-spec_up[i-1][j])/h);
 
	} else if ((j == NTAU) && (i < NQ)) {
	  flux_bulk[i][j]=-(spec_up[i][j]-spec_up[i][j-1])/(3*htau) + vel_tau*(spec_up[i][j]-(1./3)*(spec_up[i+1][j]-spec_up[i][j])/h);

	} else {

	  flux_bulk[i][j]=-(spec_up[i][j+1]-spec_up[i][j])/(3*htau) + vel_tau*(spec_up[i][j]-(1./3)*(spec_up[i+1][j]-spec_up[i][j])/h);

	}

      }

    }

/* ####################################################################################### */
/* Check the normalization of spec_low and spec_up at any optical depth */
/* ####################################################################################### */

    norm_low=0; norm_int=0, norm_up=0, eps=0.05;


    for (j=0; j< NTAU; j++) {
      norm_low=norm_low+integ_low[j]*htau;
      norm_int=norm_int+integ_int[j]*htau;
      norm_up=norm_up+integ_up[j]*htau;
    }


    if ((norm_int/norm_up) > 1-eps &&   (norm_int/norm_up) < 1+eps  && (norm_up  > 0))  { 
      if (m==iter_val+1) {
	count_conv++;
      }
 

    } 

    iter_val=m;
 
    m++;

  }


  for(i=0;i<=NQ;i++){
    flux_tot[i]=flux_bulk[i][NTAU];
  }


/* ################################################################################ */
/* This is the normalization constant which allows for photon conservation number */
/* ############################################################################### */

  phtot[NTAU]=taumax*pow(ktbb,3)*exp(gammln(3)) * zeta(3) /(norm_up);

  for (ie = 0; ie <=Nflux; ie++) {

    ii=ie-1;
    x[ie]=(energy[ie])/kte;

/* ####################################################################################### */

    if (x[ie] < xnum[0]) {

      flux_comp[ie]=flux_comp[0];
      
    } else if (x[ie] > xnum[NQ]) {

      flux_comp[ie]=flux_comp[NQ];

    } else {

      flux_comp[ie]=interp_funct(xnum, flux_tot, NQ, x[ie]);

    }
 
/* ####################################################################################### */

    flux_comp[ie]=flux_comp[ie]/energy[ie]; 

    if (ie > 0) {

      delta_ene=(energy[ie]-energy[ie-1]);

      flux[ii]=0.5*(flux_comp[ie]+flux_comp[ie-1])*delta_ene ;

      flux_bb[ii]=0.5*(bb2d(ratio_temp, x[ie]) +  bb2d(ratio_temp, x[ie-1]))/energy[ie]*delta_ene ;


    }

  }

/* ##################################################################### */

 free_dmatrix(a,0,NQ,0,NTAU);
 free_dmatrix(b,0,NQ,0,NTAU);
 free_dmatrix(c,0,NQ,0,NTAU);
 free_dmatrix(d,0,NQ,0,NTAU);
 free_dmatrix(e,0,NQ,0,NTAU);
 free_dmatrix(f,0,NQ,0,NTAU);
 free_dmatrix(s,0,NQ,0,NTAU);
 free_dmatrix(S_hat,0,NQ,0,NTAU);
 free_dmatrix(L_tilde,0,NQ,0,NTAU);
 free_dmatrix(K_tilde,0,NQ,0,NTAU);
 free_dmatrix(S_tilde,0,NQ,0,NTAU);
 free_dmatrix(L,0,NQ,0,NTAU);
 free_dmatrix(K,0,NQ,0,NTAU);
 free_dmatrix(spec_low,0,NQ,0,NTAU);
 free_dmatrix(spec_int,0,NQ,0,NTAU);
 free_dmatrix(spec_up,0,NQ,0,NTAU);
 free_dmatrix(flux_bulk,0,NQ,0,NTAU);

 free_dvector(xnum, 0,NQ); 
 free_dvector(flux_tot,0,NQ);
 free_dvector(flux_comp,0,Nflux);  
 free_dvector(vel, 0,NTAU);
 free_dvector(tau,0,NTAU);
 free_dvector(phtot,0,NTAU);
 free_dvector(bbnorm,0,NTAU);
 free_dvector(integ_int,0,NTAU);
 free_dvector(integ_low,0,NTAU);
 free_dvector(integ_up,0,NTAU);
 free_dvector(fact,0,NTAU);
 free_dvector(flux_bb,0,Nflux); 
 free_dvector(x,0,Nflux);  
 free_dvector(q, 0,NQ); 

 return(flux);

}


/* #################################################################################### */
/* Blackbody seed energy spectrum */
/* #################################################################################### */

double bb2d(double R, double z) {
  double value;
  value=pow(z,3)/(exp(R*z)-1);
   
   if (value > 1e-30) {
     return(value);
   } else{
     return(0.);
   }

}

/* ############################################################################# */

double check_conv(double previous_index, double index) {

  static int iter=0, count=0, previous_count=0, flag=0;

  int nmax=200;
  double eps=1e-5;

  iter++;

  if (previous_index/index > 1-eps  &&  previous_index/index < 1+eps) {
    count=iter;

    if (count==previous_count+1) { 

      flag++;
    } 

    previous_count=count;
 
  }

/* ################################################################## */
/* If the index value remain within tolerance=eps after */
/* nmax consecutive iterations stop the process */
/* ################################################################## */

  if (flag==nmax) {

    flag=0;
    iter=0;
    count=0;
    previous_count=0;


    return(1);
  } else {

    return(2);

  }

  previous_count=count;
}

/* ###################################################################### */

double slope(double energies[], double **spectrum, int qmin, double h, double ktbb, double kte, int NQ) {

  double alpha=0.2, emin,emax;
  int Nmin, Nmax;

/* ################################################################### */
/* Start from energies > 7 ktbb */
/* ################################################################### */

  emin=7*ktbb;
  emax=30*ktbb;

/* ################################################################### */

  Nmin=(log(emin/kte)-qmin)/h;
  Nmax=(log(emax/kte)-qmin)/h;


  if (spectrum[Nmin][0] < 0 || spectrum[Nmax][0] < 0) {

    alpha=1.0; 
    return(alpha);

  } else {

    alpha=(log(spectrum[Nmax][0])-log(spectrum[Nmin][0]))/(log(emax)-log(emin)); 

  }

 /* if (alpha > 0) alpha=-0.01;  */
  return(-alpha);

}

