/*############################################################################*/
/* Numerical solution of the thermal plus bulk Comptonization equation        */
/* reported in Titarchuk, Mastichiadis & Kylafis (1997)                       */
/* using a finite-difference method and neglecting the second order           */            
/* bulk term (vb/c)^2.                                                        */                 
/* The output spectral shape is formally given by                             */                                                         
/* F(E)=Cn/(A+1) [BB + A G*BB]                                                */                                                    
/* where  G*BB is the convolution between the Green's function                */
/* and the seed BB-like spectrum.                                             */
/*############################################################################*/

#include "numeric.h"

#define N 4000

Real bbody (Real R, Real x, Real bb_par);

Real * xscomptb(Real* energy, int Nflux, Real* parameter, int spectrum, Real*  flux, Real* fluxError, char* init) {

/* ################################################################### */

Real *x=dvector(0,Nflux);
Real *flux_comp=dvector(0,Nflux); 
Real *unnorm_flux=dvector(0,Nflux); 
Real *flux_bb=dvector(0,Nflux); 
Real *bb=dvector(0,Nflux); 

/* ########################################################################## */

 Real ktbb, bb_par, alpha, delta, kte, A, h,  delta_ene, cn, R;
Real fb, gamma, tmin, tmax, phtot, integ=0;
 int ie, ii, n;

static Real a[N+1], b[N+1], c[N+1], s[N+1], p[N+1], g[N+1], r[N+1];
static Real t[N+1], L[N+1], K[N+1], spec[N+1], xnum[N+1];

/* ########################################################################## */
/* Model parameters */
/* ########################################################################## */

ktbb=parameter[0];
bb_par=parameter[1];
alpha=parameter[2];
delta=parameter[3];
kte=parameter[4];
A=pow(10,parameter[5]);

/* ########################################################################## */
/* In TMK97 the Fokker-Planck equation includes also the term */
/* fb=1+(vb/c)^2/(3*theta), but here fb=1 */
/* ########################################################################## */

 fb=1.0;

/* ########################################################################## */
/* Blackbody normalization */
/* ########################################################################## */

cn=8.0525*pow(kte/ktbb,bb_par)*pow(ktbb, bb_par-4);

/* ########################################################################## */

gamma=alpha*(alpha+delta+3)+3*delta;
 
tmin=-10; tmax=10;
h=(tmax-tmin)/(N);
R=kte/ktbb;

/* ########################################################################## */

for (n=0; n<N; n++) {

t[n]=tmin+n*h;
xnum[n]=exp(t[n]);

p[n]=fb;
g[n]=(exp(t[n])-3*fb-delta);
r[n]=(exp(t[n])-gamma+3*delta);

s[n]=-bbody(R, exp(t[n]), bb_par);

if (s[n] > -pow(10,-15)) s[n]=0; 

a[n]=p[n]/pow(h,2);
b[n]=-(2.*p[n]/pow(h,2)+g[n]/h-r[n]); 
c[n]=(p[n]/pow(h,2)+g[n]/h);

}

/* ########################################################################## */
/* Value of the coefficients and of the function at the boundaries            */
/* ########################################################################## */

K[0]=0; L[0]=-c[0]/b[0];

a[0]=0; b[0]=1; c[0]=0; s[0]=0;
a[N]=0; b[N]=1; c[N]=0; s[N]=0;


/* ########################################################################## */

for (n=1; n<N; n++) {
L[n]=-c[n]/(a[n]* L[n-1] + b[n]);
K[n]=(s[n]- a[n]*K[n-1])/(a[n]* L[n-1] + b[n]);

}

/* ########################################################################## */

spec[N]=0;

for (n=N; n>0; n--) {
spec[n-1]=L[n-1]*spec[n]+K[n-1]; 

integ=integ+spec[n-1]*h;
if (spec[n-1] < 0) {
 printf("Warning: with this parameters flux is negative\n\n");
  }
}


/* ########################################################################## */

for (ie = 0; ie <=Nflux; ie++) {

ii=ie-1;
x[ie]=(energy[ie])/kte;

 
flux_comp[ie]=interp_funct(xnum, spec, N, x[ie])/energy[ie];
bb[ie]=bbody(R, x[ie], bb_par)/energy[ie];

if (ie > 0) {

delta_ene=(energy[ie]-energy[ie-1]);
unnorm_flux[ii]=0.5*(flux_comp[ie]+flux_comp[ie-1])*delta_ene ;
flux_bb[ii]=0.5*(bb[ie] +  bb[ie-1])*delta_ene;

}

}

/* ################################################################################*/
/* Normalization constant which allows to keep conserved the photon number         */
/* ############################################################################### */

phtot=pow(ktbb,bb_par)* exp(gammln(bb_par))  * zeta(bb_par) /(integ*pow(kte, bb_par)); 

for (ii = 0; ii < Nflux; ii++) {
flux[ii]= cn/(A+1) * (flux_bb[ii] + A* phtot*unnorm_flux[ii]);
}


 free_dvector(x,0,Nflux);
 free_dvector(flux_comp,0,Nflux); 
 free_dvector(unnorm_flux,0,Nflux); 
 free_dvector(flux_bb,0,Nflux); 
 free_dvector(bb,0,Nflux); 

return(flux);

}


/* ##################################################################*/
/* Seed photon distribution, blackbody or modified blackbody         */
/* ##################################################################*/

Real bbody(Real R, Real z, Real bb_par) {

Real value;

value=pow(z,bb_par)/(exp(R*z)-1);

return(value);

}

