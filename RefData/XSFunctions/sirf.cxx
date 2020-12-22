#include "xsTypes.h"
#include <XSFunctions/functionMap.h>
#include <cmath>
#include <iostream>
#include <fstream>

// written in 2006-2008 by P. Abolmasov,
// see Abolmasov, Karpov \& Kotani (2009)

// just maximum of two values:

Real max(Real x, Real y)
  {
    if(x>y)return(x); else return(y);
  }              

// just minimum of two values:

Real min(Real x, Real y)
  {
    if(x<y)return(x); else return(y);
  }            

// makes a logarythmically-spaced array of n elements, spanning the range from xmin (inclusive) till xmax (non-inclusive)

RealArray logarray(Real xmin, Real xmax, unsigned int n)
{
  Real lxmax=log(xmax), lxmin=log(xmin);
  RealArray y(n);

  for(unsigned int k=0;k<n;++k)
    {
      y[k]=exp((lxmax-lxmin)*(Real)k/(Real)n)*xmin;
      //     std::cerr << y[k] << "\n";
    }
  return y;
}

// auxiliary function f_\alpha (section 2.3 of the paper)

Real funval(Real alpha, Real gamma)
{
  return pow(3.*(1.+alpha),2./(1.+alpha))*(3.+alpha)*pow(2.+alpha,-1.+(gamma-1.)*(2.+alpha));
}

// black body (  x^4/(exp(x)-1))
// caution -- the result is E L_E !!

Real L_bb(Real x)
{
  return (x*x*x*x/(exp(x)-1.));
}

// Y(r) function (formula 21 of the paper)

Real yfun(Real r, Real alpha, Real gamma)
{
  return pow(r,-(2.+alpha)*gamma)*pow(1.-exp(-1./r/(2.+alpha)),1.-(gamma-1.)*(2.+alpha));
}

// the same Y(r) for an array. 

RealArray yfun(RealArray r, Real alpha, Real gamma)
{
  return pow(r,-(2.+alpha)*gamma)*pow(1.-exp(-1./r/(2.+alpha)),1.-(gamma-1.)*(2.+alpha));
}

// kernel used when calculating self-irradiation by the funnel walls
// in principle, it should be iterated

RealArray irradkernel(RealArray  x, Real theta, unsigned int n)
{
      Real cos_theta = cos(M_PI*theta/180.0), sin_theta = sin(M_PI*theta/180.0);

      RealArray a(n);
      for(unsigned int k=0;k<n;++k)a[k]=1.-2.*x[k]*cos_theta*cos_theta+x[k]*x[k];
      RealArray  b(n);
      for(unsigned int k=0;k<n;++k)b[k]=-2.*x[k]*sin_theta*sin_theta;
      RealArray  y(n);
      for(unsigned int k=0;k<n;++k)y[k]=(1.-fabs(1.-x[k])*pow(a[k]-b[k],-1.5)*(a[k]-2.*b[k]))*cos_theta*cos_theta/sin_theta/sin_theta/4.;

      return y;
}

// self-irradiation for a given array of distances r and seed fluxes f0
// theta is the half-opening angle, rin -- inner radius in R_{sph} units. 
// niter is the number of iterations

RealArray selfirrad(RealArray r, RealArray f0, Real theta, Real rin, unsigned int niter)
 { 
   Real dr, i=0., i_tmp, i0=0., i0_tmp;
   // dr = dR/R

   unsigned int nx=r.size();
   RealArray f(nx);
   RealArray f1(nx);

   f1=f=f0;
   unsigned int k, j, q;
   RealArray irrad(nx);
   //      while(r<rmax)
   Real f_tmp;

   for(q=0;q<niter;++q) 
     { // one iteration:
       f1=f;
       for(j=0;j<nx;++j) // j marks the point where flux is calculated
	 {
	   irrad=irradkernel(r / r[j], theta, nx);   
	   i=i0=0.;
	   for(k=0;k<nx-1;++k) // k for the source annulus 
	     { 
	       dr = r[k+1]-r[k];
	       if(q==0)f_tmp=f0[k]; else f_tmp=f[k]-f0[k];
	       f_tmp=f_tmp/f0[j];
	       if(k!=0)
		 {
		   if(k%2==0)
		     {
		       i_tmp = f_tmp * irrad[k]*dr*2./3.; 
		       i0_tmp = irrad[k]*dr*2./3.; 
		     }
		   else
		     {
		       i_tmp = f_tmp * irrad[k]*dr*4./3.;
		       i0_tmp = irrad[k]*dr*4./3.;
		     }
		 } 
	       else
		 {
		   i_tmp = f_tmp * irrad[k]*dr/3.; 
		   i0_tmp = irrad[k]*dr/3.; 
		 }
	       //	       if(i0_tmp>1.e-4*(Real)nx)
	       i += i_tmp;
	       i0 += i0_tmp;
	     }
	   f1[j] *= (i / r[j] / r[j] + 1.);
	 }
       f0=f;
       f=f1;
     }

   return f;
 }

// the function calculates the size of the visible fraction of the annulus:

Real obscur(Real x, Real k)
  {
    Real y=1.;
    Real sina=0.;

    if(k > 0.)sina=(k*k+1.)/(k*2.)+(1.-k*k)/(k*2.*x);

    if(k >= 1.)
     {
      if(sina <= -1.)
       {
        y=0.;
       }
      else
       {
        y=.5+asin(sina)/M_PI;
        if(y<0)y=0.;
       }       
     }
    return y; 
  }

// E L_E for funnel walls:

Real L_w(RealArray r, RealArray tar, Real x,  Real theta, Real tgrat, RealArray v)
{
  unsigned int nx=r.size();
  Real rmax=r[nx-1];
  Real dr, rc, i=0., i_tmp, tau;
  Real cos_theta = cos(M_PI*theta/180.0);
  Real sin_theta = sin(M_PI*theta/180.0);
  Real delta=1.; // -- Doppler factor

  for(unsigned int k=0;k<nx-1;++k)
    {
      dr = r[k+1]-r[k];
      rc = 0.5*(r[k+1]+r[k]);
      tau = (tar[k+1]+tar[k])/2.;
      // let's not bother we relativistic effects at low flow velocities (\beta < 0.01):
      if(v[k]<0.01)delta=1.; else delta=1./((1.-v[k]*cos_theta)/sqrt(1.-v[k]*v[k]));
      i_tmp = cos_theta*sin_theta*rc*dr/(exp(x/tau/delta) - 1); 
      if(v[k]>0.01)i_tmp*=delta*delta*delta;
      i_tmp *= obscur(r[k]/rmax, tgrat); 
      if(k!=0) // numerical integration (Simpson) :
	{
	  if(k%2==0)
	    {
	      i_tmp *= 2./3.; 
	    }
	  else
	    {
	      i_tmp *= 4./3.;
	    }
	} 
      i += i_tmp;
    }

  //  i*=2.*M_PI;

  return x*x*x*x*i;
}

// E L_E for the bottom
// sometimes it is obscured, sometimes visible. 
// We consider the bottom point-like, so it is visible if and only if i < theta_f
// cos_incl is not used, we multiply by this Lambert factor in the main function (see below)

Real L_bottom(Real x, Real theta, Real incl)
{
    Real cos_theta = cos(M_PI*theta/180.0);
    Real sin_theta = sin(M_PI*theta/180.0);
    Real cos_incl = cos(M_PI*incl/180.0);

    Real y;

    if(fabs(cos_incl) >= fabs(cos_theta))
     { 
       y=L_bb(x); //x*x*x*(1. - cos_theta)/(exp(x) - 1.);
     }
    else
     {
      y=0.;
     } 
    y*=sin_theta*sin_theta/2.;

    return y;
}

// photosphere projection surface area

Real surph(Real theta, Real incl)
{
  Real cos_theta = cos(M_PI*theta/180.0), sin_theta = sin(M_PI*theta/180.0);
  Real cos_incl = cos(M_PI*incl/180.0), sin_incl = fabs(sin(M_PI*incl/180.0));

  if(sin_incl <= cos_theta)
    {
      return((-sin_theta*cos_incl+1.)/2.); 
    }
  else
    {
      return(M_PI*cos_theta*cos_theta);
//(1.-2.*theta/180.)+sin(2.*M_PI*theta/180.0));
    }

  //  std::cerr << "surph: error!"; 

  return 0.;
}

// calculates temperature radial profile with irradiation by iterative method. 

RealArray tempcalc(RealArray r, Real theta, Real val, Real gamma, int air, Real rin, Real* tmin)
{
  // air = niter -- number of iterations

    Real cos_theta = cos(M_PI*theta/180.0);
    Real sin_theta = sin(M_PI*theta/180.0);

    unsigned int nx=r.size();

    RealArray x(nx);
    RealArray f0(nx);
    RealArray tar(nx);
    RealArray dbot(nx);

    x=r/rin;

    f0=yfun(r, val, gamma) * pow(x,1.+val); // undisturbed (unirradiated) normal flux profile
    *tmin=pow(f0[nx-1]/f0[0],0.25);
    dbot=(x-1.)*sin_theta*sin_theta*sin_theta*cos_theta*cos_theta / 2. / (1. + x*x - 2.*x*cos_theta)/ (1. + x*x - 2.*x*cos_theta);
    if(air>0)f0+=dbot*rin*rin;
    if(air>0)f0=selfirrad(x,f0,theta,rin,air);

    tar=pow(f0/f0[0],0.25);

    return(tar);
}

// the main function

void sirf(const RealArray& energy, const RealArray& parameter, 
	  int spectrum, RealArray& flux, RealArray& fluxError, 
	  const string& init)
{
  unsigned int n=energy.size()-1;
  flux.resize(n);

  Real tin=parameter[0], rin=parameter[1], rout=parameter[2], theta=parameter[3], incl=parameter[4], val=parameter[5],  gamma=parameter[6], mdot=parameter[7];

  // mdot affects the velocity profile. At lower mdot v is higher, v \propto mdot^{-1/2}

  Real tmin0[1];
  Real *tmin;
  tmin=tmin0;

  int air=(int)parameter[8];

  unsigned int nx=256; // number of intergation points along the flow

  RealArray tar(nx);
  RealArray r(nx);

  r=logarray(rin,rout,nx);

  RealArray v(nx);
  Real cos_theta=cos(theta/180.*M_PI);

  v=sqrt(cos_theta/mdot)/6.*pow(r,val); // the only place where mdot is really needed
  if (v.max() > 0.99)
    {
      for(unsigned int j=0;j<v.size();++j)v[j]=min(v[j],0.99);
    }

  tar=tempcalc(r, theta, val,gamma, air,rin,tmin);

  Real cos_incl=cos(incl/180.*M_PI);
  Real tgrat=tan(M_PI*incl/180.0)/tan(M_PI*theta/180.0); // tg(i)/tg(\theta)

  Real ec, ed, sur=surph(theta,incl);
  Real delta_bot=1.;
  if(v[0]>0.01)delta_bot=sqrt(1.-v[0]*v[0])/(1.-v[0]*cos_incl); // Doppler boost factor for the bottom (significant only if the motions are relativistic)

  for(unsigned int jl=0;jl<flux.size();++jl)
    {
       //       ++enar2;
      ec=0.5*(energy[jl]+energy[jl+1]);
      ed=fabs(-energy[jl]+energy[jl+1]);
      flux[jl]=( ( L_w(r, tar, ec/tin, theta, tgrat,v) + // L_.. calculates EF_E, but dE n_E (photon number) is needed!
		delta_bot  * delta_bot* L_bottom(ec/tin/delta_bot,theta,incl)*rin*rin ) * cos_incl + 
		L_bb(ec/tin / *tmin) * sur * *tmin * *tmin * *tmin * *tmin * rout * rout) * ed/ec/ec;
      //  fluxError[jl]+=sqrt(flux[jl]);
    }

}
