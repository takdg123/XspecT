#ifndef g77Fortran
#define g77Fortran 1
#endif

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "bwcycl.h"

/*Max terms in series*/
#define MAX_TERM 70
/*Relative precision in series and integrals*/
#define TOL 1e-5

#define FILEOUTPUT 0


#define locEPS (1000.0*GSL_DBL_EPSILON)


#define askInfo(a,b) _askInfo(a,b,sizeof(b))


/*Useful constants*/
static const double electron_mass_kev = 510.998918;
static const double kelvin2keV = GSL_CONST_CGSM_BOLTZMANN/GSL_CONST_CGSM_ELECTRON_VOLT/1e3;
static const double kev2erg = GSL_CONST_CGSM_ELECTRON_VOLT*1e3; /*1.6e-9*/

/*Declaration of global variables which depend only on model parameters*/
static double alpha, w, kappa, sigma_par, sigma_med, tau_th, T_th, tau_max, chi_abs;
static double Xn_vec[MAX_TERM], An_vec[MAX_TERM], cyc_ser_vec[MAX_TERM], cyc_ser_vec_sp[MAX_TERM], Green_ser_vec[MAX_TERM], FF_ser_vec[MAX_TERM], lambda_vec[MAX_TERM], mu_vec[MAX_TERM], gn_vec[MAX_TERM];
static double norm_Green, norm_FF, norm_cyc, norm_BB;


/*In the implementation the parameter  xi is called zeta */

/*Functions from A.S. ask library for keybord I/O*/
void  _askInfo(char *Query, char *Answ,int MaxLen)
{

	char MAXCHARS = 121;

    char Answer[MAXCHARS];



    char           *RetChar;
    int AnsLen;

        /*-------------------------------------------
         * printout the query and the default answer
         *-------------------------------------------*/


    printf("\n%s ==> ", Query);

    if (Answ[0] != '\0')
        printf("[%s] ", Answ);  /* default value */




    
    if (fgets(Answer, MAXCHARS, stdin)==NULL)
      printf("\n askInfo: ERROR!");


    if ((RetChar = strchr(Answer, '\n')) != NULL)
      *RetChar = '\0';


    /*-----------------------------------------------------
     * check whether the length of the answer was too long
     *-----------------------------------------------------*/

    if ((AnsLen=strlen(Answer))==(MAXCHARS-1))
      printf("\n WARNING: Input string too long was truncated to %d chars!\n",MAXCHARS-1);

    /*-----------------------------------TOL
     *  remove the ending newline char
     *----------------------------------*/




    if (Answer[0] == '\0')
        return;                 /* return the default */


    /*------------------------------------------
     * New january 2001: check for max lenght
     *------------------------------------------*/

    if (AnsLen>=MaxLen)
      {
	printf("\n WARNING: Input string too long was truncated to %d chars!\n",MaxLen-1);
	memcpy(Answ,Answer,MaxLen-1);
	Answ[MaxLen-1]='\0';
      }
    else
      strcpy(Answ, Answer);

//printf("\n AnsLen=%d; MaxLen=%d ",AnsLen,MaxLen);

    return;
}

void  waitkey(void)
{
        askInfo(" hit the enter key to continue ...","");
}


/*Modification of the GSL 2F1 function to compute the 2F2, not complete,but it works */

int
hyperg_2F2_series_e(const double a, const double b, const double c, const double d,
                  const double x, 
                  gsl_sf_result * result
                  )
{
  double sum_pos = 1.0;
  double sum_neg = 0.0;
  double del_pos = 1.0;
  double del_neg = 0.0;
  double del = 1.0;
  double k = 0.0;
  int i = 0;

  if(fabs(c) < GSL_DBL_EPSILON) {
    result->val = 0.0; /* FIXME: ?? */
    result->err = GSL_NAN;
    GSL_ERROR ("error", GSL_EDOM);
  }

  do {
    if(++i > 30000) {
      result->val  = sum_pos - sum_neg;
      result->err  = del_pos + del_neg;
      result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
      result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k)+1.0) * fabs(result->val);
      GSL_ERROR ("error", GSL_EMAXITER);
    }
    del *= (a+k)*(b+k) * x / ((c+k) * (d+k) * (k+1.0));  /* Gauss series */

    if(del > 0.0) {
      del_pos  =  del;
      sum_pos +=  del;
    }
    else if(del == 0.0) {
      /* Exact termination (a or b was a negative integer).
       */
      del_pos = 0.0;
      del_neg = 0.0;
      break;
    }
    else {
      del_neg  = -del;
      sum_neg -=  del;
    }

    k += 1.0;
  } while(fabs((del_pos + del_neg)/(sum_pos-sum_neg)) > GSL_DBL_EPSILON);

  result->val  = sum_pos - sum_neg;
  result->err  = del_pos + del_neg;
  result->err += 2.0 * GSL_DBL_EPSILON * (sum_pos + sum_neg);
  result->err += 2.0 * GSL_DBL_EPSILON * (2.0*sqrt(k) + 1.0) * fabs(result->val);

  return GSL_SUCCESS;
}

int
gsl_sf_hyperg_2F2_e(double a, double b, const double c, const double d,
                       const double x,
                       gsl_sf_result * result)
{

   //printf("2F2 %f %f %f %f %f\n", a, b, c, d, x);
  //const double e = d + c - a - b;  //?????
  const double rinta = floor(a + 0.5);
  const double rintb = floor(b + 0.5);
  const double rintc = floor(c + 0.5);
  const double rintd = floor(d + 0.5);
  const int a_neg_integer = ( a < 0.0  &&  fabs(a - rinta) < locEPS );
  const int b_neg_integer = ( b < 0.0  &&  fabs(b - rintb) < locEPS );
  const int c_neg_integer = ( c < 0.0  &&  fabs(c - rintc) < locEPS );
  const int d_neg_integer = ( d < 0.0  &&  fabs(d - rintd) < locEPS );

  result->val = 0.0;
  result->err = 0.0;

//   if(x < -1.0 || 1.0 <= x) {
//     DOMAIN_ERROR(result);
//   }

  if(c_neg_integer) {
    if(! (a_neg_integer && a > c + 0.1)) DOMAIN_ERROR(result);
    if(! (b_neg_integer && b > c + 0.1)) DOMAIN_ERROR(result);
  }

  if(d_neg_integer) {
    if(! (a_neg_integer && a > d + 0.1)) DOMAIN_ERROR(result);
    if(! (b_neg_integer && b > d + 0.1)) DOMAIN_ERROR(result);
  }

  if(fabs(c-a) < locEPS){
    return gsl_sf_hyperg_1F1_e(b,d,x,result);
  }

  if( fabs(c-b) < locEPS) {
    return gsl_sf_hyperg_1F1_e(a,d,x,result);
  }
  if(fabs(d-a) < locEPS){
    return gsl_sf_hyperg_1F1_e(b,c,x,result);
  }

  if( fabs(d-b) < locEPS) {
    return gsl_sf_hyperg_1F1_e(a,c,x,result);
  }



  //if(a >= 0.0 && b >= 0.0 && c >=0.0 && x >= 0.0 && x < 0.995) {
    /* Series has all positive definite
     * terms and x is not close to 1.
     */
    return hyperg_2F2_series_e(a, b, c, d, x, result);
}


/*-*-*-*-*-*-*-*-*-* GSL Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

double gsl_sf_hyperg_2F2(double a, double b, double c, double d, double x)
{
  EVAL_RESULT(gsl_sf_hyperg_2F2_e(a, b, c, d, x, &result));
}

/*end of the hack from GSL*/

/*Beginning the model*/

int alphaf_e(double zeta, double R_star, double M_star, gsl_sf_result * result )
{
    
    double c=GSL_CONST_CGSM_SPEED_OF_LIGHT;
    result->err = 0.0;
    result->val = 32.0 * sqrt(3.0) * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * M_star / 49.0 / log(7./3.) / R_star / c/c * zeta;
    return GSL_SUCCESS;
}

double alphaf(double zeta, double R_star, double M_star)
{
    EVAL_RESULT(alphaf_e( zeta,  R_star,  M_star, &result ));
}

double wf(double zeta)
{

	return sqrt(9.0+12.0*zeta*zeta);
}

double kappaf(double delta)
{

	return 0.5*(delta+4.0);
}

double lambdaf(double n)
{
	return (4.0*n*w + w +3.0)/2.0;
}


int muf_e(double delta, double zeta, double n, gsl_sf_result *result)
{

    double tmp=(3.0-delta)*(3.0-delta) + 4.0*delta*lambda_vec[(int)n];
    if(tmp <0)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("muf_e :: sqrt argument should be positive", GSL_EDOM);
    }
    result -> val=0.5 * sqrt( tmp );
	return GSL_SUCCESS;
}

int Xn_e(double zeta, double n, gsl_sf_result * result)
{
    if(n<0)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("Xn_e :: factorial argument should be positive", GSL_EDOM);
    }
    gsl_sf_result result_gamma;
    int status=GSL_SUCCESS;
    status=gsl_sf_gamma_e(n+0.5, &result_gamma);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val = GSL_NAN;
       GSL_ERROR ("Gamma function", status);
    }
    
    double t1 = result_gamma.val;
	double t2 = pow((3.0-w),n-1.0);
	double t3 = (3.0-w-4.0*n*w);
	double t4 = gsl_sf_fact((unsigned int)n);
	double t5 = pow(alpha,1.5);
	double t6 = pow ( (3.0+w), n+1.5);
    result->val = 2. * t1*t2*t3/t4/t5/t6;
	return GSL_SUCCESS;

}

int An_factor_e(double zeta, double n, double m, gsl_sf_result *result)
{
    
    if(m<0 || (n-m) < 0)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("Xn_e :: factorial argument should be positive", GSL_EDOM);
    }
    
    gsl_sf_result result_gamma;
    int status=GSL_SUCCESS;
    status=gsl_sf_gamma_e(n+0.5, &result_gamma);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("Gamma function", status);
    }
    
    double t1 = result_gamma.val;
    
	double t2 = pow((2.0*w),m);
	double t3 = pow( (3.0-w), -m);
    status=gsl_sf_gamma_e(m+0.5, &result_gamma);
    if(status!= GSL_SUCCESS && status != GSL_EOVRFLW)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("Gamma function", status);
    }
	double t4 = result_gamma.val;
    
	double t5 = 2.0*gsl_sf_fact((int)m);
	double t6 = gsl_sf_fact((int)(n-m));
	
    double t7_1= alpha*(w-3.)/4.;
    status=gsl_sf_gamma_inc_e ( m , t7_1 * tau_th * tau_th, &result_gamma);
    if(status)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("Inc Gamma function", status);
    }
    double t7 = result_gamma.val;
    status=gsl_sf_gamma_inc_e ( m , t7_1 * tau_max * tau_max, &result_gamma);
    if(status)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("Inc Gamma function", status);
    }
	double t8 = result_gamma.val;
    
    if (gsl_isinf(t4))
    {
        result->val=0;
    }
    else
    {
        result->val= t1*t2*t3/t4/t5/t6*(t7-t8);
    }
            
    
    
    return GSL_SUCCESS;
    
}

int An_e(double zeta,  double n, gsl_sf_result *result )
{
	double m=0;
	double An_m = 0;
	gsl_sf_result result_An;
    int status=0;
	for(m=0;m<=n;m++)
	{
        status=An_factor_e(zeta, n, m, &result_An);
        if(status)
        {
            result->val = GSL_NAN;
            GSL_ERROR ("An factor error", status);
        }
        An_m += result_An.val;
	}
    result->val=An_m;
    
	return GSL_SUCCESS;
}


int An_sp_e(double zeta,  double n, double DT, gsl_sf_result *result )
{
    double tmp=log(7./3.) / sqrt(3) / alpha / zeta;
    if(tmp <0)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("AN_sp_e :: sqrt argument should be positive", GSL_EDOM);
        
    }
	double tau_sp = sqrt( tmp);
	//sonic point
    gsl_sf_result result_gn;
    
    int status= gn_e(tau_sp,n, &result_gn);
    if(status)
    {
        result->val = GSL_NAN;
        GSL_ERROR ("Error in gn_e", status);
        
    }
    double tmp1=result_gn.val;
    result->val = exp(1.5*alpha*tau_sp*tau_sp) * tmp1 * DT ;
    
    return GSL_SUCCESS;
}

double An_sp(double zeta,  double n, double DT)
{
    
    EVAL_RESULT( An_sp_e(zeta, n, DT, &result));
}


double	sigma_parf(double r0, double Mdot, double zeta)
{
	double t0= ( M_PI * r0 * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / Mdot / zeta) / sqrt(GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
	return t0*t0;
}


double T_thf(double Mdot, double r0)
{
	return  2.32e3 * pow(Mdot, 2./5.) * pow(r0, -2./3.);
}

double tau_thf(double Mdot, double r0, double R_star, double M_star, double zeta, double *T_th)
{
	*T_th = T_thf(Mdot, r0);
    /*Possible errors are handled at parameter sanity check level*/
	return 2.64e28 * Mdot * (R_star) / (M_star) / pow(r0,1.5) / pow(*T_th, 7./4.) / zeta;
}

double z_thf(double Mdot, double r0, double R_star, double M_star, double zeta, double *T_th)
{
	*T_th = T_thf(Mdot, r0);
	return 5.44e15 * Mdot * (R_star) / (M_star) / r0 / pow(*T_th, 7./2.) / zeta / sigma_par;
}

double tau_maxf(double Mdot, double r0, double R_star, double M_star, double zeta)
{
	
	double r_sig_p_o = ( M_PI * r0 * GSL_CONST_CGSM_MASS_PROTON * GSL_CONST_CGSM_SPEED_OF_LIGHT / Mdot / zeta) / (GSL_CONST_CGSM_THOMSON_CROSS_SECTION);
	/*square root of the ratio between optical depths parallel / orthogonal*/
		

	double C1 = 4. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT *  r0 * zeta / alpha *(M_star/ gsl_pow_2(GSL_CONST_CGSM_SPEED_OF_LIGHT*R_star)) / r_sig_p_o;
	
	double z_max = R_star / 2.0;
 
	if(fabs(C1)<1e-4)
		z_max *= C1 / 2.;
	else
		z_max *= ( sqrt(1.0+C1) -1.0 );
	
	return sqrt(r_sig_p_o) * sqrt( 2.0 * z_max / alpha / zeta / r0);

}

int ser_cyc_factor_e(double e, double n, double *params, gsl_sf_result *result)
{

	double e_cyc = params[4];
	double Te = params[6];

	double mu = mu_vec[(int)n];

	double chi_min = e;
	double chi_max = e;
	double t5, t6;
	if(e_cyc < chi_min)
		chi_min = e_cyc;

	if(e_cyc > chi_max)
		chi_max = e_cyc;
	
	chi_min /= Te;
	chi_max /= Te;

	double t2 = 0.5+mu-kappa;

	double t3 = 1+2*mu;
	

    gsl_sf_result result_sf;
    
    int status=gsl_sf_hyperg_1F1_e(t2,t3,chi_min, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_1F1_e", status);
    }
    double t_5_1=result_sf.val;
    
	t5 = exp(-0.5*chi_min) * pow(chi_min,mu+0.5) * t_5_1;
	/*whittaker MU (kappa, mu, energy, e_cyc)*/
    status=gsl_sf_hyperg_U_e(t2,t3,chi_max, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_1F1_e", status);
    }
    double t_6_1=result_sf.val;
	t6 = exp(-0.5*chi_max) * pow(chi_max,mu+0.5) * t_6_1;
	/*whittaker W (kappa, mu, energy, e_cyc)*/
    result->val=cyc_ser_vec[(int)n]*t5*t6;
    
    return GSL_SUCCESS;
    
}


double ser_cyc_factor_sp(double e, double n, double *params, double DT)
{

	double e_cyc = params[4];
	double Te = params[6];

	double mu = mu_vec[(int)n];

	double chi_min = e;
	double chi_max = e;
	double t5, t6;
	if(e_cyc < chi_min)
		chi_min = e_cyc;

	if(e_cyc > chi_max)
		chi_max = e_cyc;
	
	chi_min /= Te;
	chi_max /= Te;

	double t2 = 0.5+mu-kappa;

	double t3 = 1+2*mu;
	


	t5 = exp(-0.5*chi_min) * pow(chi_min,mu+0.5) * gsl_sf_hyperg_1F1(t2,t3,chi_min);
	/*whittaker MU (kappa, mu, energy, e_cyc)*/
	t6 = exp(-0.5*chi_max) * pow(chi_max,mu+0.5) * gsl_sf_hyperg_U(t2,t3,chi_max);
	/*whittaker W (kappa, mu, energy, e_cyc)*/

	return cyc_ser_vec_sp[(int)n]*t5*t6;
}

int ser_cyc_e(double e, double *params, gsl_sf_result *result)
{
	double n=0;
    gsl_sf_result result_sf;
    int status=ser_cyc_factor_e(e,n,params,&result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in ser_cyc_factor_e", status);
    }
    double curr = result_sf.val;
	double new=curr;

	for(n=1;n<MAX_TERM;n++)
	{
        status=ser_cyc_factor_e(e,n,params,&result_sf);
        if(status)
        {
            result->val=GSL_NAN;
            GSL_ERROR ("Error in ser_cyc_factor_e", status);
        }
        
		new += result_sf.val;
		if( fabs(new/curr)-1 <= TOL  || new == curr)
		{
            result->err=fabs(new-curr);
			break;
		}
		else
		{
			curr=new;
		}
	}


	if(n>=MAX_TERM-1)
	{
		fprintf(stderr, "Warning Cyclotron series reached the maximum number of terms\n");
	}
    result->val=new;
    
	return GSL_SUCCESS;
}


double ser_cyc_sp(double e, double *params, double DT)
{
	double n=0;
	double curr = ser_cyc_factor_sp(e,n,params,DT);
	double new=curr;

	for(n=1;n<MAX_TERM;n++)
	{
		new += ser_cyc_factor_sp(e,n,params,DT);
		if( fabs(new/curr)-1 <= TOL  || new == curr)
		{
			break;
		}
		else
		{
			curr=new;
		}
	}


	if(n>=MAX_TERM-1)
	{
		fprintf(stderr, "Warning Cyclotron series reached the maximum number of terms\n");
	}
	
	return new;
}


double Hf(double chi)
{
	double t0 = 0.41;
	
	if(chi<7.5 && chi >=0 )
		t0=0.15*sqrt(chi);

	return t0;
}

int cyc_e(double e, double *params, gsl_sf_result *result)
{
	/*cyclotron emisssion*/
	double e_cyc = params[4];
	double Te = params[6];
	double t1 =(pow(e*kev2erg,kappa-2) / (pow(e_cyc*kev2erg,kappa+1.5)) / exp(0.5*(e+e_cyc)/Te));
    gsl_sf_result result_sf;
    int status =  ser_cyc_e(e, params, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in ser_cyc_e", status);
    }
    
    double t2 = result_sf.val;
    
    result->val=t2 * (norm_cyc * t1);
    
	return GSL_SUCCESS;
}

double cyc_sp(double e, double *params, double DT)
{
	/*cyclotron emisssion*/
	double e_cyc = params[4];
	double Te = params[6];
	double t1 =(pow(e*kev2erg,kappa-2) / (pow(e_cyc*kev2erg,kappa+1.5)) / exp(0.5*(e+e_cyc)/Te));
	double t2 =  ser_cyc_sp(e, params, DT);
	return t2 * (norm_cyc * t1);
}



double tauf(double z, double zeta, double r0)
{
	//zeta is the greek zeta
	//z is the altitude
	double t1 = sigma_par / GSL_CONST_CGSM_THOMSON_CROSS_SECTION;
	return pow( t1 , 0.25) * sqrt(2.0*z/alpha/zeta/r0);

}


double zf(double tau, double zeta, double r0)
{
	//zeta is the greek zeta
	//z is the altitude
	double t1=GSL_CONST_CGSM_THOMSON_CROSS_SECTION/sigma_par;
	return 0.5*tau*tau * sqrt(t1)*alpha*zeta*r0;
}


double rhof(double tau, double zeta, double r0, double mdot)
{
	//zeta is the greek zeta
	//z is the altitude

	return mdot/r0/r0/GSL_CONST_CGSM_SPEED_OF_LIGHT/M_PI/alpha/tau;
}




int gn_e(double tau, int n, gsl_sf_result *result)
{

	double tau2=tau*tau;
	double t1 = exp(-alpha * (3.0+w) * tau2/4.0);
    gsl_sf_result result_sf;
    int status =gsl_sf_laguerre_n_e(n, -0.5, alpha*w*tau2/2.0, &result_sf);
    /*
     if(status == GSL_EUNDRFLW)
    {
        printf("Underflow at %g %d\n", tau, n);
    }
     */
    if(status != GSL_EUNDRFLW && status != GSL_SUCCESS)
    {
        printf("Values %g +/- %g\n", result_sf.val, result_sf.err);
        result->val=0;
        GSL_ERROR ("Error in gsl_sf_laguerre_n_e", status);
    }
	double t2 = result_sf.val;

    result->val=t1*t2;
	return GSL_SUCCESS;
}


int ser_Green_factor_e(double tau0, double e0, double e, double n, double *params, gsl_sf_result *result)
{
	double Te = params[6];

	double mu = mu_vec[(int)n];
	double t1 = 0.5+mu-kappa;
	double t2 = 1.0+2.0*mu;


	double emin=e0;
	double emax=e0;

	double t4,t5,t6;

	mu+=0.5;

	t4 = gn_vec[(int)n];	

	if(e<e0)
		emin=e;

	if(e>e0)
		emax=e;

	emin/=Te;
	emax/=Te;
	/*mu is  mu+0.5 !!!*/
    gsl_sf_result result_sf;
    int status =gsl_sf_hyperg_1F1_e(t1,t2,emin,&result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_1F1", status);
    }
	t5 = exp(-0.5*emin) * pow(emin,mu) * result_sf.val;
    
    status=gsl_sf_hyperg_U_e(t1,t2,emax,&result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_U", status);
    }
	t6 = exp(-0.5*emax) * pow(emax,mu) * result_sf.val;

    result->val= Green_ser_vec[(int)n]*t4*t5*t6;
    
	return GSL_SUCCESS;
}


int ser_Green_e(double tau0, double e0, double e, double *params, gsl_sf_result *result)
{
	double n=0;
    gsl_sf_result result_sf;
    int status = ser_Green_factor_e(tau0,e0,e,n,params,&result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in ser_Green_factor_e", status);
    }
    double curr = result_sf.val;
	double new=curr;

	for(n=1;n<MAX_TERM;n++)
	{
        status = ser_Green_factor_e(tau0,e0,e,n,params,&result_sf);
        if(status)
        {
            result->val=GSL_NAN;
            GSL_ERROR ("Error in ser_Green_factor_e", status);
        }
		new += result_sf.val;
        //printf("%.0f %g\n", n, fabs(new-curr)/curr);
        
        result->err=fabs(new-curr);
        
		if( fabs(new-curr)/curr <= TOL || curr == new)
		{
            /*printf("Stop Green function at %f\n", n);*/
			break;
		}
		else
		{
			curr=new;
		}
	}
	if(n>=MAX_TERM-1)
	{
		fprintf(stderr, "Warning Green series reached the maximum number of terms\n");
	}
	
    result->val=new;
	return GSL_SUCCESS;
}

int Green_e(double tau0, double e0, double e, double *params, gsl_sf_result *result)
{
    gsl_sf_result result_sf;
    int status =ser_Green_e(tau0, e0, e, params,&result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in ser_Green_e", status);
    }
	/*column intergrated Green function*/
    double res = result_sf.val;
	double Te = params[6];
	if(res>0)
		res*=norm_Green * pow(e0, -kappa) * pow(e,kappa-2.0) * exp(0.5*(e0-e)/Te) * exp(1.5*alpha*tau0*tau0);;
    result->val=res;
	return GSL_SUCCESS;

}

double bb_integrand(double e0, void *params)
{
	double *int_params = (double *)params;
	double bb = e0*e0/gsl_sf_expm1(e0/(T_th*kelvin2keV));
    gsl_sf_result result_sf;
    int status = Green_e(tau_th, e0, int_params[8], params, &result_sf);
    if(status)
    {
        gsl_error ("Error in Green_e", __FILE__, __LINE__, status);
        bb=0;
    }
    else
    {
        bb *= result_sf.val;
    }
    
	return bb;
	
}


double BB(double e, double *params)
{
	/*Black body reprocessed emission*/
	double r0 = params[7];
	double D  = params[8];

	gsl_function F;
	double result, error, T_th_kev ;
	int err_code;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	
	/*Use it to store*/
	params[8] = e;

	F.function = &bb_integrand;
	F.params   = params;
	
	T_th_kev = T_th*kelvin2keV;

	err_code = gsl_integration_qag (&F, T_th_kev/20., T_th_kev*20.0, 1e-6, TOL*10, 1000,GSL_INTEG_GAUSS41, w, &result, &error);


	switch (err_code)
	{
		case GSL_EMAXITER : 
		{
			fprintf(stderr, "BB :: the maximum number of subdivisions was exceeded\n");
		}
		break;
		case GSL_EROUND :
		{
			fprintf(stderr, "BB :: cannot reach tolerance because of roundoff error, or roundoff  error was detected in the extrapolation table\n");
		} break;
		case GSL_ESING : 
		{
			fprintf(stderr, "BB :: a non-integrable singularity or other bad integrand behavior was found in the integration interval.\n");
		} break;
		case GSL_EDIVERGE : 
		{
			fprintf(stderr, "BB :: the integral is divergent, or too slowly convergent to be integrated numerically.\n");
		} break;
		case 0 :
		{
			;/*printf("BB :: Integration succesful\n");*/
        
		} break;
		default:
		{
			fprintf(stderr, "BB :: Unknown error in the integration.\n");
			fprintf(stderr, "Error code is %d '%s'\n", err_code ,gsl_strerror (err_code));

		}
	}


	gsl_integration_workspace_free(w);
	/*printf("%g\t%g\n", norm, result);*/
	params[8] = D;
	/* 	printf("BB = %g +/- %g\n", result, error);*/
	return norm_BB*r0*r0*result;
}


/*Free-Free*/


double chi_abs_int(double tau, void *int_params)
{
	/*Te is in keV transform in Kelvin !*/
	double *params = (double *)int_params;

	double Te = params[6];

	double norm = 6.08e12*pow(Te/kelvin2keV, -7./4.);
	double rho = params[5] / (params[7] * params[7] * GSL_CONST_CGSM_SPEED_OF_LIGHT );
	rho /= ( M_PI * alpha * tau );

	//params[7] is r0
	//params[5] is Mdot
	return norm * sqrt(rho);

}


double chi_absf(double *params)
{

	double t1 = chi_abs_int(tau_th, params);
	double t2 = chi_abs_int(tau_max, params);
	

	
	if(t1 >0 && t2>0)
		return sqrt(t1*t2);//return 2. * t1 * t2 / (t1 + t2);
	else
	{
		fprintf(stderr, "Error in computation of chi_abs\n");
        return 0;
	}
	
}

double chi_absf_int(double *params)
{


	double result, error;
	int err_code;

	gsl_function F;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	F.function = &chi_abs_int;
	F.params = params;
	

	err_code = gsl_integration_qag (&F,0,tau_max, 1e-6, 0.1*TOL, 1000, GSL_INTEG_GAUSS41, w, &result, &error);

	switch (err_code)
	{
		case GSL_EMAXITER : 
		{
			fprintf(stderr, "chi_absf :: the maximum number of subdivisions was exceeded\n");
		}
		break;
		case GSL_EROUND : 
		{
			fprintf(stderr, "chi_absf :: cannot reach tolerance because of roundoff error, or roundoff  error was detected in the extrapolation table\n");
		} break;
		case GSL_ESING : 
		{
			fprintf(stderr, "chi_absf :: a non-integrable singularity or other bad integrand behavior was found in the integration interval.\n");
		} break;
		case GSL_EDIVERGE : 
		{
			fprintf(stderr, "chi_absf :: the integral is divergent, or too slowly convergent to be integrated numerically.\n");
		} break;
		case 0 :
		{
			;/*printf("chi_absf :: Integration succesful\n");*/
        
		} break;
		default:
		{
			fprintf(stderr, "chi_absf :: Unknown error in the integration.\n");
		}
	}

	gsl_integration_workspace_free(w);
	return result/(tau_max-tau_th);

}

double Bn_int(double chi0, void *int_params)
{
	double *params = (double *) int_params;
	
	double mu = params[1];
	
	double chi_min,chi_max, t1,t2,t3, t4, t5;

	t4 =  0.5+mu-kappa;
	t5 =  1.+2.0*mu;
	mu += 0.5;

	if(chi0<params[0])
	{
		chi_min=chi0;
		chi_max=params[0];
	}
	else
	{
		chi_min=params[0];
		chi_max=chi0;
	}

	/*mu is mu+0.5*/
	t1 = pow(chi0, -(1.0+kappa)) * exp(-0.5*chi0);
	t2 = exp(-0.5*chi_min) * pow(chi_min,mu) * gsl_sf_hyperg_1F1(t4,t5,chi_min);
	/*whittaker MU (kappa, mu, energy, e_cyc)
	// 	printf("t5 %g\n", t5);*/
	t3 = exp(-0.5*chi_max) * pow(chi_max,mu) * gsl_sf_hyperg_U(t4,t5,chi_max);
	/*whittaker W (kappa, mu, energy, e_cyc)*/

	return t1*t2*t3;
}

double Bn_num(double chi, int n)
{

	double params[2];
	double result, error;
	
	gsl_function F;
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	int err_code;

	params[0]=chi;
	params[1]=mu_vec[n];

	F.function = &Bn_int;
	F.params = params;

	err_code = gsl_integration_qagiu (&F, chi_abs, TOL, TOL, 1000, w, &result, &error);


	switch (err_code)
	{
		case GSL_EMAXITER : 
		{
			fprintf(stderr, "Bn :: the maximum number of subdivisions was exceeded\n");
		}
		break;
		case GSL_EROUND : 
		{
			fprintf(stderr, "Bn :: cannot reach tolerance because of roundoff error, or roundoff  error was detected in the extrapolation table\n");
		} break;
		case GSL_ESING : 
		{
			fprintf(stderr, "Bn :: a non-integrable singularity or other bad integrand behavior was found in the integration interval.\n");
		} break;
		case GSL_EDIVERGE : 
		{
			fprintf(stderr, "Bn :: the integral is divergent, or too slowly convergent to be integrated numerically.\n");
		} break;
		case 0 :
		{
			;/*printf("Bn :: Integration succesful\n");*/
        
		} break;
		default:
		{
			fprintf(stderr, "Bn :: Unknown error in the integration.\n");
		}
	}

	gsl_integration_workspace_free(w);
	

	return result;

}



int L_M_e(double z, double mu, gsl_sf_result *result)
{
	double t1 = mu - kappa +0.5;
    if(t1 ==0 )
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in L_M_e", GSL_EDOM);
    }
    
    gsl_sf_result result_sf;
    int status = gsl_sf_hyperg_2F2_e(mu+kappa+0.5, t1, 1+2*mu,mu-kappa+1.5,-z, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_2F2_e", status);
    }
    result->val=pow(z, t1)/t1 *  result_sf.val;
    
	return GSL_SUCCESS;
}


int I_M_ser_e(double z, double n, double mu, gsl_sf_result *result)
{
	double t1=mu-kappa+0.5;
	
	double t2=1+2*mu;
	
    double t3,t4,t5,t6,t7;
    
    gsl_sf_result result_sf;
    int status = gsl_sf_gamma_e(t1+n, &result_sf);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    t3=result_sf.val;

    
    status = gsl_sf_gamma_inc_e(t1+n,z, &result_sf);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_inc", status);
    }
    t4=result_sf.val;
    
    status = gsl_sf_fact_e(n, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_fact", status);
    }
    t5=result_sf.val;
    
    status = gsl_sf_poch_e(t1,n, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_poch", status);
    }
    t6=result_sf.val;
    
    status = gsl_sf_poch_e(t2,n, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_poch", status);
    }
    t7=result_sf.val;
    
    result->val=(t3 - t4 )/ t5  * t6 / t7;
    
	return GSL_SUCCESS;
	
}


int I_M_e(double z, double mu, gsl_sf_result *result)
{
	
	int n=0;
    gsl_sf_result result_sf;
    int status = I_M_ser_e(z,(double)n,mu, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in I_M_ser_e", status);
    }
    double curr = result_sf.val;
	double new=curr;
	
	int max_n=1000;
	for(n=1;n<max_n;n++)
	{
        status = I_M_ser_e(z,(double)n,mu, &result_sf);
        if(status)
        {
            result->val=GSL_NAN;
            GSL_ERROR ("Error in I_M_ser_e", status);
        }
		new += result_sf.val;
		if( fabs(new/curr)-1 <= TOL || curr == new)
		{
            result->err=fabs(new-curr);
			break;
		}
		else if (n<999)
		{
			curr=new;
		}
	}
	if(n>=max_n-1)
	{
		fprintf(stderr, "Warning I_M series reached the maximum number of terms\n");
		fprintf(stderr, "Values are %f %f (%f) at iteraction %d\n", curr, new, curr/new -1, n);
		
	}
	
    result->val=new;
	return GSL_SUCCESS;
    
}


int L_W_e(double z, double mu, gsl_sf_result *result )
{
	double e1 = mu-kappa+0.5;
	double e2 = -mu-kappa+0.5;
    gsl_sf_result result_sf;
    int status = gsl_sf_gamma_e(e1,&result_sf);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double g1=result_sf.val;
    
    status = gsl_sf_gamma_e(e2,&result_sf);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double g2=result_sf.val;
    
    status = gsl_sf_gamma_e(1-2*kappa,&result_sf);
    if(status != GSL_SUCCESS && status != GSL_EOVRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double g3=result_sf.val;
    
	double t1 = g1 * g2/g3;
    if (gsl_isinf(g3))
        t1=0;
        
//to avoid domain error for negative integers
    
    status = gsl_sf_gamma_e(-2*mu,&result_sf);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double g4=result_sf.val;
    
    status = gsl_sf_hyperg_2F2_e(mu+kappa+0.5, e1, 1+2*mu,mu-kappa+1.5,-z, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_2F2_e", status);
    }
    double h1=result_sf.val;
    
    
    double t2 = g4/g2*pow(z, e1)/e1 * h1;
	
    status = gsl_sf_gamma_e(2*mu,&result_sf);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double g5=result_sf.val;
    
    status = gsl_sf_hyperg_2F2_e(-mu+kappa+0.5, e2, 1-2*mu,-mu-kappa+1.5,-z, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_2F2_e", status);
    }
    double h2=result_sf.val;
    
    
    double t3 = g5/g1*pow(z, e2)/(-e2) * h2;
		

    result->val=t1-t2+t3;
    return GSL_SUCCESS;
}



int I_W_ser_e(double z, double n, double mu, gsl_sf_result *result)
{
	double t1=mu-kappa+0.5;
// 	
 	double t2=2*mu-n;
    gsl_sf_result result_sf;
    int status = gsl_sf_hyperg_U_e(t1, t2,z,&result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double U=result_sf.val;
    
    status = gsl_sf_gamma_e(t1+n,&result_sf);
    if(status != GSL_SUCCESS && status != GSL_EOVRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double g1=result_sf.val;
    
    result->val=U/pow(z,n)/g1;
    if(gsl_isinf(g1))
        result->val=0;
    
    return GSL_SUCCESS;
}

int I_W_e(double z, double mu, gsl_sf_result *result)
{
    gsl_sf_result result_sf;
    int status = gsl_sf_gamma_e(mu-kappa+0.5, &result_sf);
    if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_gamma_e", status);
    }
    double g1=result_sf.val;
    
	double norm = pow(z, mu -kappa-0.5) *exp(-z) * g1;
	
	double n=0;
    
    status=I_W_ser_e(z,n,mu, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in I_W_ser_e", status);
    }
    double curr = result_sf.val;
	double new=curr;

	for(n=1;n<MAX_TERM;n++)
	{
        status=I_W_ser_e(z,n,mu, &result_sf);
        if(status)
        {
            result->val=GSL_NAN;
            GSL_ERROR ("Error in I_W_ser_e", status);
        }
		new += result_sf.val;
		if( fabs(new/curr)-1 <= TOL || curr == new)
		{
            result->err=fabs(new-curr);
			break;
		}
		else
		{
			curr=new;
		}
	}
	if(n>=MAX_TERM-1)
	{
		fprintf(stderr, "Warning I_W series reached the maximum number of terms\n");
	}
	
    result->val=norm*new;
    return GSL_SUCCESS;
}

int II_M_e(double z, double mu, gsl_sf_result *result)
{
    
    int status =GSL_SUCCESS;
	if(z<20)
    {
        status=L_M_e(z,mu, result);
        if(status)
        {
            GSL_ERROR ("Error in L_M_e", status);
        }
        
    
    }
	else
    {
        status=I_M_e(z,mu, result);
        if(status)
        {
            GSL_ERROR ("Error in I_M_e", status);
        }

    }
    
    return GSL_SUCCESS;

}

int II_W_e(double z, double mu, gsl_sf_result *result)
{
    int status=GSL_SUCCESS;
	if(z<15)
    {
        status=L_W_e(z,mu, result);
		
    }
	else
    {
		status= I_W_e(z,mu, result);
    }
    
    if(status)
    {
        GSL_ERROR ("Error in I_W_e or L_W_e", status);
    }
    
    return status;
    

}

int Bn_e(double chi, int n, gsl_sf_result *result)
{
    gsl_sf_result result_sf;
    int status=GSL_SUCCESS;
    
 	double mu = mu_vec[n];
	
    double t2,t3, t4, t5, h1;

	t4 =  0.5+mu-kappa;
	t5 =  1.+2.0*mu;
    status= gsl_sf_hyperg_1F1_e(t4,t5,chi, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_1F1", status);
    }
    h1=result_sf.val;
    
	if(chi<chi_abs)
	{
        status=II_W_e(chi,mu, &result_sf);
        if(status)
        {
            result->val=GSL_NAN;
            GSL_ERROR ("Error in II_W_e", status);
        }
        t3 = result_sf.val;
		mu += 0.5;
        
		t2 = exp(-0.5*chi) * pow(chi,mu) * h1;
        
        result->val=t3* t2;
        
        return GSL_SUCCESS;
	}


    status=II_M_e(chi,mu, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in II_M_e", status);
    }
    double i0 = result_sf.val;

    status=II_M_e(chi_abs,mu, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in II_M_e", status);
    }
    
	i0 -= result_sf.val;
    
    
    status=II_W_e(chi, mu, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in II_W_e", status);
    }
    double i1 = result_sf.val;

	mu += 0.5;
	/*mu is mu+0.5*/
    
	t2 = exp(-0.5*chi) * pow(chi,mu) * h1;
	/*whittaker MU (kappa, mu, energy, e_cyc)*/
    
    status=  gsl_sf_hyperg_U_e(t4,t5,chi,&result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in gsl_sf_hyperg_U_e", status);
    }
    h1=result_sf.val;
    
	t3 = exp(-0.5*chi) * pow(chi,mu) * h1;
	/*whittaker W (kappa, mu, energy, e_cyc)*/
    result->val=t3*i0+t2*i1;
    return GSL_SUCCESS;
}


int ser_FF_factor_e(double chi, double n, double *params, gsl_sf_result *result)
{

    gsl_sf_result result_sf;
    int status=Bn_e(chi, (int)n, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in Bn_e", status);
    }
    double t1 = result_sf.val;
    
	result->val= FF_ser_vec[(int)n]*t1;
    
    return GSL_SUCCESS;

}

int ser_FF_e(double chi,double *params, gsl_sf_result *result)
{
	double n=0;
    gsl_sf_result result_sf;
    int status=ser_FF_factor_e(chi,n,params, &result_sf);
    if(status)
    {
        result->val=GSL_NAN;
        GSL_ERROR ("Error in ser_FF_factor_e", status);
    }
    double curr = result_sf.val;
	double new=curr;

	for(n=1;n<MAX_TERM;n++)
	{
        status=ser_FF_factor_e(chi,n,params, &result_sf);
        if(status)
        {
            result->val=GSL_NAN;
            GSL_ERROR ("Error in ser_FF_factor_e", status);
        }
        new += result_sf.val;
		if( fabs(new/curr)-1 <= TOL || curr == new)
		{
            result->err=fabs(new-curr);
			break;
		}
		else
		{
			curr=new;
		}
	}
	if(n>=MAX_TERM-1)
	{
		fprintf(stderr, "Warning FF series reached the maximum number of terms\n");
	}
	
    result->val=new;
	return GSL_SUCCESS;
}

double FF(double e, double *params)
{
	/*chi abs is computed separately only once*/

	double Te = params[6];
	double chi = e / Te;
	double norm =  norm_FF * pow(e*kev2erg,kappa-2)*exp(-0.5*e/Te);
    gsl_sf_result result_sf;
    int status=GSL_SUCCESS;
    status=ser_FF_e( chi, params, &result_sf);
    if (status) {
        gsl_error ("Error in ser_FF_e", __FILE__, __LINE__, status);
        //waitkey();
        result_sf.val=0;
    }
    return norm * result_sf.val;

}

int init_series(double *params, double DT)
{
	double zeta = params[2];
	double delta = params[3];
	double t1=0,g1=0,g2=0;
	int n=0;
    gsl_sf_result result_sf;
    int status=GSL_SUCCESS;
    char message[1024];
	for(n=0;n<MAX_TERM;n++)
	{

		lambda_vec[n] = lambdaf(n);
        status=muf_e( delta,  zeta, n, &result_sf);
        if(status)
        {
            mu_vec[n] = GSL_NAN;
            GSL_ERROR ("Error in muf_e", status);
            
        }
        mu_vec[n] = result_sf.val;
		
        status= Xn_e(zeta,n, &result_sf);
        if(status)
        {
            Xn_vec[n] = GSL_NAN;
            GSL_ERROR ("Error in Xn_e", status);
            
        }
        Xn_vec[n] = result_sf.val;
		
        
        status= An_e(zeta,n, &result_sf);
        if(status)
        {
            An_vec[n] = GSL_NAN;
            GSL_ERROR ("Error in An_e", status);
            
        }
        An_vec[n] = result_sf.val;
// 		An_vec_sp[n] = An_sp(zeta, n, DT);

        status= gn_e(tau_th, n, &result_sf);
        /*This is a typical error, catched above*/
        if(status == GSL_EUNDRFLW )
        {
            
            sprintf(message, "Underflow in gn_e at values %g %d",tau_th, n );
            gsl_error (message, __FILE__, __LINE__, status);
            gn_vec[n] =0;
        }
        else if (status)
        {
            gn_vec[n] = GSL_NAN;
            
            sprintf(message, "Error in gn_e at values %g %d",tau_th, n );
            GSL_ERROR (message, status);
        }
        else
        {
            gn_vec[n] = result_sf.val;
        }

		/*cyc_ser_vec*/
        status= gsl_sf_gamma_e(mu_vec[n]-kappa+0.5, &result_sf);
        if(status != GSL_SUCCESS && status != GSL_EUNDRFLW)
        {
            
            if(status == GSL_EOVRFLW)
            {
                sprintf(message, "Overflow error in computing series due to gsl_sf_gamma_e at element %d with argument %g",
                        n, mu_vec[n]-kappa+0.5);
                gsl_stream_printf ("ERROR", __FILE__, __LINE__, message);
                return status;
            }
            else
            {
                GSL_ERROR (message, status);
            }
            return status;
            /*Necessary, because, we have suppressed the overflow error output*/
            
        }
        g1=result_sf.val;
        
        status= gsl_sf_gamma_e(1.0+2*mu_vec[n], &result_sf);
        if(status == GSL_EOVRFLW)
        {
            t1=0;
        }
        else if( status == GSL_SUCCESS)
        {
            g2=result_sf.val;
            double g_n=gsl_sf_gamma(n+0.5);
            if(gsl_isinf(g_n))
                t1=0;
            else
                t1 = g1 * gsl_sf_fact(n) / ( g2 * g_n);
        }
        else
        {
            
            /*Underflow for g2, causes t1 to go to infinity, we need to catch the problem*/
            
            sprintf(message, "Error in gamma at values %g %d %g",1.0+2*mu_vec[n], n , result_sf.val);
            g2=GSL_NAN;
            GSL_ERROR (message, status);
            
        }
        
        
        
        
        
		cyc_ser_vec[n] =  t1 * Xn_vec[n] * An_vec[n];
	
		cyc_ser_vec_sp[n] = t1 * Xn_vec[n] * An_sp(zeta, (double)n, DT);

		Green_ser_vec[n] = t1 * Xn_vec[n];

		FF_ser_vec[n] = cyc_ser_vec[n];

	}
    
	return GSL_SUCCESS;
}

int init_norm(double *params)
{
    int status=GSL_SUCCESS;
	double zeta = params[2];
	double delta = params[3];
	double e_cyc = params[4];
	double Mdot = params[5];
	double Te = params[6];
	double dist = 4 * M_PI *  gsl_pow_2(params[8] *1e3*GSL_CONST_CGSM_PARSEC);
	double t1 = zeta*zeta*sqrt(gsl_pow_3(alpha)*w);
	double t2 = Mdot  / sigma_med;
	t2 *= t1;
// 	printf("\t%g %g %g\n", t1, t2, sigma_med);

	norm_FF    = 2.80e-12 * t2 / pow(Te*kev2erg,kappa+0.5)  * kev2erg;
	norm_FF/=dist;

	//printf("\t%g %g %g\n", t1, t2, kev2erg);
	norm_cyc   = Hf(e_cyc/Te) * ((3.43e-16* kev2erg) * t2 )  ;
	norm_cyc/=dist;

	norm_Green = 3.0 * delta * Te * sqrt(2.0) * t1;
	

	//printf("\t%g %g %g\n", norm_FF, norm_cyc, norm_Green);

	norm_BB = 2.0*gsl_pow_2(M_PI/GSL_CONST_CGSM_SPEED_OF_LIGHT)/gsl_pow_3(GSL_CONST_CGSM_PLANCKS_CONSTANT_H/ GSL_CONST_CGSM_ELECTRON_VOLT /1e3);
	norm_BB/=dist;

    return status;
}

void my_gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
    if(gsl_errno != GSL_EUNDRFLW && gsl_errno != GSL_EOVRFLW)
    {
        gsl_stream_printf ("ERROR", file, line, reason);
        fflush (stdout);
    }
    else
    {
        ;
        /*printf("Suppressed error\n");
        gsl_stream_printf ("ERROR", file, line, reason);
         */
    }
    
}


void beckerwolff(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux)
{
    
    /*gsl_error_handler_t *old_error_handler=*/
    gsl_set_error_handler (my_gsl_error);
    
	int i=0;
	
 
	/*params[8] is used to store energy in the integration of BB*/
	
	double params[9];
	double e=0, t1=0,t2=0,t3=0;
    int calc_flag = 0;
    char flag_nan=0;
	
    static char *parameter_names[9]={"R_star", "M_star", "Csi", "delta", "B", "Mdot", "Te", "r0", "D"};

	/*local copy of parameters in cgs units*/

	params[0]=parameter[0] *1e5; /*R_star from km to cm*/
	params[1]=parameter[1] *GSL_CONST_CGSM_SOLAR_MASS; /*//M_star*/

	/* This patch is to avoid Domain Errors in SF*/
	if(parameter[2]==floor(parameter[2]) || parameter[2] == ceil(parameter[2]))
	{
		params[2]=parameter[2]+1e-6*parameter[2]; /* zeta*/
	}
	else
		params[2]=parameter[2]; /* zeta*/

	if(parameter[3]==floor(parameter[3]) || parameter[3] == ceil(parameter[3]))
	{
		params[3]=parameter[3]+1e-6*parameter[3]; /* delta*/
	}
	else
		params[3]=parameter[3]; /* delta*/
	

	if( fabs(params[3] - params[2]) == floor(params[3] - params[2]))
		params[3]*=1.000001;
	
	/*This is not necessary*/

	/*end of patch*/

	params[4]=parameter[4]*11.57; /*//e_cyc from  B*/
	params[5]=parameter[5]*1e17; /*// Mdot*/
	params[6]=parameter[6]; /*// Te*/
	params[7]=parameter[7]*1e2; /*// r0*/
	params[8]=parameter[8]; /*//D*/
    
    

	/*
	for(i=0;i<9;i++) printf("%g\t", params[i]);
	printf("\n About to compute model for %d energy bins\n", Nflux);
	*/

	/*Global variables*/
	alpha      = alphaf(params[2],params[0],params[1]);
	w          = wf(params[2]);
	kappa      = kappaf(params[3]);
	sigma_par  = sigma_parf( params[7],  params[5], params[2]);
	sigma_med  = alpha/3./(params[6]/electron_mass_kev) * sigma_par/params[3];
	tau_th     = tau_thf(params[5],  params[7],  params[0],  params[1],  params[2], &T_th);
	tau_max    = tau_maxf(params[5],  params[7],  params[0],  params[1],  params[2]);
	chi_abs    = chi_absf(params);

	if (init_series(params,-parameter[10]))
    {
        printf("Error in initializing series\n");
        flag_nan=1;
        calc_flag=0;
        
    }
	if (init_norm(params))
    {
        printf("Error in initializing norm\n");
        flag_nan=1;
        calc_flag=0;
    }
	
	/*
	makes an equally spaced grid in energy with 20 points per decade
	and computes the model on it.
	Limit the model in the range 0.5 - 200 keV 
	*/
	static int n_points=54;
	static double step = 0.0490955;/*/(max_energy - min_energy)/((double)n_points -1);*/

	static int n_params=11;

	static double static_parameter[11];
	static int n_call=0;
	static double static_spec[54];
	static double static_en[54];

	static double static_log_en_min = 10;
	static double static_log_en_max = -10;

	

	if(n_call==0)
	{
		for(i=0;i<n_points;i++)
		{
			static_en[i] = -0.301030 + i * step; //min_energy + i* step;
			static_spec[i] = 0;
		}
		for(i=0;i<n_params;i++)
		{
			static_parameter[i] = 0;
		}
	}

	
	/*Store parameters*/
	for(i=0;i<n_params;i++)
	{
		
  		if(static_parameter[i] != parameter[i])
		{
			calc_flag =1;
		}
		
		static_parameter[i] = parameter[i];
	}

	double log_en_min = log10(energy[0]);
	double log_en_max = log10(energy[Nflux]);
	
/*
 * double ec = parameter[11];
 * double sc = parameter[12];
 * double dc = parameter[13];
 */
    /* Samity check for parameters*/
    if(parameter[7] <=0 || parameter[5] <=0 || parameter[2] <=0)
    {
        printf("\nError :: Parameters Mdot, csi, r0 must be positive, please correct your allowed limits !\n");
        for(i=0;i<9;i++) printf("%g\t", params[i]);
        printf("\n\n");
        waitkey();
        for(i=0;i<n_points;i++)
            static_spec[i]=0;
        calc_flag=0;
    }
    
    
	if(calc_flag == 1 && flag_nan == 0)
	{

		/* Output aux parameters

		double beta_th = 7.86e10 * pow(T_th, -7./4.) * params[5] * pow(params[7], -1.5) / GSL_CONST_CGSM_SPEED_OF_LIGHT;

		double ratio_sigma=sigma_par/GSL_CONST_CGSM_THOMSON_CROSS_SECTION;

		printf("alpha, w ,kappa, sigma_par/sigma_T, sigma_med/sigma_T, T_th, v_th/c, tau_th, tau_max, chi_abs\n");
 		printf("%g %g %g %g %g %g %g %g %g %g\n", alpha, w ,kappa, sigma_par/GSL_CONST_CGSM_THOMSON_CROSS_SECTION, sigma_med/GSL_CONST_CGSM_THOMSON_CROSS_SECTION,T_th, beta_th, tau_th, tau_max, chi_abs); // *kelvin2keV
		double z_th=zf(tau_th, params[2], params[7]);
		double z_max=zf(tau_max, params[2], params[7]);
		
		printf("T_th = %g keV\n", T_th*kelvin2keV);
		printf("Xi, delta, M , R, r0, Mdot, Te, B\n");
		printf("%g %g %g %g %g %g %g %g\n", params[2], params[3], params[1], params[0], params[7], params[5], params[6], params[4]/11.57);  ///kelvin2keV
		//printf("E_c = %g s_c = %g d_c =%g\n", ec, sc, dc);
		double z_sp=log(7./3.) /2/sqrt(3) * params[7] * sqrt(GSL_CONST_CGSM_THOMSON_CROSS_SECTION/sigma_par);
		double tau_sp=sqrt(log(7./3.) /sqrt(3)/alpha/params[2]);
		
		double beta_sp = pow(ratio_sigma, 0.25) * sqrt( 2*alpha*z_sp/params[2]/params[7]);
		double beta_max = pow(ratio_sigma, 0.25) * sqrt( 2*alpha*z_max/params[2]/params[7]);
	
		printf("z_th=%g    z_sp=%g    z_max=%g\n", z_th, z_sp, z_max);
		printf("tau_th=%g  tau_sp=%g  tau_max=%g\n", tau_th, tau_sp, tau_max);
		printf("rho_th=%g  rho_sp=%g  rho_max=%g\n", rhof(tau_th,params[2],params[7],params[5]),  rhof(tau_sp,params[2],params[7],params[5]),  rhof(tau_max,params[2],params[7],params[5]) );
		printf("beta_th=%g    beta_sp=%g    beta_max=%g\n", beta_th, beta_sp, beta_max);
		printf(" CGS Units !\n");
 		//printf("About to compute detailed model for %d energy bins at call %d\n", n_points, n_call);
		///*/
		
#if FILEOUTPUT
        
		char fname[512];
		sprintf(fname, "output-%03d.qdp", spectrum);
		FILE *output=fopen(fname, "w");
		//fprintf(output, "read\ncpd /xw\nr y 1e-7 10\nr x 1.7 110\nlw 2\n");
		fprintf(output, "read\ncpd /xw\nr y -7 1\nr x 0.25 2\nlw 2\n");
		fprintf(output, "grid x 7,2\ngrid y 8,2\ngrid on\n");
		//fprintf(output, "log x on\nlog y on\n");
		//fprintf(output, "lab x Energy [keV]\nlab y dN/dE [s\\u-1\\dcm\\u-2\\dkeV\\u-1\\d]\n");
		fprintf(output, "lab x log10(Energy [keV])\nlab y log10(dN/dE [s\\u-1\\dcm\\u-2\\dkeV\\u-1\\d])\n");
		fprintf(output, "time off\nlab T Becker and Wolff (2007)\n");
		fprintf(output, "lab file \\gc=%4.2f \\gd=%4.2f B=%4.2fe12 G Mdot=%2ge17 g/s T\\de\\u=%3.1f keV r\\d0\\u=%4.1f m\n",
		parameter[2], parameter[3], parameter[4], parameter[5], parameter[6], parameter[7]);
#endif
        
        gsl_sf_result result_sf;
        int status=GSL_SUCCESS;

		for(i=0;i<n_points;i++)
		{
			e = static_en[i];
			if( (e < static_log_en_min ||  e > static_log_en_max ) && n_call )
				continue;
			e = pow(10,e);
			//double de = (e - ec)/sc;
			//double abs_c = 1.;
			//if(fabs(de) <= 5)
			//	abs_c -= dc / sc / sqrt(2.*M_PI) *exp(- 0.5*de*de);
			/*parameter 10 may disable the computation of BB emission and scale it*/
			static_spec[i]=0;
			if(parameter[9] > 0)
			{
				t1 = BB(e,params) * parameter[9];// * abs_c; 
				if(!gsl_finite(t1)|| t1>1e6)
				{
                    flag_nan=1;
                    //printf("\nBB not finite\n");
					/*for(i=0;i<9;i++) printf("%g\t", params[i]);
					printf("\nEnd BB\n");
					waitkey();*/
				}
				static_spec[i] += t1;
				t1=log10(t1);
					
			}
			/*parameter 11 may disable the cmputation of CYC emission and scale it*/
			if(parameter[10] > 0)
			{
                status=cyc_e(e,params, &result_sf) ;// * abs_c;
                if(status)
                {
                    gsl_error ("Error in cyc_e", __FILE__, __LINE__, status);
                    t2=0;
                    flag_nan=1;
                }
                else
                {
                    t2 = result_sf.val * parameter[10];
                }
                
                
				if(!gsl_finite(t2)|| t2>1e6)
				{
                    flag_nan=1;
					/*printf("Cyc not finite\n");
					for(i=0;i<9;i++) printf("%g\t", params[i]);
					printf("\nEnd Cyc\n");
					waitkey();*/
				}
				static_spec[i] +=t2;
				t2=log10(t2);
			}
			else if(parameter[10] <0)
			{
				t2 = cyc_sp(e, params, -parameter[10]);// this is for a small layer of material around the sonic point
                
				if(!gsl_finite(t2)|| t2>1e6)
				{
                    flag_nan=1;
					/*printf("Cyc SP not finite\n");
					for(i=0;i<9;i++) printf("%g\t", params[i]);
					printf("\nEnd Cyc\n");
					waitkey();*/
				}
				static_spec[i] +=t2;
				t2=log10(t2);

			}
			t3 = FF(e,params) ;// * abs_c; 
			if(!gsl_finite(t3) || t3>1e7)
			{
                flag_nan=1;
				/*
                 printf("FF not finite\n");
                
				for(i=0;i<9;i++) printf("%g\t", params[i]);
				printf("\nEnd FF\n");
				waitkey();
                 */
			}
			
			static_spec[i] += t3;
			t3=log10(t3);

			if(static_spec[i]>0)
				static_spec[i]=log10(static_spec[i]);
			else
				static_spec[i]=GSL_NAN;
			
			if(!gsl_finite(static_spec[i]))
			{
                flag_nan=1;
				/*
                 printf("Log not finite %g %g \n", e, static_spec[i]);
				
                for(i=0;i<9;i++) printf("%g\t", params[i]);
				printf("\nEnd Log\n");
				waitkey();
                */
			}
#if FILEOUTPUT
			fprintf(output, "%g\t%g\t%g\t%g\t%g\n", log10(e), t1, t2, t3, static_spec[i] );
#endif
		} /*end of energy cycle*/
#if FILEOUTPUT
		fclose(output);
#endif

	} /*end if calc_flag*/

	

	if(log_en_min < static_log_en_min)
	{
		static_log_en_min = log_en_min;
	}

	if(log_en_max > static_log_en_max)
	{
		static_log_en_max = log_en_max;
	}

	n_call++;


    /*Stops her in case of NAN*/
    if(flag_nan)
    {
        printf("NAN values inserted in the spectrum due to parameter space boundary violations\n Model Parameters: ");
        for(i=2;i<9;i++) printf("%s=%g\t", parameter_names[i], parameter[i]);
        printf("\nEnd\n");
        for(i=0;i<Nflux;i++)
        {
            flux[i] = GSL_NAN;
        }
        return;
    }
	/*
	prepares spline interpolatoin
	*/
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, n_points);
    gsl_spline_init (spline, static_en, static_spec, n_points);

	/*
	makes a numerical integration
	using five points per bin
	*/

	for(i=0;i<Nflux;i++)
	{
		if(energy[i] < 0.5 || energy[i+1] > 200.0)
			flux[i] = 0;
		else
		{
			flux[i] = integralbw(energy+i, spline, acc);
            /*
			if(!gsl_finite(flux[i]) || flux[i]>1e7)
			{
				printf("model not finite %g %g\n", energy[i], flux[i]);
				
                
				for(i=0;i<9;i++) printf("%g\t", params[i]);
				printf("\nEnd\n");
				waitkey();
                 
			}
             */
		}
	}

	/*
	frees variables
	*/

    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

/*
	printf("\n stop interpolation\n");
	printf("Done\n");*/

  	

	return;
}

double integralbw(const double *e_a, gsl_spline *spline, gsl_interp_accel *acc)
{
        /*5-points Newton-Cotes integration formula, we assume energy bins are not too large*/
        double e1 = e_a[0];
        double h  = (e_a[1]-e1)/4.;
	
	
        double pp[5];
        int i;
        for(i=0;i<5;i++)
        {
            pp[i]=pow(10,gsl_spline_eval (spline, log10(e1+h*i), acc));
        }

        return 2.0*h*(7.0*(pp[0]+pp[4])+32.0*(pp[1]+pp[3])+12.0*pp[2])/45.0;
}


/*
 * You should run the following command to link the gsl without recompiling XSPEC
g++ -shared -o libbwmod.so  bw.o   bwmod.o gausabs.o bwmodFunctionMap.o /usr/lib/libgsl.a /usr/lib/libgslcblas.a
or
g++ -shared -o libbwmod.so  bw.o   bwmod.o gausabs.o bwmodFunctionMap.o `gsl-config --libs`
 
 Or 
 
 see the section
 
 Third-Party Libraries In Local Models Build
 
 In the page
 https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixLocal.html
*/

