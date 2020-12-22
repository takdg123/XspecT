/* H+					
				
	Title:	ismatten.c
	Author: Pat Jelinsky
	Date:	10/20/88
	Synopsis: ismatten containes several subroutines that are useful
		  in the calculation of the attenuation from the
		  interstellar medium
		  "tauhe" calculates the optical depth due to neutral helium
		  "tauh" calculates the optical depth for any hydrogenic atom
		  "atten" calculates the transmission of neutral hydrogen,
		          neutral helium, and once ionized helium
	Keywords:  scs_section, major_pgm_name, subject	
	Revisions:		
	05/05/29	patj	changed he cross sections to a new
				data provided by Stefan Vennes.
	03/16/89	patj	changed he cross sections to a linear
				interpolation of measured values.
	06/08/90	patj	split out subroutines to a separate file
        07/22/92        toddr   added he autoionization resonances  
   H-						*/		

/* U+					
 	Usage:	tau = tauhe(wav,heicol);
	Input:	double	wav 		wavelength (angstroms)
		double	heicol		neutral helium column density (cm**2)
	Output: double	tau		optical depth
   U-						*/

/* U+					
	Usage:	tau = tauh(wav,hcol,zee);
	Input:	double	wav 		wavelength (angstroms)
		double	hcol		neutral hydrogen column density (cm**2)
		double	zee		number of protons in atom
	Output:	double	tau		optical depth
   U-						*/

/* U+					
 	Usage:	transmission = atten(wav,hcol,heicol,heiicol);
	Input:	double	wav		wavelength (angstroms)
		double	hcol		neutral hydrogen column density (cm**2)
		double	heicol		neutral helium column density (cm**2)
		double	heiicol 	ionized helium column density (cm**2)
	Output: double	transmission	transmission of ism
   U-						*/

#include	<stdio.h>
#include	<math.h>


/* E+
	tauh returns the ism optical depth 
	for any hydrogenic atom as per Spitzer (Physical processes in
	the interstellar medium)  Page 105.
	The inputs are the wavelength, wav, in angstroms; the
	column density, hcol, in cm**-2; and the charge of the
	ion, zee.
   E-									*/
double	tauh(wav,hcol,zee)
double	wav;
double	hcol;
double	zee;
{
	double	z,sigma,tau,ratio;
	
	ratio = zee * zee * wav / 911.75;
	if(ratio < 1.0) {
		z = sqrt(ratio/(1.0-ratio));
		sigma = 3.44e-16 * pow(ratio,4.0) * exp(-4.0 * z * atan(1/z))/
		((1.0 - exp(-6.283185308*z))*zee*zee);
		tau = hcol * sigma;
	}
	else {
		tau = 0.0;
	}
	return(tau);
}

/* E+
	atten returns the attenuation due to the ISM
	given the wavelength, wav, in angstroms; the neutral hydrogen
	column density, hcol, in cm**-2; the neutral helium column
	density, heicol, in cm**-2; and the singly ionized helium
	column density, heiicol, in cm**-2.
   E-							*/

double	atten(wav,hcol,heicol,heiicol)
double	wav;				/* Wavelength (Angstroms) */
double	hcol;				/* Neutral Hydrogen column Density */
double	heicol;				/* Neutral Helium column Density */
double	heiicol;			/* Ionized Helium column Density */
{
	double	tau,trans;
	extern	double	tauh();
	extern	double	tauhe();

	tau = tauh(wav,hcol,1.0) + tauh(wav,heiicol,2.0) + tauhe(wav,heicol);
	trans = exp(-tau);
	return(trans);
}

/*
	From experimental data compiled by Marr & West (1976)
	Atomic Data and Nuclear Data Tables, 18, 497
	(from Frits Paerels 12/90) (then from Stefan Vennes 4/92)
						*/

/*
	polynomial coeffients for He I cross section to use for wavelengths 
	greater than or equal to 46 A.
									*/
static double c1[] = {-2.953607e+01, 7.083061e+00, 8.678646e-01,-1.221932e+00,
		 4.052997e-02, 1.317109e-01, -3.265795e-02, 2.500933e-03};

/*
	polynomial coeffients for He I cross section to use for wavelengths 
	less than 46 A.
									*/
static double c2[] = {-2.465188e+01, 4.354679e+00, -3.553024e+00, 5.573040e+00,
		 -5.872938e+00, 3.720797e+00, -1.226919e+00, 1.576657e-01};

/*      parameters of autoionization resonances for 4 strongest for helium
        Numbers are from Oza (1986), Phys Rev. A, 33, 824 -- nu and gamma
        and Fernley et al., J. Phys. B., 20, 6457, 1987 -- q
                                                                */

static double fano_q[]  = {2.81, 2.51, 2.45, 2.44};
static double fano_nu[] = {1.610, 2.795, 3.817, 4.824};
static double fano_gamma[] = {2.64061e-03, 6.20116e-04, 2.56061e-04,
				1.320159e-04 };

/* E+
	tauhe returns the neutral helium optical depth. The
	cross sections are from experimental data compiled by Marr & West (1976)
	Atomic Data and Nuclear Data Tables, 18, 497
	Code based on fortran code from (Frits Paerels 12/90)
	I got the code from Stefan Vennes (4/1992)
	the inputs are the wavelength, wav, in angstroms, 
	and the neutral helium column density, hecol, in cm**-2.
	tauhe returns the optical depth from a log-log polynomial 
	interpolation of the measured cross sections.
   E-								*/
double tauhe(lambda,hecol)
double	lambda;	/* wavelength in Angstroms	*/
double	hecol;	/* column density of neutral helium (cm**-2) */
{
	double	x;
	double	y;

	extern double fano();

	if(lambda > 503.97) {	/* if wavelength is above ionization limit */
		return(0.0);	/* then there is no absorption	*/
	}

	x=log10(lambda);	/* polynomial fits use log of lambda */

	if(lambda < 46.0) {	/* if wavelength < 46 use c2 */
      		y=c2[0]+x*(c2[1]+x*(c2[2]+x*(c2[3]+x*(c2[4]+x*(c2[5]+
					x*(c2[6]+x*c2[7]))))));
	}
	else {			/* if wavelength > 46.0 use c1 */
		y=c1[0]+x*(c1[1]+x*(c1[2]+x*(c1[3]+x*(c1[4]+x*(c1[5]+
					x*(c1[6]+x*c1[7]))))));
		y += log10(fano(fano_q[0],fano_nu[0],fano_gamma[0],lambda));
		y += log10(fano(fano_q[1],fano_nu[1],fano_gamma[1],lambda));
		y += log10(fano(fano_q[2],fano_nu[2],fano_gamma[2],lambda));
		y += log10(fano(fano_q[3],fano_nu[3],fano_gamma[3],lambda));
	}
	
	return(hecol*pow(10.0,y));
}


/* E+
	double fano( q, nu, gamma,lambda) returns a Fano line profile for use 
        by tauhe.  The form of this line shape is taken from:
	
	Fernley, Taylor and Seaton; Journal of Physics (B), 20: 6457 
        (1987).
   E-								*/
double fano(q, nu, gamma, lambda)
double	q;	
double	nu;	
double gamma;
double lambda;

{ 
  double epsilon = 911.2671/lambda;   /* energy in rydbergs */
  double esubi  = 3.0 - 1.0/(nu*nu) + 1.807317; 
  double x = 2.0*((epsilon-esubi)/gamma);

  return(((x-q)*(x-q))/(1+x*x));
}
