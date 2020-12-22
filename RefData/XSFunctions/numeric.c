#include "numeric.h"

/* ############################################################################# */
/* Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
/* ############################################################################# */

double **dmatrix(long nrl, long nrh, long ncl, long nch)

{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate pointers to rows */
  m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
  /* check_alloc(m); */
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) {
    printf("Memory allocation failure in routine dmatrix");
    exit(1);
  }

  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
  /* return pointer to array of pointers to rows */
  return m;
}

/* ###################################################################### */
/* Allocate a double vector with subscript range v[nl..nh] */
/* ############################################################################# */

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) {
    printf("Memory allocation failure in routine dvector");
    exit(1);
  }

  return v-nl+NR_END;
}

/* ###################################################################### */

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)

{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

/* ###################################################################### */

void free_dvector(double *v, long nl, long nh)

{
  free((FREE_ARG) (v+nl-NR_END));
}

/* ############################################################################# */

double interp_funct(double *array_x, double *array_y, int N, double x) {

  int i;
  double m, q, interp_y;

  for (i=0; i<N; i++) {

    if (array_x[i] <= x && array_x[i+1] >=x) {

      m=(array_y[i+1]-array_y[i])/(array_x[i+1]-array_x[i]);	
      q=array_y[i]-m*array_x[i];
      interp_y=m*x+q;
      return(interp_y);
    }
    
  }

  return 0.0;
}

/* ##################################################################*/
/* Logarithm of the gamma function */
/* ##################################################################*/

double gammln(double xx)

{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

/* ##################################################################*/
/* Riemann zeta function */
/* ##################################################################*/


double zeta(double xx) {

  int n;
  double value=0;
  
  for (n=1; n<=10; n++) {
    value=value+1/(pow(n,xx));
  }

  return(value);

}
