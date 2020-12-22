
#include "cfortran.h"
#include "XSFunctions/Utilities/xsFortran.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "exsource.h"

#define MMIN(a,b) ((a)<(b)?(a):(b))
#define MMAX(a,b) ((a)>(b)?(a):(b))

#define NCONV 8192
#define HC 12.398

void rgsxsrc
(float *ear,int ne,float *param,int ifl,float *photar, float *photer);

void convlv(float data[],int n, float respns[], int m, int isign,float ans[]);
void four1(float data[], int nn, int isign);
void hunt(float xx[],int n, float x, int *jlo);
void nrerror(char error_text[]);
float *nrvector(long nl, long nh);
void free_nrvector(float *v, long nl, long nh);
void realft(float data[], int n, int isign);
void twofft(float data1[], float data2[], float fft1[], float fft2[], int n);

FCALLSCSUB6(rgsxsrc,RGSXSRC,rgsxsrc,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)



void rgsxsrc
(float *ear,int ne,float *param,int ifl,float *photar, float *photer)
{
  /* this routine will convolve the input spectrum by an *angular* structure
     function, according to the exsource code. this is really a test on how
     XSPEC can fake up an angular distribution function so that matrix 
     generation can be bypassed for quick & dirty applications.
     to do this we need to redistribute the photar[] array onto a grid
     that is uniformly spaced in (lamda) and then convolve by a kernel.
     normally the energy bin range is sufficient to describe the wavelength
     range of the spectrometer. the grid should be uniformly spaced in
     (lamda) because 
    
         delta(lamda)_eff = d/m * sin(alpha) * F/L * phi 

     (beta is held constant).
*/
  static exsource_projection *exs=NULL;
  int i,j,c1,c2;
  float minlam,maxlam,minene,maxene;

  float order=param[0];
  float *lamarr=(float*)malloc((NCONV+1)*sizeof(float));
  float *tmparr=(float*)malloc(NCONV*sizeof(float));

  if (exs==NULL) { /* first invocation */

    exs=(exsource_projection*)malloc(sizeof(exsource_projection));

    if (exs==NULL) {
      fprintf(stderr,"can't allocate exsource projection\n");
      return;
    }
    
    exs->initialized=0;
    exs->fullreshist.values=NULL;
    exs->hist.values=NULL;
    exs->image="RGS_XSOURCE_IMAGE";
    exs->bores="RGS_XSOURCE_BORESIGHT";
    exs->extraction="RGS_XSOURCE_EXTRACTION";
    exs->image_val=NULL;
    exs->bores_val=NULL;
    exs->extraction_val=NULL;

  }

  {
    char s1[2048],s2[2048],s3[2048],s4[2048],s5[2048];
    char line[2048];
    char sIn[10240];
    char *parfile;
    int  nconv;
    FILE *fp;

    parfile = (char*)malloc(2048*sizeof(char));
    parfile = FGMSTR("RGS_XSOURCE_FILE");
    if (strlen(parfile)!=0) {
      if ((fp=fopen(parfile,"r"))!=NULL) {

	while (fgets(line,2048,fp)) {

	  nconv=sscanf(line,"%s %s %s %s %s\n",s1,s2,s3,s4,s5);

	  switch (nconv) {
	  case 5:
	    sprintf(sIn,"%s %s %s %s",s2,s3,s4,s5);
	    break;
	  case 4:
	    sprintf(sIn,"%s %s %s",s2,s3,s4);
	    break;
	  case 3:
	    sprintf(sIn,"%s %s",s2,s3);
	    break;
	  case 2:
	    sprintf(sIn,"%s",s2);
	    break;
	  default:
	    break;
	  }

	  switch (nconv) {
	  case 1:
	    /* do nothing */
	    break;
	  case 5:	  case 4:	  case 3:	  case 2:
	    /* set variable as usual
	       check first to see if the value is being changed. */
	    if ( strcmp(s1,exs->image) == 0 ) {
              if ( exs->image_val != NULL ) {
		if ( strcmp(exs->image_val,sIn) != 0 ) {
		  fprintf(stderr,"value in file %s supercedes previous value!!\n"
		          "setting variable %s to %s\n was %s\n",parfile,s1,sIn,
                          exs->image_val);
                  exs->image_val = malloc(sizeof(sIn));
                  strcpy(exs->image_val,sIn);
	          exs->initialized=0;
                }
              } else {
                exs->image_val = malloc(sizeof(sIn));
                strcpy(exs->image_val,sIn);
	        exs->initialized=0;
              }
            } else if ( strcmp(s1,exs->bores) == 0 ) {
              if ( exs->bores_val != NULL ) {
		if ( strcmp(exs->bores_val,sIn) != 0 ) {
		  fprintf(stderr,"value in file %s supercedes previous value!!\n"
		          "setting variable %s to %s\n was %s\n",parfile,s1,sIn,
                          exs->bores_val);
                  exs->bores_val = malloc(sizeof(sIn));
                  strcpy(exs->bores_val,sIn);
	          exs->initialized=0;
                }
              } else {
                exs->bores_val = malloc(sizeof(sIn));
                strcpy(exs->bores_val,sIn);
	        exs->initialized=0;
              }
            } else if ( strcmp(s1,exs->extraction) == 0 ) {
              if ( exs->extraction_val != NULL ) {
		if ( strcmp(exs->extraction_val,sIn) != 0 ) {
		  fprintf(stderr,"value in file %s supercedes previous value!!\n"
		          "setting variable %s to %s\n was %s\n",parfile,s1,sIn,
                          exs->extraction_val);
                  exs->extraction_val = malloc(sizeof(sIn));
                  strcpy(exs->extraction_val,sIn);
	          exs->initialized=0;
                }
              } else {
                exs->extraction_val = malloc(sizeof(sIn));
                strcpy(exs->extraction_val,sIn);
	        exs->initialized=0;
              }
	    }
	    break;
	  case EOF:
	  default:
	    break;
	  }
	}
	fclose(fp);
      } else {
	fprintf(stderr,"%s is set to %s but can't open that file!\n",
                "RGS_XSOURCE_FILE", parfile);
        return;
      }
    } else {
	fprintf(stderr,
		"RGS_XSOURCE_FILE has not been set - rgsxsrc doing nothing.\n");
        return;
    }
  }

  /* otherwise exs is presumably initialized. */

  if (!lamarr || !tmparr) {
    fprintf(stderr,"can't allocate within rgsxsrc\n");
    return;
  }

  i=NCONV;
  while (i--) tmparr[i] = 0.0;

  minlam=HC/MMAX(ear[0],ear[ne]);
  maxlam=HC/MMIN(ear[0],ear[ne]);

  for (i=0;i<=NCONV;i++) 
    lamarr[i]=minlam+i*((maxlam-minlam)/(float)(NCONV));

  if (1) {
    /* now go through and redistribute */
    c1 = c2 = 0;
    for (i=0;i<ne;i++) {
      minlam=HC/MMAX(ear[i],ear[i+1]);
      maxlam=HC/MMIN(ear[i],ear[i+1]);
      hunt(lamarr-1,NCONV,minlam,&c1);
      hunt(lamarr-1,NCONV,maxlam,&c2);
      if (c1==NCONV) c1--;
      if (c1) c1--; /* since it's numrec */
      if (c2==NCONV) c2--;
      if (c2) c2--; /* since it's numrec */
      
      /* now 
	 lamarr[c1] < minlam < lamarr[c1+1] 
	 or
	 lamarr[c1] > minlam > lamarr[c1+1]  
	 etc.. 
      */

      for (j=MMIN(c1,c2);j<=MMAX(c1,c2);j++) {
	tmparr[j] += photar[i] 
	  * (MMIN(maxlam,MMAX(lamarr[j],lamarr[j+1]))
	     -MMAX(minlam,MMIN(lamarr[j],lamarr[j+1])))/(maxlam-minlam);
      }
    }

    /* ready to convolve lamarr[] in wavelength space */
    
    {
      float *conv_func,*conv,*kernel;
      int   N=NCONV,M=NCONV/2-1;
      kernel=   (float*)malloc(NCONV*sizeof(float));
      conv=     (float*)malloc(2*NCONV*sizeof(float));
      conv_func=(float*)malloc(NCONV*sizeof(float));
      if (!kernel || !conv || !conv_func) {
	fprintf(stderr,"can't allocate within rgsxsrc\n");
	return;
      }
      j=NCONV;
      while (j--) kernel[j] = conv[j] = 0.0;

      /*      (maxlam-minlam)/((float)NCONV) */
      /*	/(15491.867/order * 7500.0/6700.0 * 0.02751) * M/2.0; */

      exs->hist.limit[0] = -1*fabs(lamarr[0]-lamarr[NCONV])*(M)/(2.0*NCONV)
	/(15491.867/order * 7500.0/6700.0 * 0.02751);
      exs->hist.limit[1] = -exs->hist.limit[0];
      exs->hist.nch      = M;

      if (exsource(exs)!=0) return;
      
      j=NCONV;
      while (j--) conv_func[j]=tmparr[j];

      j=0;
      while (j<M) {
	kernel[j]=exs->hist.values[((M-1)/2+j)%M];
	j++;
      }

      /* do convolution */
      convlv(conv_func-1,N,kernel-1,M,1,conv-1);
      j=NCONV;
      while (j--) tmparr[j]=conv[j];
      /* complete. ready to deallocate and redistribute */
      free(conv);
      free(conv_func);
      free(kernel);
    }

    /* ready to redistribute back into photar[] */
    /* clear out the photar[] first */

    for (i=0;i<ne;i++) photar[i]=0.0;
    c1 = c2 = 0;
    for (j=0;j<NCONV;j++) {
      minene=HC/MMAX(lamarr[j],lamarr[j+1]);
      maxene=HC/MMIN(lamarr[j],lamarr[j+1]);
      hunt(ear-1,ne+1,minene,&c1);
      hunt(ear-1,ne+1,maxene,&c2);
      if (c1==ne+1) c1--;
      if (c1) c1--; /* since it's numrec */
      if (c2==ne+1) c2--;
      if (c2) c2--; /* since it's numrec */
      for (i=MMIN(c1,c2);i<=MMAX(c1,c2);i++) {
	photar[i]+=tmparr[j]
	  * (MMIN(maxene,MMAX(ear[i],ear[i+1]))
	     -MMAX(minene,MMIN(ear[i],ear[i+1])))/(maxene-minene);
      }
    }

  } else {
    /* do nothing; perform no calculation. 
       free up allocated arrays and return */
  }

  free(lamarr);
  free(tmparr);
  return;
}

void convlv(float data[],int n,float respns[], int m, int isign, float ans[])
{
	int i,no2;
	float dum,mag2,*fft,*nrvector();

	fft=nrvector(1,2*n);
	for (i=1;i<=(m-1)/2;i++)
		respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
		respns[i]=0.0;
	twofft(data,respns,fft,ans,n);
	no2=n/2;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		} else if (isign == -1) {
 		        if ((mag2=ans[i-1]*ans[i-1]+ans[i]*ans[i]) == 0.0)
				nrerror("Deconvolving at response zero in CONVLV");
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		} else nrerror("No meaning for ISIGN in CONVLV");
	}
	ans[2]=ans[n+1];
	realft(ans,no2,-1);
	free_nrvector(fft,1,2*n);
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(float data[], int nn, int isign)
{
	int n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=2*mmax;
		theta=6.28318530717959/(isign*mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

#undef SWAP

void hunt(float xx[], int n, float x, int *jlo)
{
	int jm,jhi,inc,ascnd;

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if ((x >= xx[*jlo]) == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while ((x >= xx[jhi]) == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while ((x < xx[*jlo]) == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if ((x > xx[jm]) == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
}

void nrerror(char error_text[])
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *nrvector(long nl, long nh)
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in nrvector()");
	return v-nl;
}

void free_nrvector(float *v, long nl, long nh)
{
	free((char*) (v+nl));
}

void realft(float data[], int n, int isign)
{
	int i,i1,i2,i3,i4,n2p3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) n;
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	n2p3=2*n+3;
	for (i=2;i<=n/2;i++) {
		i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n,-1);
	}
}

void twofft(float data1[], float data2[], float fft1[], float fft2[], int n)
{
	int nn3,nn2,jj,j;
	float rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	four1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}
