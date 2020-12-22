/* exsource.c - routine library for generating kernel input arrays for
   computing the effects of extended sources on the line spread function. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#define EXSOURCE_CODE
#include "exsource.h"

#define  min(a,b) (((a)<(b))?(a):(b))
#define  max(a,b) (((a)>(b))?(a):(b))

int exsource (exsource_projection *exs) {
  char *a,*b,*c;
  aspect asp;
  float extraction_halfwidth;

  if (exs->initialized == 0) {
    if ((a=exs->image_val)==NULL) {
      char errstr[2048];
      sprintf(errstr,"%s not set. exsource will do nothing without this.\n"
              "It should be set in RGS_XSOURCE_FILE.\n", exs->image); 
      exsource_complain(errstr);
      return 1;
    }

    if ((b=exs->bores_val)==NULL) {
      char errstr[2048];
      sprintf(errstr,"%s not set. exsource will do nothing without this.\n"
              "It should be set in RGS_XSOURCE_FILE.\n", exs->bores); 
      exsource_complain(errstr);
      return 1;
    }

    if ((c=exs->extraction_val)==NULL) {
      char errstr[2048];
      sprintf(errstr,"%s not set. exsource will do nothing without this.\n"
              "It should be set in RGS_XSOURCE_FILE.\n", exs->extraction); 
      exsource_complain(errstr);
      return 1;
    }

    
    /* extraction_halfwidth is converted to radians. */
    extraction_halfwidth=M_PI*atof(c)/180.0/60.0/2.0;

    fprintf(stderr,
      "-------------------------------\n"
      "RGS XSRC parameters            : (author=arasmus@astro.columbia.edu)\n"
      "rgs xsrc image file            : %s\n"
      "rgs xsrc boresight (ra,dec,pa) : %s\n"
      "rgs xsrc x-disp extraction     : %s arcmin\n"
      "-------------------------------\n",a,b,c);
    str2coords (b,&asp); 

    if (exs->fullreshist.values!=NULL) {
      free(exs->fullreshist.values);
      exs->fullreshist.values=NULL;
    }

    if (exs->hist.values!=NULL) {
      free(exs->hist.values);
      exs->hist.values=NULL;
    }

    /* open the image and read the pixels. */

    if (img2exhist(a,&asp,extraction_halfwidth,&exs->fullreshist)!=0) return 1;
    exs->initialized=1;
  }

  /* rebin the hist */
  if (exs->hist.values != NULL) free(exs->hist.values);
  exs->hist.values=(float*)malloc(exs->hist.nch*sizeof(float));
  hist2hist(&exs->hist,&exs->fullreshist);
  return 0;
}

void hist2hist (hgram *to_hg, hgram *fm_hg) {
  int i,j;
  double chmin,chmax,factor,offset,slope;

  for (j=0;j<to_hg->nch;j++) to_hg->values[j]=0.0;

  offset = to_hg->nch*
    (fm_hg->limit[0]-to_hg->limit[0])/(to_hg->limit[1]-to_hg->limit[0]);
  slope  = (to_hg->nch*(fm_hg->limit[1]-fm_hg->limit[0]))
    /(fm_hg->nch*(to_hg->limit[1]-to_hg->limit[0]));

  for (i=0;i<fm_hg->nch;i++) {
    if (slope>0) {
      chmin=offset+i*slope;
      chmax=offset+(i+1)*slope;
    } else {
      chmax=offset+i*slope;
      chmin=offset+(i+1)*slope;
    }

    /* and distribute the power in fm_hg->values[i] into to_hg->values[] */
    for (j=max(0,floor(chmin));j<min(to_hg->nch,ceil(chmax));j++) {
      factor=(min(chmax,j+1)-max(chmin,j))/(chmax-chmin);
      to_hg->values[j] += factor * fm_hg->values[i];
    }
  }
  if (0) {
    hgram *h;
    FILE  *fp;

    fp=fopen("dork.qdp","w");
    h=fm_hg;
    for (i=0;i<h->nch;i++) 
      fprintf(fp,"%g %g\n",
	      h->limit[0]+(i+0.5)*(h->limit[1]-h->limit[0])/(float)h->nch,
	      h->values[i]);
    fprintf(fp,"no no\n");
    h=to_hg;
    for (i=0;i<h->nch;i++) 
      fprintf(fp,"%g %g\n",
	      h->limit[0]+(i+0.5)*(h->limit[1]-h->limit[0])/(float)h->nch,
	      h->values[i]);
    fclose(fp);
  }
}

/* The following routine does all the FITS file access. Originally this used */
/* the wcssubs library. It has been modified to use cfitsio and wcslib       */

int img2exhist (char *file,aspect *asp,float exhw,hgram *hg) {

  long nxpix, nypix;

  /* wcslib variables */
  int nwcs = 0;
  struct wcsprm *wcss;
  int nreject = 0;
#define NCOORD 1
#define NELEM 2
  double pixcrd[NELEM], imgcrd[NELEM], world[NELEM];
  double phi[NCOORD], theta[NCOORD];
  int stat[NCOORD];

  /* cfitsio variables */
  fitsfile *fptr = 0;
  char *fitshdr;
  float *fitsimg;
  int status = 0;
  int lhead = 0;
  int nfound, anynull;
  long naxes[2], npixels, fpixel=1;
  float nullval=0.;

  /* open the image fits file */

  fits_open_file(&fptr, file, READONLY, &status);
  if (status) {
    exsource_complain("Failed to open image file!?");
    return 1;
  }

  /* convert the header into a string that wcslib can use */

  fits_hdr2str(fptr, 1, NULL, 0, &fitshdr, &lhead, &status);
  if (status) {
    exsource_complain("Failed to read FITS image header!?");
    return 1;
  }

  /* load the wcsprm structure */

  status = wcspih(fitshdr, lhead, 0, 2, &nreject, &nwcs, &wcss);
  if (status) {
    exsource_complain("Failed to load wcsprm structure!?");
    return 1;
  }

  /* and initialize it for use */

  status = wcsset(wcss);
  if (status) {
    exsource_complain("Failed to initialize wcsprm structure!?");
    return 1;
  }

  /* read the NAXIS1 and NAXIS2 keywords to get the image size */
  fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status);
  if (status) {
    exsource_complain("Failed to read image NAXIS keywords!?");
    return 1;
  }
  nxpix = naxes[0];
  nypix = naxes[1];
  npixels = nxpix * nypix;

  /* allocate memory for the FITS image */
  fitsimg = (float*)malloc(npixels*sizeof(float));
  
  /* read the FITS image array */

  fits_read_img(fptr, TFLOAT, fpixel, npixels, &nullval, fitsimg, &anynull, &status);
  if (status) {
    exsource_complain("Failed to read FITS image!?");
    return 1;
  }

  {
    rotmat rm1,rm2,tm;
    vec    v,xv,yv,zv,xv1,yv1,zv1;
    double x,y;
    float nphot;
    
    matrix321(&rm1,M_PI,0,0);
    matrix123(&rm2,
	      asp->c.ra*M_PI/180.0,
	      asp->c.dec*M_PI/180.0,
	      asp->roll*M_PI/180.0);

    /*    fprintf(stderr,"bs matrix det: %lf\n",matrix_determinant(&rm1)); */
    /*    fprintf(stderr,"at matrix det: %lf\n",matrix_determinant(&rm2)); */
    matrix_matrix_multiply(&tm,&rm2,&rm1);

    /*    fprintf(stderr,"tm matrix det: %lf\n",matrix_determinant(&tm)); */

    xv.comp[0]=1;  xv.comp[1]=0;  xv.comp[2]=0;
    yv.comp[0]=0;  yv.comp[1]=1;  yv.comp[2]=0;
    zv.comp[0]=0;  zv.comp[1]=0;  zv.comp[2]=1;

    matrix_vector_multiply(&xv1,&tm,&xv);
    matrix_vector_multiply(&yv1,&tm,&yv);
    matrix_vector_multiply(&zv1,&tm,&zv);

    /* print_vec(&xv1); */
    /* print_vec(&yv1); */
    /* print_vec(&zv1); */
    
    {
    /* determine appropriate array size and limits for the fullres histogram. */
      double ycomp,zcomp,ycmin,ycmax,zcmin,zcmax;

      zcmin=ycmin=1e8;
      zcmax=ycmax=-1e8;
      for (y=0;y<=nypix;y+=nypix) {
	pixcrd[1] = y;
	for (x=0;x<=nxpix;x+=nxpix) {
	  pixcrd[0] = x;
	  wcsp2s(wcss,NCOORD,NELEM,pixcrd,imgcrd,phi,theta,world,stat);
	  ccoords2vec(world[0]*M_PI/180.0,world[1]*M_PI/180.0,&v);
	  /* negative z component corresponds to positive delta alpha.  */
	  ycomp=asin(dotproduct(&yv1,&v));
	  zcomp=asin(dotproduct(&zv1,&v));
	  ycmin=min(ycmin,ycomp);	  ycmax=max(ycmax,ycomp);
	  zcmin=min(zcmin,zcomp);	  zcmax=max(zcmax,zcomp);
	}
      }
      hg->limit[0]=zcmin;
      hg->limit[1]=zcmax;
      hg->nch=floor(sqrt(pow(nxpix,2.0)+pow(nypix,2.0)));
      if (hg->values != NULL) free(hg->values);
      hg->values=(float*)malloc(hg->nch*sizeof(float));
      {
	int i;
	i=hg->nch;
	while (i--) hg->values[i]=0.0;
      }
    }

    for (y=0;y<nypix;y++) {
      for (x=0;x<nxpix;x++) {
	nphot=(float)floor(fitsimg[(int)(x + y * nxpix)]);
	
	if (nphot>0) {
	  int chn,xi,yi;
	  double zc[4],yc[4],zcmin,zcmax,ycmin,ycmax; 
	  double chmin,chmax;

	  /* compute the channel number from the z component of the sky pixel.*/

	  for (xi=0;xi<2;xi++) {
	    pixcrd[0] = x + xi;
	    for (yi=0;yi<2;yi++) {
	      pixcrd[1] = y + yi;
	      wcsp2s(wcss,NCOORD,NELEM,pixcrd,imgcrd,phi,theta,world,stat);
	      ccoords2vec(world[0]*M_PI/180.0,world[1]*M_PI/180.0,&v);
	      /* negative z component corresponds to positive delta alpha.  */
	      yc[yi*2+xi]=asin(dotproduct(&yv1,&v));
	      zc[yi*2+xi]=asin(dotproduct(&zv1,&v));
	    }
	  }


	  ycmin=zcmin=1e8;
	  ycmax=zcmax=-1e8;

	  for (xi=0;xi<2;xi++) {
	    for (yi=0;yi<2;yi++) {
	      ycmin=min(ycmin,yc[yi*2+xi]);
	      ycmax=max(ycmax,yc[yi*2+xi]);
	      zcmin=min(zcmin,zc[yi*2+xi]);
	      zcmax=max(zcmax,zc[yi*2+xi]);
	    }
	  }
	  
	  if (fabs((ycmin+ycmax)/2.0) < exhw) {
	    chmin=(zcmin-hg->limit[0])/(hg->limit[1]-hg->limit[0])*hg->nch;
	    chmax=(zcmax-hg->limit[0])/(hg->limit[1]-hg->limit[0])*hg->nch;
	    
	    for (chn=max(0,floor(chmin));chn<min(hg->nch,ceil(chmax));chn++) {
	      /* increment the channel */
	      hg->values[chn] += 
		nphot * (min(chmax,chn+1)-max(chmin,chn))/fabs(chmax-chmin);
	    }
	  }
	}
      }
    }
    /* now normalize the histogram. */
    {
      double total=0.0;
      int i;
      for (i=0;i<hg->nch;i++) total+=hg->values[i];
      for (i=0;i<hg->nch;i++) hg->values[i] /= total;
    }
  }
  free(fitsimg);
  free(fitshdr);
  return 0;
}

double dotproduct (vec *v1,vec *v2) {
  return(v1->comp[0]*v2->comp[0]+
	 v1->comp[1]*v2->comp[1]+
	 v1->comp[2]*v2->comp[2]);
}

double matrix_determinant(rotmat *m) {
  double det;
  det=(m->elem[0][0]*(m->elem[1][1]*m->elem[2][2]-m->elem[1][2]*m->elem[2][1])+
       m->elem[0][1]*(m->elem[1][2]*m->elem[2][0]-m->elem[1][0]*m->elem[2][2])+
       m->elem[0][2]*(m->elem[1][0]*m->elem[2][1]-m->elem[1][1]*m->elem[2][0]))
    ;
  return(det);
}

void print_vec(vec *v) {
  fprintf(stderr,"(x,y,z)=(%f,%f,%f)\n",
	  v->comp[0],v->comp[1],v->comp[2]);
}

vec *ccoords2vec ( double ra, double dec, vec *v ) {
  double tmp;
  tmp=cos(dec);
  v->comp[0]=tmp*cos(ra);
  v->comp[1]=tmp*sin(ra);
  v->comp[2]=sin(dec);
  return(v);
}

vec *matrix_vector_multiply ( vec *prod, rotmat *m, vec *v ) {
  int i,j;

  for (i=0;i<3;i++) {
    prod->comp[i]=0.0;
    for (j=0;j<3;j++) prod->comp[i] += m->elem[i][j] * v->comp[j];
  }
  return(prod);
}

rotmat *matrix_matrix_multiply ( rotmat *prod, rotmat *m1, rotmat *m2 ) {
  int i,j,k;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++) {
      prod->elem[i][j]=0.0;
      for (k=0;k<3;k++) {
	prod->elem[i][j] += m1->elem[i][k] * m2->elem[k][j];
      }
    }
  return(prod);
}

rotmat *matrix321 ( rotmat *m, double psi, double theta, double phi ) {
  rotmat psmat,thmat,phmat,tmpmat;
  double tmp;
  int i,j;

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      psmat.elem[i][j]=phmat.elem[i][j]=thmat.elem[i][j]=0.0;
    }
    psmat.elem[i][i]=phmat.elem[i][i]=thmat.elem[i][i]=1.0;
  }

  tmp=cos(psi);   psmat.elem[0][0]=tmp;  psmat.elem[1][1]=tmp;
  tmp=sin(psi);   psmat.elem[0][1]=-tmp; psmat.elem[1][0]=tmp;

  tmp=cos(theta); psmat.elem[0][0]=tmp;  psmat.elem[2][2]=tmp;
  tmp=sin(theta); psmat.elem[0][2]=tmp;  psmat.elem[2][0]=-tmp;

  tmp=cos(phi);   psmat.elem[1][1]=tmp;  psmat.elem[2][2]=tmp;
  tmp=sin(phi);   psmat.elem[1][2]=-tmp; psmat.elem[2][1]=tmp;

  matrix_matrix_multiply(&tmpmat,&thmat,&phmat);
  matrix_matrix_multiply(m,&psmat,&tmpmat);
  return(m);
}

rotmat *matrix123 ( rotmat *m, double ra, double dec, double roll ) {
  rotmat ramat,decmat,rollmat,tmpmat;
  double tmp;
  int i,j;

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      ramat.elem[i][j]=decmat.elem[i][j]=rollmat.elem[i][j]=0.0;
    }
    ramat.elem[i][i]=decmat.elem[i][i]=rollmat.elem[i][i]=1.0;
  }

  tmp=cos(ra+M_PI);   ramat.elem[0][0]=tmp;  ramat.elem[1][1]=tmp;
  tmp=sin(ra+M_PI);   ramat.elem[0][1]=-tmp; ramat.elem[1][0]=tmp;

  tmp=cos(dec); decmat.elem[0][0]=tmp;  decmat.elem[2][2]=tmp;
  tmp=sin(dec); decmat.elem[0][2]=tmp;  decmat.elem[2][0]=-tmp;

  tmp=cos(roll);   rollmat.elem[1][1]=tmp;  rollmat.elem[2][2]=tmp;
  tmp=sin(roll);   rollmat.elem[1][2]=-tmp; rollmat.elem[2][1]=tmp;

  matrix_matrix_multiply(&tmpmat,&decmat,&rollmat);
  matrix_matrix_multiply(m,&ramat,&tmpmat);
  return(m);
}


void str2coords (char *s,aspect *a) {
  char s1[2048],s2[2048];
  double rah,ram,ras;
  double decd,decm,decs;

  rah=ram=ras=0.0;
  decd=decm=decs=0.0;
  sscanf(s,"%s %s %lf",s1,s2,&a->roll);
  sscanf(s1,"%lf:%lf:%lf",&rah,&ram,&ras);
  sscanf(s2,"%lf:%lf:%lf",&decd,&decm,&decs);
  if (decd<0.0) {
    a->c.dec=decd-(decm+decs/60)/60;
  } else {
    a->c.dec=decd+(decm+decs/60)/60;
  }
  a->c.ra=15.0*(rah+(ram+ras/60)/60);
}

void exsource_complain (char *s) {
  fprintf(stderr,"%s\n",s);
}




