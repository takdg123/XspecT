#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <fitsio.h>
#include "xslimbb.h"
#ifndef DEBUG
#include "XSFunctions/Utilities/xsFortran.h"
#include "cfortran.h"
#endif


#define XWRITE(A,B) CCALLSFSUB2(XWRITE,xwrite,STRING,INT,A,B)  
// const Real* energy
// energy array (size Nflux+1)
//
// int Nflux
// size of flux array
//
// const Real* parameter
// parameter values
//
// int spectrum
// spectrum number of model component being calculated
//
// Real* flux
// output flux array
//
// Real* fluxError
// output flux error array (optional)
//
// const char* init
// initialization string (see below)

//typedef Real double;

#define erg2kev         6.241507e+08         // erg->keV (erg/(1e3*electronvolt))

#define DEFAULT_SLIMBB_TABLE   "slimbb-full.fits"

#define CHECKERR(cond,message) if(!(cond)){sprintf(errorstr,message); goto error;}
#define ave(a, b, w) ((1.0-(w))*(a) + (w)*(b))
#define logave(a, b, w) (exp((1.0-(w))*log(a) + (w)*log(b)))

#ifdef DEBUG
    #define FGMSTR(x) "\0"
#endif


void indexxx(float* array, long n, float value, long* index, double* w, int* status)
{
    *index=0; 
    *w=0.0;
    if ((*status)!=0) return;  // return on error from previous call

    if ((n<=0) || (value<array[0]) || (value>array[n-1])) {*status=1; return;}   // return error if out of range
    if (n==1) return;  // return OK and index=0,w=0.0 if only one element in the array
 
    int i = 0;
    while ((i<n-2) && (value > array[i+1])) i++;
    *w = (value-array[i])/(array[i+1]-array[i]);
    *index = i;
}



void slimbb_log(FILE* logfile, const char *template, ...)
//**************************************************
{
    if (logfile) {
        va_list ap;
        va_start (ap, template);
        vfprintf(logfile, template, ap);
        va_end (ap);
    }
}





//extern “C” 
void slimbbmodel(const double* energy, int Nflux, const double* parameter, int spectrum, double* flux, double* fluxError, const char* init)
{
// energy ... array of boundaries of observed energy bins [keV]; dimensions [0:Nflux]
// Nflux  ... size of flux array
// flux   ... output array of observed energy bins [photons/cm^2/s]; dimensions [0:Nflux-1]
// param  ... parameters of the model
// spectrum  ... no meaning
// init      ... no meaning

    char errorstr[256];    
    fitsfile *fptr;
    int status;
    int i;
    float *spin_grid=NULL, *lumi_grid=NULL, *incl_grid=NULL, *ener_grid=NULL, *alph_grid=NULL;
    float *f=NULL;

    // open log file if in debug mode
    int debug = (atof(FGMSTR("SLIMBB_DEBUG"))==1);
    #ifdef DEBUG
        debug = 1;
    #endif
    FILE *logfile = (debug) ? fopen("slimbb.log", "w") : NULL;

    // read model parameters
    double param_mass =  parameter[0];
    double param_spin =  parameter[1];
    double param_lumi =  parameter[2];
    double param_alph =  parameter[3];
    double param_incl =  parameter[4]/180.*M_PI;
    double param_dist =  parameter[5];
    double param_hard =  parameter[6];
    int    param_limb = (parameter[7] > 0.0);
    int    param_vert = (parameter[8] > 0.0);

    // read fits table location and name
    //char* tmpstr;
    char slimbb_path[256];
    strncpy(slimbb_path, FGMODF(), 255);
    slimbb_path[255] = 0;
    char slimbb_table[256] = DEFAULT_SLIMBB_TABLE;
    char slimbb_table_full_path[512];
    if (strlen(FGMSTR("SLIMBB_DIR"))  >0) strncpy(slimbb_path, FGMSTR("SLIMBB_DIR"), 255);
    if (strlen(FGMSTR("SLIMBB_TABLE"))>0) {
       strncpy(slimbb_table, FGMSTR("SLIMBB_TABLE"), 255);
       slimbb_table[255]=0;
    }
    // compile full path to fits table
    sprintf(slimbb_table_full_path, "%s/%s", slimbb_path, slimbb_table);

    // which flux column to use
    int flux_column;
    if (param_hard > 0.0)
        flux_column = (param_limb) ? (param_vert?FLUX_COL_FLUX_LV:FLUX_COL_FLUX_L0) : (param_vert?FLUX_COL_FLUX_IV:FLUX_COL_FLUX_I0);
    else
        flux_column = param_vert?FLUX_COL_FLUX_SV:FLUX_COL_FLUX_S0;

    // log input parameters
    slimbb_log(logfile, "--------------\n");
    slimbb_log(logfile, "params: mass = %.3f\n", param_mass);
    slimbb_log(logfile, "params: spin = %.5f\n", param_spin);
    slimbb_log(logfile, "params: lumi = %.5f\n", param_lumi);
    slimbb_log(logfile, "params: alph = %.3f\n", param_alph);
    slimbb_log(logfile, "params: incl = %.3f (%.5f rad)\n", param_incl/M_PI*180., param_incl);
    slimbb_log(logfile, "params: dist = %.3f\n", param_dist);
    slimbb_log(logfile, "params: hard = %.3f\n", param_hard);
    slimbb_log(logfile, "params: limb = %d\n", param_limb);
    slimbb_log(logfile, "params: vert = %d\n", param_vert);
    slimbb_log(logfile, "params: table path = %s\n", slimbb_table_full_path);
    slimbb_log(logfile, "params: flux column = %d\n", flux_column);


    // open data table 
    status = 0;
    fits_open_table(&fptr, slimbb_table_full_path, READONLY, &status);
    if (status!=0) fptr = NULL;
    CHECKERR(status==0, "failed to open data file");

    // reading META table
    long nspin, nlumi, nincl, nener, nalph;
    float ref_mass, ref_dist;
    fits_movnam_hdu(fptr, BINARY_TBL, "META", 0, &status);
    CHECKERR(status==0, "failed to find META extension");
    {
        // read dimensions of grid arrays
        fits_read_col(fptr, TLONG, META_COL_N, META_ROW_INCL_GRID, 1, 1, NULL, &nincl, NULL, &status);
        fits_read_col(fptr, TLONG, META_COL_N, META_ROW_LUMI_GRID, 1, 1, NULL, &nlumi, NULL, &status);
        fits_read_col(fptr, TLONG, META_COL_N, META_ROW_SPIN_GRID, 1, 1, NULL, &nspin, NULL, &status);
        fits_read_col(fptr, TLONG, META_COL_N, META_ROW_ENER_GRID, 1, 1, NULL, &nener, NULL, &status);
        fits_read_col(fptr, TLONG, META_COL_N, META_ROW_ALPH_GRID, 1, 1, NULL, &nalph, NULL, &status);
        CHECKERR(status==0, "failed to read grid dimensions");

        // allocation of grids
        spin_grid = (float*)calloc(nspin, sizeof(float));
        lumi_grid = (float*)calloc(nlumi, sizeof(float));
        incl_grid = (float*)calloc(nincl, sizeof(float));
        ener_grid = (float*)calloc(nener, sizeof(float));
        alph_grid = (float*)calloc(nalph, sizeof(float));
    
        // read grid values
        fits_read_col(fptr, TFLOAT, META_COL_VALUE, META_ROW_SPIN_GRID, 1, nspin, NULL, spin_grid, NULL, &status);
        fits_read_col(fptr, TFLOAT, META_COL_VALUE, META_ROW_LUMI_GRID, 1, nlumi, NULL, lumi_grid, NULL, &status);
        fits_read_col(fptr, TFLOAT, META_COL_VALUE, META_ROW_INCL_GRID, 1, nincl, NULL, incl_grid, NULL, &status);
        fits_read_col(fptr, TFLOAT, META_COL_VALUE, META_ROW_ENER_GRID, 1, nener, NULL, ener_grid, NULL, &status);
        fits_read_col(fptr, TFLOAT, META_COL_VALUE, META_ROW_ALPH_GRID, 1, nalph, NULL, alph_grid, NULL, &status);
        fits_read_col(fptr, TFLOAT, META_COL_VALUE, META_ROW_REF_MASS,  1, 1, NULL, &ref_mass, NULL, &status);
        fits_read_col(fptr, TFLOAT, META_COL_VALUE, META_ROW_REF_DIST,  1, 1, NULL, &ref_dist, NULL, &status);
        CHECKERR(status==0, "failed to read grid values");

        if (debug) {
            slimbb_log(logfile, "--------------\n");
            slimbb_log(logfile, "meta: nspin = %ld\n", nspin);
            slimbb_log(logfile, "meta: nlumi = %ld\n", nlumi);
            slimbb_log(logfile, "meta: nincl = %ld\n", nincl);
            slimbb_log(logfile, "meta: nener = %ld\n", nener);
            slimbb_log(logfile, "meta: nalph = %ld\n", nalph);
            slimbb_log(logfile, "meta: spin_grid ="); for (i=0;i<nspin;i++) slimbb_log(logfile, " %.3f", spin_grid[i]); slimbb_log(logfile, "\n");
            slimbb_log(logfile, "meta: lumi_grid ="); for (i=0;i<nlumi;i++) slimbb_log(logfile, " %.3f", lumi_grid[i]); slimbb_log(logfile, "\n");
            slimbb_log(logfile, "meta: incl_grid ="); for (i=0;i<nincl;i++) slimbb_log(logfile, " %.3f", incl_grid[i]); slimbb_log(logfile, "\n");
            slimbb_log(logfile, "meta: ener_grid ="); for (i=0;i<nener;i++) slimbb_log(logfile, " %.3e", ener_grid[i]); slimbb_log(logfile, "\n");
            slimbb_log(logfile, "meta: alph_grid ="); for (i=0;i<nalph;i++) slimbb_log(logfile, " %.3f", alph_grid[i]); slimbb_log(logfile, "\n");
            slimbb_log(logfile, "meta: ref_mass = %.3e\n", ref_mass);
            slimbb_log(logfile, "meta: ref_dist = %.3e\n", ref_dist);
        }
    }

    // calculate grid indices 
    long ialph, ispin, ilumi, iincl;
    double walph, wspin, wlumi, wincl;
    indexxx(alph_grid, nalph, param_alph, &ialph, &walph, &status);
    indexxx(spin_grid, nspin, param_spin, &ispin, &wspin, &status);
    indexxx(lumi_grid, nlumi, param_lumi, &ilumi, &wlumi, &status);
    indexxx(incl_grid, nincl, param_incl, &iincl, &wincl, &status);
    slimbb_log(logfile, "--------------\n");
    slimbb_log(logfile, "grid indices: alph = %ld (%ld) / w=%.3f\n", ialph, nalph-1, walph);
    slimbb_log(logfile, "grid indices: spin = %ld (%ld) / w=%.3f\n", ispin, nspin-1, wspin);
    slimbb_log(logfile, "grid indices: lumi = %ld (%ld) / w=%.3f\n", ilumi, nlumi-1, wlumi);
    slimbb_log(logfile, "grid indices: incl = %ld (%ld) / w=%.3f\n", iincl, nincl-1, wincl);
    CHECKERR(status==0, "failed to constrain grid range");


    // read flux tables ...
    int ia;
    float* alpha_flux[2];
    alpha_flux[0] = calloc(nener, sizeof(float));
    alpha_flux[1] = calloc(nener, sizeof(float));
    for (ia=0; ia<=1; ia++) {
        char table_name[16];
        sprintf(table_name, "FLUX-%.3f", (nalph>1) ? alph_grid[ialph+ia] : alph_grid[0]);
        slimbb_log(logfile, "--------------\n");
        slimbb_log(logfile, "flux table: reading %s\n", table_name);
        fits_movnam_hdu(fptr, BINARY_TBL, table_name, 0, &status);
        CHECKERR(status==0, "failed to open model table for alpha");

        float *f000=NULL, *f010=NULL, *f001=NULL, *f011=NULL, *f100=NULL, *f110=NULL, *f101=NULL, *f111=NULL;
        long rowindex000, rowindex100, rowindex010, rowindex001, rowindex101, rowindex011, rowindex110, rowindex111;
        // read a box of spectra from fits table: fxxx = f(spin)(lumi)(incl)
        f000 = (float*)calloc(nener, sizeof(float));
        f100 = (float*)calloc(nener, sizeof(float));
        f010 = (float*)calloc(nener, sizeof(float));
        f001 = (float*)calloc(nener, sizeof(float));
        f101 = (float*)calloc(nener, sizeof(float));
        f011 = (float*)calloc(nener, sizeof(float));
        f110 = (float*)calloc(nener, sizeof(float));
        f111 = (float*)calloc(nener, sizeof(float));
        rowindex000 = (ispin+0)*nlumi*nincl + (ilumi+0)*nincl + (iincl+0);
        rowindex100 = (ispin+1)*nlumi*nincl + (ilumi+0)*nincl + (iincl+0);
        rowindex010 = (ispin+0)*nlumi*nincl + (ilumi+1)*nincl + (iincl+0);
        rowindex001 = (ispin+0)*nlumi*nincl + (ilumi+0)*nincl + (iincl+(int)(nincl>1)); // special case when inclination grid only contains one element
        rowindex101 = (ispin+1)*nlumi*nincl + (ilumi+0)*nincl + (iincl+(int)(nincl>1));
        rowindex011 = (ispin+0)*nlumi*nincl + (ilumi+1)*nincl + (iincl+(int)(nincl>1));
        rowindex110 = (ispin+1)*nlumi*nincl + (ilumi+1)*nincl + (iincl+0);
        rowindex111 = (ispin+1)*nlumi*nincl + (ilumi+1)*nincl + (iincl+(int)(nincl>1));
        slimbb_log(logfile, "flux table: rowindex = %ld %ld %ld %ld %ld %ld %ld %ld\n", rowindex000, rowindex100, rowindex010, rowindex001, rowindex101, rowindex011, rowindex110, rowindex111);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex000+1, 1, nener, NULL, f000, NULL, &status);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex100+1, 1, nener, NULL, f100, NULL, &status);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex010+1, 1, nener, NULL, f010, NULL, &status);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex001+1, 1, nener, NULL, f001, NULL, &status);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex101+1, 1, nener, NULL, f101, NULL, &status);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex011+1, 1, nener, NULL, f011, NULL, &status);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex110+1, 1, nener, NULL, f110, NULL, &status);
        fits_read_col(fptr, TFLOAT, flux_column, rowindex111+1, 1, nener, NULL, f111, NULL, &status);
        CHECKERR(status==0, "failed to read flux array");

        // cosine interpolation
        //wspin = 0.5*(1.-cos(M_PI*wspin));
        //wlumi = 0.5*(1.-cos(M_PI*wlumi));
        //wincl = 0.5*(1.-cos(M_PI*wincl));

        // calculate spectrum using trilinear interpolation of log values
        // result in units [log(#photons cm^-2 s^-1 keV^-1)]
        for (i=0; i<nener; i++)    alpha_flux[ia][i] = 
            log(f000[i]) * (1.0-wspin) * (1.0-wlumi) * (1.0-wincl) +
            log(f100[i]) *      wspin  * (1.0-wlumi) * (1.0-wincl) +
            log(f010[i]) * (1.0-wspin) *      wlumi  * (1.0-wincl) +
            log(f001[i]) * (1.0-wspin) * (1.0-wlumi) *      wincl  +
            log(f101[i]) *      wspin  * (1.0-wlumi) *      wincl  +
            log(f011[i]) * (1.0-wspin) *      wlumi  *      wincl  +
            log(f110[i]) *      wspin  *      wlumi  * (1.0-wincl) +
            log(f111[i]) *      wspin  *      wlumi  *      wincl  ;

        // free arrays
        free(f000);
        free(f010);
        free(f001);
        free(f011);
        free(f100);
        free(f110);
        free(f101);
        free(f111);
    }

    // interpolate spectrum in alpha (linear) and unlog result
    f = (float*)calloc(nener, sizeof(float));
    for (i=0; i<nener; i++) {
        f[i]  = exp(ave(alpha_flux[0][i], alpha_flux[1][i], walph));
        f[i] *= erg2kev/ener_grid[i];   // convert erg -> #photons
    }
    free(alpha_flux[0]);
    free(alpha_flux[1]);

    // apply transformations to the spectrum
    slimbb_log(logfile, "--------------\n");
    slimbb_log(logfile, "transformations: by mass = %.4f\n", param_mass/ref_mass);
    slimbb_log(logfile, "transformations: by dist = %.4f\n", param_dist/(ref_dist/1e3));
    slimbb_log(logfile, "transformations: by hard = %.4f\n", param_hard);
    for (i=0; i<nener; i++) {
        f[i] *= param_mass/ref_mass; ener_grid[i] *= pow(param_mass/ref_mass,-.5);                      // accounts for mass
        f[i] *= pow(param_mass/ref_mass,.5); ener_grid[i] *= pow(param_mass/ref_mass,.25);              // accounts for luminosity (which is \propto mass)
        f[i] *= 1./pow(param_dist/(ref_dist/1e3),2.);                                                   // accounts for distance
        if (param_hard>0.0) { f[i] *= 1./pow(fabs(param_hard),2.); ener_grid[i] *= fabs(param_hard); }  // accounts for hardening
    }

    //for (i=0; i<nener-1; i+=nener/200) printf("%e  %e  %e\n", ener_grid[i], ener_grid[i+1], logave(f[i],f[i+1],0.5));


    // project spectrum to the supplied energy grid
    // output in [photons/cm^2/s]
    slimbb_log(logfile, "--------------\n");
    slimbb_log(logfile, "spectrum: en1, en2, flux1, flux2, integrated flux\n");
    for (i=0; i<Nflux; i++) {
        long iflux;
        double wflux;
        status = 0;
        indexxx(ener_grid, nener, energy[i], &iflux, &wflux, &status);
        double f1 = (status==0) ? logave(f[iflux],f[iflux+1],wflux) : ((energy[i]<ener_grid[0])?f[0]*pow(energy[i]/ener_grid[0], -2./3.):1e-50) ;
        if (isnan(f1)) fprintf(stderr, "NaN - f1=%e  status=%d  en=%e\n", f1, status, energy[i]);

        status = 0;
        indexxx(ener_grid, nener, energy[i+1], &iflux, &wflux, &status);
        double f2 = (status==0) ? logave(f[iflux],f[iflux+1],wflux) : ((energy[i]<ener_grid[0])?f[0]*pow(energy[i]/ener_grid[0], -2./3.):1e-50) ;
        if (isnan(f2)) fprintf(stderr, "NaN - f2=%e  status=%d  en=%e  F=%e-%e\n", f2, status, energy[i], f[iflux], f[iflux+1]);

        flux[i] = logave(f1,f2,0.5) * (energy[i+1]-energy[i]);
        if (isnan(flux[i])) flux[i] = 1e-50;
        if (logfile) fprintf(logfile, "spectrum: %e  %e    %e  %e    %e    %e\n", energy[i], energy[i+1], f1, f2, logave(f1,f2,0.5), flux[i]);
    }



    goto finish;
    
    // [exception block] report errors 
    error:
    fits_report_error(stderr, status);
    fprintf(stderr, "ERROR: %s\n", errorstr);
    slimbb_log(logfile, "ERROR: %s\n", errorstr);

    // [final block] free all dynamic memory
    finish:    

    status = 0;
    if (fptr) fits_close_file(fptr, &status);
    //CHECKERR(status==0, "failed to close data file");

    free(f);
    free(incl_grid);
    free(lumi_grid);
    free(spin_grid);
    free(ener_grid);
    free(alph_grid);

    if (logfile) fclose(logfile);
}




