#ifndef __BW_H__
#define __BW_H__

#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

#undef __BEGIN_DECLS
#undef __END_DECLS

typedef double Real;

#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif


__BEGIN_DECLS

/* copied from heacore/gsl/specfunc/error.h */
#define DOMAIN_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("domain error", GSL_EDOM); } while(0)

/* copied from heacore/gsl/specfunc/eval.h */
#define EVAL_RESULT(fn) \
   gsl_sf_result result; \
   int status = fn; \
   if (status != GSL_SUCCESS) { \
     GSL_ERROR_VAL(#fn, status, result.val); \
   } ; \
   return result.val;

int hyperg_2F2_series(const double a, const double b, const double c, const double d,
                  const double x, 
                  gsl_sf_result * result
                  );
int gsl_sf_hyperg_2F2_e(double a, double b, const double c, const double d,
                       const double x,
                       gsl_sf_result * result);

double gsl_sf_hyperg_2F2(double a, double b, double c, double d, double x);

int alphaf_e(double zeta, double R_star, double M_star, gsl_sf_result * result );
double alphaf(double zeta, double R_star, double M_star);
double wf(double zeta);
double kappaf(double delta);
double lambdaf(double n);
int muf_e(double delta, double zeta, double n, gsl_sf_result *result);
int Xn_e(double zeta, double n, gsl_sf_result *result);
int An_factor_e(double zeta, double n, double m, gsl_sf_result *result);
int An_e(double zeta, double n, gsl_sf_result *result);
double sigma_parf(double r0, double Mdot, double zeta);
double T_thf(double Mdot, double r0);
double tau_thf(double Mdot, double r0, double R_star, double M_star, double zeta, double *T_th);
double z_thf(double Mdot, double r0, double R_star, double M_star, double zeta, double *T_th);
double tau_maxf(double Mdot, double r0, double R_star, double M_star, double zeta);
int ser_cyc_factor_e(double e, double n, double *params, gsl_sf_result *result);
double ser_cyc(double e, double *params);
double Hf(double chi);
int cyc_e(double e, double *params, gsl_sf_result *result);
/*double tauf(double z, double zeta, double r0);
double zf(double tau, double zeta, double r0);*/
int gn_e(double tau, int n, gsl_sf_result *result);
int ser_Green_factor_e(double tau0, double e0, double e, double n, double *params, gsl_sf_result *result);
int ser_Green_e(double tau0, double e0, double e, double *params, gsl_sf_result *result);
int Green_e(double tau0, double e0, double e, double *params, gsl_sf_result *result);

double bb_integrand(double e0, void *params);
double BB(double e, double *params);
double chi_abs_int(double tau, void *int_params);
double chi_absf(double *params);
double Bn_int(double chi0, void *int_params);
double Bn(double chi, int n);
double ser_FF_factor(double chi, double n, double *params);
double ser_FF(double chi,double *params);
double FF(double e, double *params);

void beckerwolff(const Real* energy, int Nflux, const Real* parameter, int spectrum, Real* flux);
double integralbw(const double *e_a, gsl_spline *spine, gsl_interp_accel *acc);
__END_DECLS

#endif /* __BW_H__ */
