#ifndef __MATHS_H__
#define __MATHS_H__
#include 	<stdlib.h>
#include 	<math.h>
#include 	<complex.h>
#include 	<float.h>
#include 	<fftw3.h>
#include 	"mt19937.h"

/** Definition of some constants (computed with GNU bc).
 */
#define C_PI  		3.14159265358979323844	///< Constant \f$\pi\f$.
#define C_2PI  		6.28318530717958647688	///< Constant \f$2 \pi\f$.
#define C_SQRT1_2  	0.70710678118654752440	///< Constant \f$1/\sqrt(2)\f$.
#define C_SQRT1_2PI	0.39894228040143267794	///< Constant \f$1/\sqrt(2 \pi)\f$.
#ifndef MINDOUBLE
#define MINDOUBLE	4.94065645841246544e-324 ///< Defined in <tt>float.h<\tt> but problem with cross-compilation for windows.
#endif
#ifndef MAXDOUBLE
#define MAXDOUBLE	1.79769313486231570e+308
#endif
double pp(const double x,const double y);
double evaluate_poly(const int degree,const double *a,const double x);
double evaluate_dpoly(const int degree, const double *a,const double x);
/*function from void
 */
typedef double (fun_void_R)(void);
typedef double (pfun_void_R)(const void *p);
typedef double (fun_R_R)(double);
typedef complex (fun_R_complex)(double);
typedef complex (pfun_R_complex)(const double x, const void *p);
double gauss_density(const double x);
double abso(const double x);
double gaussian();
double gaussians();
double simulate_student(const double t);
double gaussian_cdf(const double x);
double gaussian_inv_cdf(const double p);
double bessel_I1(const double x);
double bessel_K1(const double x);
double max(double a,double b);
double gammas(const double xx);
double student_cdf(const double t1,const double x);
double student_inv_cdf(const double t1,const double x);
double density_chi2(const double t1,const double x);
double** inv_mat(double **a,int n);
#endif 
