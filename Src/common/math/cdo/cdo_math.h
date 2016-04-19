#ifndef __CDO_MATHS_H__
#define __CDO_MATHS_H__
#include        <stdlib.h>
#include        <math.h>
#include        <float.h>

#include "pnl/pnl_random.h"
#include "pnl/pnl_complex.h"

/** Definition of some constants (computed with GNU bc).
 */
#ifndef MINDOUBLE
#define MINDOUBLE       4.94065645841246544e-324 ///< Defined in <tt>float.h<\tt> but problem with cross-compilation for windows.
#endif
#ifndef MAXDOUBLE
#define MAXDOUBLE       1.79769313486231570e+308
#endif
double evaluate_poly(int degree, const double *a, double x);
double evaluate_dpoly(int degree, const double *a, double x);
/*function from void
 */
typedef double (fun_void_R)(void);
typedef double (pfun_void_R)(const void *p);
typedef double (fun_R_R)(double);
typedef dcomplex (fun_R_complex)(double);
typedef dcomplex (pfun_R_complex)(const double x, const void *p);
double student_cdf(double t1, double x);
double student_inv_cdf( double t1, double x);
double simulate_student(double t);
#endif 
