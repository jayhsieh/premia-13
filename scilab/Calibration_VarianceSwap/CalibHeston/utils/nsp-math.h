#ifndef NSP_INC_MATH 
#define NSP_INC_MATH 

/*
 * This Software is GPL (Copyright ENPC 1998-2005) 
 * Jean-Philippe Chancelier Enpc/Cermics         
 */

#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "nsp-machine.h"

#ifndef TRUE
#define TRUE (1)
#endif 

#ifndef FALSE 
#define FALSE (0)
#endif 

#define OK 0
#define FAIL 1

/* sur sun solaris finite est dans ieeefp **/

#if (defined(sun) && defined(SYSV)) 
#include <ieeefp.h>
#endif 

#ifndef HAVE_FINITE
#if defined(WIN32) && !(defined __CYGWIN32__) && !(defined __ABSC__)
#include <float.h>
#define finite(x) _finite(x) 
#else 
/* should do something better here */
#define finite(x) FINITE_IS_UNDEFINED
#endif
#endif 

#ifndef HAVE_ISNAN 
#if defined(WIN32) && !(defined __CYGWIN32__) && !(defined __ABSC__)
#include <float.h>
#define isnan(x) _isnan(x)
#else 
/* should do something better here */
#define isnan(x)  ISNAN_IS_UNDEFINED
#endif 
#endif 

/* just a synonym */

#define ISNAN(x) isnan(x)

#ifndef Abs 
#define Abs(x) ( ( (x) >= 0) ? (x) : -( x) )
#endif

#ifndef Min
#define Min(x,y)	(((x)<(y))?(x):(y))
#endif 

#ifndef Max 
#define Max(x,y)	(((x)>(y))?(x):(y))
#endif 

extern double Mini();  /* XXXX a mettre ailleurs **/
extern double Maxi();  /* XXXX a mettre ailleurs **/

#define PI0 (int *) 0
#define PD0 (double *) 0
#define SMDOUBLE 1.e-200 /* Smalest number to avoid dividing by zero */

#define linint(x) ((int) floor(x + 0.5 )) 
#define inint(x) ((int) floor(x + 0.5 ))  

/* a revoir precisement un jour ou l''autre XXXXXX **/
/* nearest : **/
#define anint(x) rint(x) 
/* partie entiere **/
#define aint(x) ((x>= 0 ) ? floor(x)  : ceil(x))
/* from fortran but no pointer here */
#define d_nint(x) (x)>=0 ? floor(x + .5) : -floor(.5 - x)


/* missing prototype */
double tgamma(double);

/* Les arguments des fonction XWindows doivent etre des int16 ou unsigned16 */

#define int16max   0x7FFF
#define uns16max   0xFFFF

#ifdef lint5
#include <sys/stdtypes.h>
#define MALLOC(x) malloc(((size_t) x))
#define FREE(x) {if (x  != NULL) { free((void *) x); x= NULL;};}
#define REALLOC(x,y) realloc((void *) x,(size_t) y)
#else
#define MALLOC(x) malloc(((unsigned) x))
#define FREE(x) {if (x  != NULL) {free((char *) x); x = NULL;}}
#define REALLOC(x,y) realloc((char *) x,(unsigned) y)
#endif

/* void name for object **/
#define NVOID ""

#if defined(WIN32) || defined(__STDC__)
# ifndef HAS_STDARG
#  define HAS_STDARG
# endif
#else 
#endif 



/* M_PI and M_E are defined in math.h */

#ifndef M_PI
#define M_PI    3.14159265358979323846 
#endif

#ifndef M_E
#define M_E     2.7182818284590452354
#endif 

#ifndef M_LOG10E
#define M_LOG10E 0.43429448190325182765
#endif

#ifdef WIN32 
extern double acosh(double);
extern double asinh(double);
extern double atanh(double);
#endif 

/*
 * 
 */

extern  double nsp_dlamch (char *cmach);

#endif /* NSP_MATH */



