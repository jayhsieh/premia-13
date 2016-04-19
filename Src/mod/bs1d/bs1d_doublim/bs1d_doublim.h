
#ifndef  _BS1D_DOUBLIM_H
#define _BS1D_DOUBLIM_H

#include "bs1d/bs1d.h"
#include "doublim/doublim.h"

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_cdf.h"
#include "numfunc.h"
#include "transopt.h"
#include "pnl/pnl_complex.h"


/* Utility function shared 
 *
 *
 */ 


extern int MOD_OPT(PutOut_KunitomoIkeda_91_lib)(double s,NumFunc_1  *L,NumFunc_1  *U,NumFunc_1  *Rebate,NumFunc_1 *PayOff,
						double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);

extern int MOD_OPT(CallOut_KunitomoIkeda_91_lib)(double s,NumFunc_1 *L,NumFunc_1 *U,NumFunc_1  *Rebate,NumFunc_1 *PayOff,
						 double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);

extern double MOD_OPT(Boundary_lib)(double s,NumFunc_1*p,double t,double r,double divid,double sigma);

/* aliases */

#define Boundary MOD_OPT(Boundary_lib)
#define CallOut_KunitomoIkeda_91 MOD_OPT(CallOut_KunitomoIkeda_91_lib)
#define PutOut_KunitomoIkeda_91 MOD_OPT(PutOut_KunitomoIkeda_91_lib)


#endif

