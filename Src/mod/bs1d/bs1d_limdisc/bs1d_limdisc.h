#ifndef  _BS1D_LIMDISC_H
#define _BS1D_LIMDISC_H

#include "bs1d/bs1d.h"
#include "limdisc/limdisc.h"

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "transopt.h"
#include "numfunc.h"
#include "pnl/pnl_complex.h"
#include "math/moments.h"


int MOD_OPT(formula_lib)(double s,double k,double r,double divid,double sigma,double t,double l, double rebate,int phi,int eta,
	    double *A,double *B,double *C,double *D,double *E,double *F);

#define formula MOD_OPT(formula_lib)


#endif
