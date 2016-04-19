#ifndef  _DUP1D_STD_H
#define _DUP1D_STD_H

#include "dup1d/dup1d.h"
#include "std/std.h"

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "numfunc.h"
#include "transopt.h"
#include "math/linsys.h"
#include <float.h>

/* Local Volatility Examples Sigma(t,x) */

extern double MOD_OPT(lib_volatility)(double t, double x,int sigma_type);
extern double MOD_OPT(lib_volatility_x)(double t, double x,int sigma_type);
extern double MOD_OPT(lib_volatility_xx)(double t, double x,int sigma_type);

#define volatility MOD_OPT(lib_volatility)
#define volatility_x MOD_OPT(lib_volatility_x)
#define volatility_xx MOD_OPT(lib_volatility_xx)

#endif
