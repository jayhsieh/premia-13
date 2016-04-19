#ifndef  _BS1D_PAD_H
#define _BS1D_PAD_H

#include "bs1d/bs1d.h"
#include "pad/pad.h"

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_cdf.h"
#include "transopt.h"
#include "numfunc.h"
#include "pnl/pnl_complex.h"
#include "math/moments.h"

extern int MOD_OPT(Analytic_KemnaVorst_lib)(double pseudo_stock, double pseudo_strike, double time_spent, NumFunc_2  *p, 
					double t, double r, double divid, double *ptprice, double *ptdelta);

#define Analytic_KemnaVorst MOD_OPT(Analytic_KemnaVorst_lib)

#endif
