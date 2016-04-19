#ifndef  _BS1D_STD_H
#define _BS1D_STD_H

#include "bs1d/bs1d.h"
#include "std/std.h"

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_random.h"
#include "numfunc.h"
#include "pnl/pnl_cdf.h"
#include "transopt.h"
#include "math/linsys.h"
#include <float.h>

#ifdef USE_ND1
static double Nd1(double s,double r,double divid,double sigma,double T,double K) 
{
  double d1=(log(s/K)+(r-divid+0.5*sigma*sigma)*T)/(sigma*sqrt(T));
  return  cdf_nor(d1); 
}
#endif

#endif
