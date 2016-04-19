#include <stdlib.h>
#include <math.h>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_fft.h"
#include "math/wienerhopf.h"
#include "cgmy1d_pad.h"

#include "pnl/pnl_cdf.h"
#include"pnl/pnl_random.h"
#include"pnl/pnl_specfun.h"

static int CHK_OPT(AP_WH_FixedLookback)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_WH_FixedLookback)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  static int first=1;

  if (first)
    {
      Met->HelpFilenameHint = "AP_FASTWH_cgmy_LookbackFixed";
      Met->Par[0].Val.V_PDOUBLE=1.0;
      Met->Par[1].Val.V_PDOUBLE=0.001;
        
      first=0;
    }
  return OK;
}

PricingMethod MET(AP_WH_FixedLookback)=
{
  "AP_FastWH",
  { {"Scale of logprice range", DOUBLE, {100}, ALLOW},
    {"Space Discretization Step",DOUBLE,{500},ALLOW},
    {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_WH_FixedLookback),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_WH_FixedLookback),
  CHK_split,
  MET(Init)
};

