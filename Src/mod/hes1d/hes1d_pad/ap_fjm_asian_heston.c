#include  "hes1d_pad.h"
#include <pnl/pnl_mathtools.h>
#include <pnl/pnl_root.h>

static int CHK_OPT(AP_FJM_ASIAN_HESTON)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_FJM_ASIAN_HESTON)(void *Opt, void *Mod, PricingMethod *Met)
{
return AVAILABLE_IN_FULL_PREMIA;
}


static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  //int type_generator;
  if ( Met->init == 0)
    {
      Met->init=1;

    }

  return OK;
}

PricingMethod MET(AP_FJM_ASIAN_HESTON)=
{
   "AP_FJM_ASIAN_HESTON",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_FJM_ASIAN_HESTON),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_FJM_ASIAN_HESTON),
  CHK_mc,
  MET(Init)
};
