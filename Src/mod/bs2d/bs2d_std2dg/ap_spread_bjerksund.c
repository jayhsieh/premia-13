#include  "bs2d_std2dg.h"
#include "pnl/pnl_cdf.h"


static int CHK_OPT(AP_Spread_Bjerksund)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_Spread_Bjerksund)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(AP_Spread_Bjerksund)=
{
  "AP_Spread_Bjerksund",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_Spread_Bjerksund),
  {{"Price",DOUBLE,{100},FORBID},{"Delta1",DOUBLE,{100},FORBID} ,{"Delta2",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_Spread_Bjerksund),
  CHK_ok,
  MET(Init)
} ;
