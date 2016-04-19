#include "locvolhw1d_std.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_root.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_root.h"

static int CHK_OPT(AP_BGM_Locvolhw)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_BGM_Locvolhw)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0) Met->init=1;
  return OK;
}

PricingMethod MET(AP_BGM_Locvolhw)=
{
  "AP_BGM_Locvolhw",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_BGM_Locvolhw),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_BGM_Locvolhw),
  CHK_ok,
  MET(Init)
};

