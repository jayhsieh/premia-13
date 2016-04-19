#include "cev1d_std.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_root.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_root.h"

static int CHK_OPT(AP_BGM_Cev)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_BGM_Cev)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0) Met->init=1;
  return OK;
}

PricingMethod MET(AP_BGM_Cev)=
{
  "AP_BGM_Cev",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_BGM_Cev),
  {{"At the Money Price",DOUBLE,{100},FORBID},{"At the Strike Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_BGM_Cev),
  CHK_ok,
  MET(Init)
};

