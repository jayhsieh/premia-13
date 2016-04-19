#include  "bergomi2d_std.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_list.h"
#include "pnl/pnl_integration.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_vector_double.h"
#include "pnl/pnl_basis.h"

static int CHK_OPT(AP_EXPANSION_OA)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_EXPANSION_OA)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->HelpFilenameHint = "AP_EXPANSION_Bergomi";
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(AP_EXPANSION_OA)=
{
  "AP_EXPANSION_OULDALY",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_EXPANSION_OA),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_EXPANSION_OA),
  CHK_ok,
  MET(Init)
};

