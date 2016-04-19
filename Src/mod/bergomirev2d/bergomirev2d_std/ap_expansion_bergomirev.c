#include  "bergomirev2d_std.h"
#include <pnl/pnl_mathtools.h>
#include <pnl/pnl_list.h>
#include <pnl/pnl_integration.h>
#include <pnl/pnl_cdf.h>
#include <pnl/pnl_random.h>
#include <pnl/pnl_finance.h>
#include <pnl/pnl_vector_double.h>
#include <pnl/pnl_basis.h>

static int CHK_OPT(AP_EXPANSION_BERGOMIREV)(void *Opt, void *Mod)
{
  return NONACTIVE;
}
int CALC(AP_EXPANSION_BERGOMIREV)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->HelpFilenameHint = "AP_EXPANSION_BERGOMIREV";
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(AP_EXPANSION_BERGOMIREV)=
{
  "AP_EXPANSION_BERGOMIREV",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_EXPANSION_BERGOMIREV),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_EXPANSION_BERGOMIREV),
  CHK_ok,
  MET(Init)
};

