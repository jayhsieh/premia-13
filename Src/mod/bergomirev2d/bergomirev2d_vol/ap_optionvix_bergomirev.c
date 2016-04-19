#include  "bergomirev2d_vol.h"

#include <pnl/pnl_mathtools.h>
#include <pnl/pnl_list.h>
#include <pnl/pnl_integration.h>
#include <pnl/pnl_cdf.h>
#include <pnl/pnl_random.h>
#include <pnl/pnl_finance.h>
#include <pnl/pnl_vector_double.h>
#include <pnl/pnl_basis.h>

static int CHK_OPT(AP_OPTIONVIX_BERGOMIREV)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_OPTIONVIX_BERGOMIREV)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->HelpFilenameHint = "AP_OPTIONVIX_BERGOMIREV";
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(AP_OPTIONVIX_BERGOMIREV)=
{
  "AP_OULDALY",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_OPTIONVIX_BERGOMIREV),
  {{"Price",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_OPTIONVIX_BERGOMIREV),
  CHK_ok,
  MET(Init)
};

