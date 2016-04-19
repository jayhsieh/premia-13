#include <stdlib.h>
#include "bsdisdiv1d_std.h"
#include "error_msg.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"

#define MALLOC_DOUBLE(n) malloc(n * sizeof(double))

static int CHK_OPT(AP_EtoreGobet)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_EtoreGobet)(void*Opt,void *Mod,PricingMethod *Met)
{
return AVAILABLE_IN_FULL_PREMIA;
 
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  Met->HelpFilenameHint = "ap_EtoreGobet";
  if ( Met->init == 0) Met->init=1;
  return OK;
}

PricingMethod MET(AP_EtoreGobet)=
{
  "AP_EtoreGobet",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_EtoreGobet),
  {{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_EtoreGobet),
  CHK_ok,
  MET(Init)
};
