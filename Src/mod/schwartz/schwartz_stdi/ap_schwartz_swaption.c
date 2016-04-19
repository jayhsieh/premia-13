#include  "schwartz_stdi.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_cdf.h"

static int CHK_OPT(AP_SCHWARTZ)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_SCHWARTZ)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  //int type_generator;
  if ( Met->init == 0)
    {
      Met->init=1;   
      Met->HelpFilenameHint ="ap_schwartz_swaption";
    }
  
  return OK;
}

PricingMethod MET(AP_SCHWARTZ)=
{
  "AP_LARSSON",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_SCHWARTZ),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_SCHWARTZ),
  CHK_ok,
  MET(Init)
};
