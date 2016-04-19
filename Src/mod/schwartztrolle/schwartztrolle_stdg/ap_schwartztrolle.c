#include  "schwartztrolle_stdg.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_cdf.h"

static int CHK_OPT(AP_SCHWARTZTROLLE)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_SCHWARTZTROLLE)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  //int type_generator;
  if ( Met->init == 0)
    {
      Met->init=1;   
      Met->HelpFilenameHint = "ap_schwartztrolle";
    }
  
  return OK;
}

PricingMethod MET(AP_SCHWARTZTROLLE)=
{
  "AP_SCHWARTZTROLLE",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_SCHWARTZTROLLE),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_SCHWARTZTROLLE),
  CHK_ok,
  MET(Init)
};
