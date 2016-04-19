#include <stdlib.h>
#include <math.h>
#include "pnl/pnl_complex.h"
#include "pnl/pnl_cdf.h"
#include "static_merton_stdr.h"


static int CHK_OPT(AP_MASDEMONTORTIZ)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_MASDEMONTORTIZ)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  static int first=1;

  if (first)
    {
      Met->HelpFilenameHint = "AP_MASDEMONTORTIZ";
      Met->Par[0].Val.V_PINT=10;        
      first=0;
    }
  return OK;
}

PricingMethod MET(AP_MASDEMONTORTIZ)=
{
  "AP_MASDEMONTORTIZ",
  { {"Approximation Scale",INT, {100}, ALLOW},
    {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_MASDEMONTORTIZ),
  {{"Value At Risk",DOUBLE,{100},FORBID},
   {"Conditional Tail Expectation ",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_MASDEMONTORTIZ),
  CHK_split,
  MET(Init)
};

