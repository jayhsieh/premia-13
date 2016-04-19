
#include  "hes1d_std.h"
#include <pnl/pnl_mathtools.h>
#include <pnl/pnl_root.h>

static int CHK_OPT(AP_Asymptotics_ImpliedVolatility)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_Asymptotics_ImpliedVolatility)(void*Opt,void *Mod,PricingMethod *Met)
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

PricingMethod MET(AP_Asymptotics_ImpliedVolatility)=
{
  "AP_Asymptotics_ImpliedVolatility",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_Asymptotics_ImpliedVolatility),
  {{"Asymptotics for Implied Volatility ",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_Asymptotics_ImpliedVolatility),
  CHK_ok,
  MET(Init)
};
