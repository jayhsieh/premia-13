#include  "doublehes1d_std.h"
#include "pnl/pnl_integration.h"
#include "pnl/pnl_finance.h"

static int CHK_OPT(CF_Carr_DoubleHeston)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(CF_Carr_DoubleHeston)(void*Opt,void *Mod,PricingMethod *Met)
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

PricingMethod MET(CF_Carr_DoubleHeston)=
{
  "CF_Carr_DoubleHeston",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(CF_Carr_DoubleHeston),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(CF_Carr_DoubleHeston),
  CHK_ok,
  MET(Init)
};

