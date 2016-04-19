#include  "scott1d_std.h"
#include "enums.h"
#include "error_msg.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_root.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_random.h"

static int CHK_OPT(MC_MultiLevel_Scott)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_MultiLevel_Scott)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;

      Met->Par[0].Val.V_LONG=1000;
      Met->Par[1].Val.V_ENUM.value=0;
      Met->Par[1].Val.V_ENUM.members=&PremiaEnumMCRNGs;
      Met->Par[2].Val.V_PDOUBLE=0.05;
    }

  return OK;
}

PricingMethod MET(MC_MultiLevel_Scott)=
{
  "MC_MultiLevel",
  {{"N iterations",LONG,{100},ALLOW},
   {"RandomGenerator",ENUM,{100},ALLOW},
   {"Accuracy",DOUBLE,{100},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_MultiLevel_Scott),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_MultiLevel_Scott),
  CHK_mc,
  MET(Init)
};
