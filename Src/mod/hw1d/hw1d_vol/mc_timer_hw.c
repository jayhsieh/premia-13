#include  "hw1d_vol.h"
#include"pnl/pnl_vector.h"
#include"pnl/pnl_random.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_cdf.h"

static int CHK_OPT(MC_Timer_HW)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_Timer_HW)(void *Opt, void *Mod, PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  //int type_generator;
  if ( Met->init == 0)
    {
      Met->init=1;
      Met->HelpFilenameHint = "mc_timer_hw";

      Met->Par[0].Val.V_LONG=15000;
      Met->Par[1].Val.V_INT=100;
      Met->Par[2].Val.V_ENUM.value=0;
      Met->Par[2].Val.V_ENUM.members=&PremiaEnumMCRNGs;
    }  
  return OK;
}

PricingMethod MET(MC_Timer_HW)=
{
  "MC_Timer_HW",
  {{"N iterations",LONG,{100},ALLOW},
   {"TimeStepNumber",LONG,{100},ALLOW},
   {"RandomGenerator",ENUM,{100},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_Timer_HW),
  {{"Price",DOUBLE,{100},FORBID},
   {"Inf Price",DOUBLE,{100},FORBID},
   {"Sup Price",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_Timer_HW),
  CHK_mc,
  MET(Init)
};
