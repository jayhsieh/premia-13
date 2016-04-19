#include  "lmm1d_cgmy_stdi.h"
#include "enums.h"
#include"pnl/pnl_vector.h"
#include"pnl/pnl_random.h"
#include"pnl/pnl_specfun.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_integration.h"

static int CHK_OPT(MC_LOGLEVY_SWAPTION)(void *Opt, void *Mod)
{
  return NONACTIVE;
}
int CALC(MC_LOGLEVY_SWAPTION)(void *Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static PremiaEnumMember skovmand_method_members[] = 
{
    { "Log Levy",1},
    { "Second Order Drift Expansion",2},
    { NULL, NULLINT }
};

static DEFINE_ENUM(skovmand_method,skovmand_method_members);

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;

      Met->Par[0].Val.V_ENUM.value=0;
      Met->Par[0].Val.V_ENUM.members=&PremiaEnumRNGs;
      Met->Par[1].Val.V_LONG=1000;
      Met->Par[2].Val.V_ENUM.value=1;
      Met->Par[2].Val.V_ENUM.members=&skovmand_method;
    }
  return OK;
}

PricingMethod MET(MC_LOGLEVY_SWAPTION)=
{
  "MC_LogLevy_Swaption",
  {
      {"RandomGenerator",ENUM,{100},ALLOW},
      {"N Simulation",LONG,{100},ALLOW},
      {"Method",ENUM,{100},ALLOW},
      {" ",PREMIA_NULLTYPE,{0},FORBID}
  },
  CALC(MC_LOGLEVY_SWAPTION),
  {
      {"Price",DOUBLE,{100},FORBID},
      {"Price Error",DOUBLE,{100},FORBID},
      {" ",PREMIA_NULLTYPE,{0},FORBID}
  },
  CHK_OPT(MC_LOGLEVY_SWAPTION),
  CHK_ok,
  MET(Init)
} ;
