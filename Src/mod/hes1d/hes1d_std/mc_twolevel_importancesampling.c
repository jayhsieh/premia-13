#include <stdlib.h>
#include  "hes1d_std.h"
#include "pnl/pnl_matrix.h"
#include "enums.h"

static int CHK_OPT(MC_TwoLevel_ImportanceSampling)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_TwoLevel_ImportanceSampling)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}


static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;

      Met->Par[0].Val.V_INT=400;
      Met->Par[1].Val.V_ENUM.value=0;
      Met->Par[1].Val.V_ENUM.members=&PremiaEnumRNGs;
      Met->Par[2].Val.V_DOUBLE= 0.95;      
    }
  return OK;
}


PricingMethod MET(MC_TwoLevel_ImportanceSampling)=
{
  "MC_TwoLevel_ImportanceSampling",
  {
	{"TimeStepNumber",LONG,{100},ALLOW},
	{"RandomGenerator",ENUM,{100},ALLOW},
	{"Confidence Value",DOUBLE,{100},ALLOW},
	{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_TwoLevel_ImportanceSampling),
  {{"Price",DOUBLE,{100},FORBID},
    {"Delta",DOUBLE,{100},FORBID} ,
   {"Inf Price",DOUBLE,{100},FORBID},
   {"Sup Price",DOUBLE,{100},FORBID} ,
    {"Inf Delta",DOUBLE,{100},FORBID},
   {"Sup Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_TwoLevel_ImportanceSampling),
  CHK_mc,
  MET(Init)
};
