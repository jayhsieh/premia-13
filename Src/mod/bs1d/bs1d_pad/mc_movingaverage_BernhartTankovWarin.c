
#include  "bs1d_pad.h"
#include "enums.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_vector.h"

static int CHK_OPT(MC_MovingAverage_BernhartTankovWarin)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_MovingAverage_BernhartTankovWarin)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;

      Met->Par[0].Val.V_LONG=50000;
      Met->Par[1].Val.V_ENUM.value=0;
      Met->Par[1].Val.V_ENUM.members=&PremiaEnumMCRNGs;
      Met->Par[2].Val.V_ENUM.value=0;
      Met->Par[2].Val.V_ENUM.members=&PremiaEnumBasis;
      Met->Par[3].Val.V_INT=3;
    }
  return OK;
}

PricingMethod MET(MC_MovingAverage_BernhartTankovWarin)=
{
  "MC_MovingAverage_BernhartTankovWarin",
  {{"N iterations",LONG,{100},ALLOW},
      {"RandomGenerator",ENUM,{100},ALLOW},
      {"Basis",ENUM,{100},ALLOW},
      {"Approximation degree",INT,{100},ALLOW},
      {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_MovingAverage_BernhartTankovWarin),
  {{"Price",DOUBLE,{100},FORBID},
      {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_MovingAverage_BernhartTankovWarin),
  CHK_ok,
  MET(Init)
};


