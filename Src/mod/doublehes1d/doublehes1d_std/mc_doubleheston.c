#include  "doublehes1d_std.h"
#include  "math/alfonsi.h"
#include "pnl/pnl_random.h"

static int CHK_OPT(MC_DoubleHeston)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_DoubleHeston)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static PremiaEnumMember DoubleHestonOrder_members[] = 
  {
    { "Exact Zhu",  1},
    { "Second Order Alfonsi", 2},
    { "Third Order Alfonsi", 3},
    { "Quadratic Exponential Martingale", 4},
    { "Euler scheme", 5},
    { NULL, NULLINT }
  };

static DEFINE_ENUM(DoubleHestonOrder,DoubleHestonOrder_members);

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
      
      Met->Par[0].Val.V_LONG=100000;
      Met->Par[1].Val.V_INT=24;
      Met->Par[2].Val.V_ENUM.value=0;
      Met->Par[2].Val.V_ENUM.members=&PremiaEnumMCRNGs;
      Met->Par[3].Val.V_ENUM.value=1;
      Met->Par[3].Val.V_ENUM.members=&DoubleHestonOrder;
    }

  return OK;
}

PricingMethod MET(MC_DoubleHeston)=
{
  "MC_DoubleHeston",
  {{"N iterations",LONG,{100},ALLOW},
   {"TimeStepNumber",LONG,{100},ALLOW},
   {"RandomGenerator",ENUM,{100},ALLOW},
   {"Scheme type",ENUM,{100},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_DoubleHeston),
  {{"Price",DOUBLE,{100},FORBID},
   {"Inf Price",DOUBLE,{100},FORBID},
   {"Sup Price",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_DoubleHeston),
  CHK_ok,
  MET(Init)
};

