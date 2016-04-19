#include  "mrc30d_stdnd.h"
#include "enums.h"
#include "error_msg.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_finance.h"
#include <string.h>

static int CHK_OPT(MC_BASKET30D)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_BASKET30D)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  //int type_generator;
  if ( Met->init == 0)
    {
      Met->init=1;

      Met->Par[0].Val.V_LONG=15000;
      Met->Par[1].Val.V_INT=10;
      Met->Par[2].Val.V_ENUM.value=0;
      Met->Par[2].Val.V_ENUM.members=&PremiaEnumMCRNGs;
    }
  
  return OK;
}

PricingMethod MET(MC_BASKET30D)=
{
  "MC_BASKET30D",
  {{"N iterations",LONG,{100},ALLOW},
   {"TimeStepNumber",LONG,{100},ALLOW},
   {"RandomGenerator",ENUM,{100},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_BASKET30D),
  {{"Price",DOUBLE,{100},FORBID},
   {"Error Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_BASKET30D),
  CHK_mc,
  MET(Init)
};
