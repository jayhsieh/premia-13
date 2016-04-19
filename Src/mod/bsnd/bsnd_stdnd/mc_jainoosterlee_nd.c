//Calcul du prix d'une option américaine sur maximum par l'algo sgm en dimension 2. Pas de
//polynomes locaux. Ce code fonctionne


#include  "bsnd_stdnd.h"
#include "pnl/pnl_random.h"
#include "pnl/pnl_basis.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_cdf.h"

static int CHK_OPT(MC_JainOosterleeND)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(MC_JainOosterleeND)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
      Met->HelpFilenameHint = "MC_JainOosterlee_ND";
      Met->Par[0].Val.V_LONG=50000;
      Met->Par[1].Val.V_ENUM.value=0;
      Met->Par[1].Val.V_ENUM.members=&PremiaEnumMCRNGs;
      Met->Par[2].Val.V_INT=6;
      Met->Par[3].Val.V_INT=10;
    }
  return OK;
}

PricingMethod MET(MC_JainOosterleeND)=
{
  "MC_JainOosterleeND",
  {{"N iterations",LONG,{100},ALLOW},
   {"RandomGenerator",ENUM,{0},ALLOW},
   {"Dimension Approximation",INT,{100},ALLOW},
   {"Number of Exercise Dates",INT,{100},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(MC_JainOosterleeND),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(MC_JainOosterleeND),
  CHK_mc,
  MET(Init)
};

