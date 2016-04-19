#include <cmath>
#include <iostream>


extern "C"{
#include  "bsnd_stdnd.h"
}
#include "math/bsnd_math/svd_bs_tools.h"

extern "C"{    
static int CHK_OPT(FD_GreedySvd)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(FD_GreedySvd)(void *Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
	  Met->HelpFilenameHint = "FD_svd_bs";
    }

  return OK;
}

PricingMethod MET(FD_GreedySvd)=
{
  "FD_Greedy",
  {{"Discretization points per dimension",PINT,{51},ALLOW}, {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(FD_GreedySvd),
  {{"Price",DOUBLE,{100},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(FD_GreedySvd),
  CHK_split,
  MET(Init)
};

}
