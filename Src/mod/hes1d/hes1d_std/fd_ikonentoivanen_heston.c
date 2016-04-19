#include  "hes1d_std.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "pnl/pnl_finance.h"



static int CHK_OPT(FD_IkonenToivanen_Heston)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}

int CALC(FD_IkonenToivanen_Heston)(void *Opt, void *Mod, PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
      Met->HelpFilenameHint ="fd_ikonentoivanen_heston";

      Met->Par[0].Val.V_INT2=1000;
      Met->Par[1].Val.V_INT2=101;
      Met->Par[2].Val.V_INT2=50;
    }

  return OK;
}

PricingMethod MET(FD_IkonenToivanen_Heston)=
{
  "FD_IkonenToivanen_Heston",
  {{"Time Step",INT2,{100},ALLOW},{"SpaceStepNumber S",INT2,{100},ALLOW},{"SpaceStepNumber V",INT2,{100},ALLOW}
,{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(FD_IkonenToivanen_Heston),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(FD_IkonenToivanen_Heston),
  CHK_ok,
  MET(Init)
};

