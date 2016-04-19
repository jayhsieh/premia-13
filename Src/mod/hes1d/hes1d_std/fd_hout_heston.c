#include  "hes1d_std.h"

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_finance.h"
#include "pnl/pnl_tridiag_matrix.h"
#include "pnl/pnl_band_matrix.h"


int CALC(FD_Hout_Heston)(void *Opt, void *Mod, PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int CHK_OPT(FD_Hout_Heston)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
      Met->HelpFilenameHint ="fd_hout_heston";

      Met->Par[0].Val.V_INT2=5000;
      Met->Par[1].Val.V_INT2=40;
      Met->Par[2].Val.V_INT2=20;
    }

  return OK;
}

PricingMethod MET(FD_Hout_Heston)=
{
  "FD_Hout_Heston",
  {{"Time Step",INT2,{100},ALLOW},{"SpaceStepNumber S",INT2,{100},ALLOW},{"SpaceStepNumber V",INT2,{100},ALLOW}
   ,{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(FD_Hout_Heston),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(FD_Hout_Heston),
  CHK_ok,
  MET(Init)
};

