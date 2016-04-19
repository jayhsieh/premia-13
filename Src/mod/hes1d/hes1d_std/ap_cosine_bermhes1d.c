#include "hes1d_std.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_fft.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_specfun.h"
#include "pnl/pnl_root.h"


static int CHK_OPT(AP_CosineBermudan)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_CosineBermudan)(void *Opt, void *Mod, PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0 )
    {
      Met->Par[0].Val.V_INT=20;
      Met->init = 1;
      Met->HelpFilenameHint = "ap_cosine_bermhes1d";
    }
  return OK;
}


PricingMethod MET(AP_CosineBermudan)=
{
  "AP_CosineBermudan",
  {{"Bermudan Steps",INT,{100},ALLOW},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_CosineBermudan),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_CosineBermudan),
  CHK_ok,
  MET(Init)
};

