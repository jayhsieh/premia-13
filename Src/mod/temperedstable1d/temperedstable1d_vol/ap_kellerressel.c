#include "temperedstable1d_vol.h"
#include"pnl/pnl_vector.h"
#include"pnl/pnl_random.h"
#include"pnl/pnl_specfun.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_integration.h"
#include "pnl/pnl_laplace.h"
#include "pnl/pnl_fft.h"

static int CHK_OPT(AP_KELLERRESSEL)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_KELLERRESSEL)(void *Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  static int first=1;

  if (first)
    {
      first=0;
      Met->HelpFilenameHint = "ap_kellerressel";
    }
  return OK;
}

PricingMethod MET(AP_KELLERRESSEL)=
{
  "AP_KELLERRESSEL",
  { {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_KELLERRESSEL),
  {   {"Price in 10000 variance points",DOUBLE,{100},FORBID},
	  {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_KELLERRESSEL),
  CHK_ok ,
  MET(Init)
} ;
