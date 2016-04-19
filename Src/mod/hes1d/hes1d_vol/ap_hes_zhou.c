#include "hes1d_vol.h"
#include "numfunc.h"
#include "pnl/pnl_mathtools.h"

static int CHK_OPT(AP_HES_VS_ZHOU)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_HES_VS_ZHOU)(void *Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  Met->HelpFilenameHint ="ap_hes_zhou";
  return OK;
}

PricingMethod MET(AP_HES_VS_ZHOU)=
{
  "AP_HES_VS_ZHOU", 
  {   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_HES_VS_ZHOU),
  {  {"Price in 10000 variance points",DOUBLE,{100},FORBID},{"Fair strike for variance swap",DOUBLE,{100},FORBID},
	 {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_HES_VS_ZHOU),
  CHK_ok ,
  MET(Init)
} ;
