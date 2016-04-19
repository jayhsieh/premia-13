
#include  "hhw4d_stdf.h"
#include <math.h>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_fft.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_integration.h"

static int CHK_OPT(AP_HHW4D_INF)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_HHW4D_INF)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  //int type_generator;
  if ( Met->init == 0)
    {
      Met->init=1;   
      Met->HelpFilenameHint = "ap_grzlek_inflation";
    }
  
  return OK;
}

PricingMethod MET(AP_HHW4D_INF)=
{
  "AP_HHW4D_INF",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_HHW4D_INF),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_HHW4D_INF),
  CHK_ok,
  MET(Init)
};
