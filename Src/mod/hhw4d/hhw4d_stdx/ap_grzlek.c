
#include  "hhw4d_stdx.h"
#include <math.h>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_fft.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_integration.h"

static int CHK_OPT(AP_HHW4D)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_HHW4D)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  //int type_generator;
  if ( Met->init == 0)
    {
      Met->init=1;   
      Met->HelpFilenameHint = "ap_grzlek";
    }
  
  return OK;
}

PricingMethod MET(AP_HHW4D)=
{
  "AP_HHW4D",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_HHW4D),
  {{"Price",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_HHW4D),
  CHK_ok,
  MET(Init)
};
