#include <stdlib.h>
#include  "kou1d_pad.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_fft.h"
#include  "math/ap_fusai_levy/DiscreteAsianFMM.h"

static int CHK_OPT(AP_Asian_FMMKOU)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_Asian_FMMKOU)(void *Opt,void *Mod,PricingMethod *Met)
{
return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
   if ( Met->init == 0)
    {
      Met->init=1;
      Met->Par[0].Val.V_INT2=52;
      Met->Par[1].Val.V_INT2=3000;   
    }
  return OK;
}

PricingMethod MET(AP_Asian_FMMKOU)=
{
  "AP_Asian_FMM_KOU",
  {{"Nb.of Monitoring Dates",INT2,{2000},ALLOW },
   {"Nb.of Integration Points ",INT2,{1000},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_Asian_FMMKOU),
  {{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID} ,{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_Asian_FMMKOU),
  CHK_ok,
  MET(Init)
};

  
