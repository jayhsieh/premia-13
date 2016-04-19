#include  "cgmy1d_stdg.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_fft.h"
#include "math/wienerhopf.h"
#include "pnl/pnl_cdf.h"
#include"pnl/pnl_random.h"
#include"pnl/pnl_specfun.h"

static int CHK_OPT(AP_CGMY_SWING_WIENERHOPF)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}

int CALC(AP_CGMY_SWING_WIENERHOPF)(void *Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  static int first=1;

  if (first) 
    {
      Met->Par[0].Val.V_PDOUBLE=0.001;
      Met->Par[1].Val.V_PDOUBLE=1.;
      Met->Par[2].Val.V_INT2=100;
      first=0;
    }

  return OK;
}

PricingMethod MET(AP_CGMY_SWING_WIENERHOPF)=
{
  "AP_CGMY_SWING_WIENERHOPF",
  {{"Space Discretization Step",DOUBLE,{500},ALLOW},{"Scale parameter",DOUBLE,{500},ALLOW},{"TimeStepNumber",INT2,{100},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_CGMY_SWING_WIENERHOPF),
  {{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_CGMY_SWING_WIENERHOPF),
  CHK_split,
  MET(Init)
};



