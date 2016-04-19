#include <stdlib.h>
#include <math.h>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_fft.h"
#include "math/wienerhopf.h"
#include "cgmy1d_stdr.h"

#include "pnl/pnl_cdf.h"
#include"pnl/pnl_random.h"
#include"pnl/pnl_specfun.h"

static int CHK_OPT(AP_VAR_FFT)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(AP_VAR_FFT)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  static int first=1;

  if (first)
    {
      Met->HelpFilenameHint = "AP_CGMY_VAR";
      Met->Par[0].Val.V_PDOUBLE=2.0;
      Met->Par[1].Val.V_PDOUBLE=0.0001;
        
      first=0;
    }
  return OK;
}

PricingMethod MET(AP_VAR_FFT)=
{
  "AP_VAR_FFT",
  { {"Scale of logprice range", DOUBLE, {100}, ALLOW},
    {"Space Discretization Step",DOUBLE,{500},ALLOW},
    {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_VAR_FFT),
  {{"Value At Risk",DOUBLE,{100},FORBID},
   {"Conditional Tail Expectation ",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_VAR_FFT),
  CHK_split,
  MET(Init)
};

