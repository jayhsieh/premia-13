#include <stdlib.h>
#include "nig1d_pad.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_fft.h"
#include  "math/ap_fusai_levy/DiscreteAsianFMM.h"

int CALC(AP_Asian_FMMNIG)(void *Opt,void *Mod,PricingMethod *Met)
{
return AVAILABLE_IN_FULL_PREMIA;
}

static int CHK_OPT(AP_Asian_FMMNIG)(void *Opt, void *Mod)
{
    if ( (strcmp(((Option*)Opt)->Name,"AsianCallFixedEuro")==0) || (strcmp( ((Option*)Opt)->Name,"AsianPutFixedEuro")==0)|| (strcmp(((Option*)Opt)->Name,"AsianCallFloatingEuro")==0) || (strcmp( ((Option*)Opt)->Name,"AsianPutFloatingEuro")==0) )
    return OK;
  return WRONG;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
   if ( Met->init == 0)
    {
      Met->init=1;
      Met->Par[0].Val.V_INT2=12;
      Met->Par[1].Val.V_INT2=3000;   
    }
  return OK;
}

PricingMethod MET(AP_Asian_FMMNIG)=
{
  "AP_Asian_FMM_NIG",
  {{"Nb.of Monitoring Dates",INT2,{2000},ALLOW },
   {"Nb.of Integration Points ",INT2,{1000},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_Asian_FMMNIG),
    {{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID} ,{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_Asian_FMMNIG),
  CHK_ok,
  MET(Init)
};

  
