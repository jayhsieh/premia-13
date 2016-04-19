#include "mer1d_std.h"
#include "error_msg.h"
#include <math.h>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

static int CHK_OPT(AP_STATICHEDGING_CARRWU)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}

int CALC(AP_STATICHEDGING_CARRWU)(void*Opt,void *Mod,PricingMethod *Met)
{
return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
      Met->Par[0].Val.V_DATE=0.5;
      Met->Res[0].Val.V_PNLVECT=NULL;
      Met->Res[1].Val.V_PNLVECT=NULL;
    }

  /* some initialisation */
  if(Met->Res[0].Val.V_PNLVECT==NULL)
    Met->Res[0].Val.V_PNLVECT=pnl_vect_create(5);
  else
    pnl_vect_resize(Met->Res[0].Val.V_PNLVECT,5);

  if(Met->Res[1].Val.V_PNLVECT==NULL)
    Met->Res[1].Val.V_PNLVECT=pnl_vect_create(5);
  else
    pnl_vect_resize(Met->Res[1].Val.V_PNLVECT,5);
  
  return OK;
}

PricingMethod MET(AP_STATICHEDGING_CARRWU)=
{
  "AP_STATICHEDGING_CARRWU",
  {{"Hedging Maturity",DATE,{100},ALLOW},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_STATICHEDGING_CARRWU),
  {{"Strikes",PNLVECT,{100},FORBID},{"Strikes Weights",PNLVECT,{1},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_STATICHEDGING_CARRWU),
  CHK_ok,
  MET(Init)
} ;
