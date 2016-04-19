#include "bsnd_stdnd.h"
#include "pnl/pnl_cdf.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"


static int CHK_OPT(AP_MultiSpread)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}

int CALC(AP_MultiSpread)(void *Opt, void *Mod, PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  TYPEOPT *opt = (TYPEOPT*)(Opt->TypeOpt);

  if ( Met->init == 0)
    {
      Met->init=1;
      Met->Res[1].Val.V_PNLVECT=NULL;
    }
  /* some initialisation */
  if(Met->Res[1].Val.V_PNLVECT==NULL)
    Met->Res[1].Val.V_PNLVECT=pnl_vect_create(opt->Size.Val.V_PINT);
  else
    pnl_vect_resize(Met->Res[1].Val.V_PNLVECT,opt->Size.Val.V_PINT);

  return OK;

}

PricingMethod MET(AP_MultiSpread)=
{
  "AP_MultiSpread",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_MultiSpread),
  {{"Price",DOUBLE,{100},FORBID},{"Deltas",PNLVECT,{1},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_MultiSpread),
  CHK_ok,
  MET(Init)
};


