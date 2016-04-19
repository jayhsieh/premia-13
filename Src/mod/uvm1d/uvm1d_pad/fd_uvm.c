#include  "uvm1d_pad.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_tridiag_matrix.h"

static int CHK_OPT(FD_UVM)(void *Opt, void *Mod)
{
  return NONACTIVE; 
}
int CALC(FD_UVM)(void*Opt,void *Mod,PricingMethod *Met)
{
  return AVAILABLE_IN_FULL_PREMIA;
}
static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;

      Met->Par[0].Val.V_INT2=100;
      Met->Par[1].Val.V_INT2=64;
      Met->Par[2].Val.V_INT2=20;
    

    }

  return OK;
}

PricingMethod MET(FD_UVM)=
{
  "FD_UVM",
  {{"TimeStepforYear",INT2,{100},ALLOW} ,{"SpaceStep S",INT2,{100},ALLOW},{"SpaceStep Z",INT2,{100},ALLOW},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(FD_UVM),
  {{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(FD_UVM),
  CHK_ok,
  MET(Init)
};
