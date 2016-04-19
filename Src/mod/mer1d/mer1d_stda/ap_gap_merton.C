extern "C"{
#include  "mer1d_stda.h"
#include "enums.h"
#include "pnl/pnl_finance.h"
}
#include "math/levy.h"
#include "math/fft.h"

extern "C"{

  static int CHK_OPT(AP_GAP_MERTON)(void *Opt, void *Mod)
  {
	return NONACTIVE; 
  }
  int CALC(AP_GAP_MERTON)(void *Opt,void *Mod,PricingMethod *Met)
  {
	return AVAILABLE_IN_FULL_PREMIA;
  }

  static int MET(Init)(PricingMethod *Met,Option *Opt)
  {
	if ( Met->init == 0)
	  {
		Met->init=1;
	  }
	return OK;
  }

  PricingMethod MET(AP_GAP_MERTON)=
  {
	"AP_GAP_MERTON",
	{{" ",PREMIA_NULLTYPE,{0},FORBID}},
	CALC(AP_GAP_MERTON),
	{{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
	CHK_OPT(AP_GAP_MERTON),
	CHK_split,
	MET(Init)
  };

}


