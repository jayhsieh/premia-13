#include  "cgmy1d_stdr.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  return OK;
}

extern PricingMethod MET(AP_VAR_FFT);
PricingMethod* MOD_OPT(methods)[]={
  &MET(AP_VAR_FFT),	
NULL
};

DynamicTest* MOD_OPT(tests)[]={
NULL
};

Pricing MOD_OPT(pricing)={
ID_MOD_OPT,
MOD_OPT(methods),
MOD_OPT(tests),
MOD_OPT(ChkMix)
};
