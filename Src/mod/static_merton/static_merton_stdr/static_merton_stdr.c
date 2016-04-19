#include  "static_merton_stdr.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  return OK;
}

extern PricingMethod MET(AP_MASDEMONTORTIZ);
PricingMethod* MOD_OPT(methods)[]={
  &MET(AP_MASDEMONTORTIZ),	
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
