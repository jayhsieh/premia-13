#include  "nig1fact1d_std.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  int status=OK;

  return status;
}

extern PricingMethod MET(AP_GOR); 

PricingMethod* MOD_OPT(methods)[]={
   &MET(AP_GOR),
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
