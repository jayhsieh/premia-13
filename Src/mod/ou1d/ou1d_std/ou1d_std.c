#include  "ou1d_std.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  int status=OK;

  return status;
}

extern PricingMethod MET(AP_FOURIERCOSINE_OU); 

PricingMethod* MOD_OPT(methods)[]={
   &MET(AP_FOURIERCOSINE_OU),
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
