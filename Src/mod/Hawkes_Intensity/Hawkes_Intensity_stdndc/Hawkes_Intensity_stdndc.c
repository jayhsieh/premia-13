#include  "Hawkes_Intensity_stdndc.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  int status=OK;

  return status;
}

extern PricingMethod MET(ErraisGieseckeGoldberg);

PricingMethod* MOD_OPT(methods)[]={
  &(MET(ErraisGieseckeGoldberg)),
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
