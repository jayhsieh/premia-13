#include "copula_stdndc.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  return OK;
}

extern PricingMethod MET(Stein);
extern PricingMethod MET(LaurentGregory);
extern PricingMethod MET(HullWhite);
extern PricingMethod MET(MonteCarlo);
extern PricingMethod MET(Saddlepoint);

PricingMethod* MOD_OPT(methods)[]={
  &MET(Stein),
  &MET(LaurentGregory),
  &MET(HullWhite),
  &MET(MonteCarlo),
  &MET(Saddlepoint),
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

