#include "dynamic_stdndc.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  return OK;
}

extern PricingMethod MET(RogersDiGraziano);
extern PricingMethod MET(Hedging_CousinFermanianLaurent);
extern PricingMethod MET(Hedging_FreyBackhaus);
extern PricingMethod MET(EberleinFreyVHammerstein);
extern PricingMethod MET(Herbertsson);

PricingMethod* MOD_OPT(methods)[]={
  &MET(RogersDiGraziano),
  &MET(Hedging_CousinFermanianLaurent),
  &MET(Hedging_FreyBackhaus),
  &MET(EberleinFreyVHammerstein),
  &MET(Herbertsson),
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

