#include  "hw1d_vol.h" 

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  int status=OK;
       
  return status;
}

extern PricingMethod MET(MC_Timer_HW);

PricingMethod* MOD_OPT(methods)[]={
  &MET(MC_Timer_HW),
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

