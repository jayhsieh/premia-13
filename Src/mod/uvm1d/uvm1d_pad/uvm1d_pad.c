#include  "uvm1d_pad.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  /* TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);
   * TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel); */
  int status=OK;

  return status;
}

extern PricingMethod MET(FD_UVM); 

PricingMethod* MOD_OPT(methods)[]={
 &MET(FD_UVM),
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
