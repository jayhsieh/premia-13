#include  "schwartz_stdi.h" 

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);
  //TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel);
  int status=OK;

  if ((strcmp(Opt->Name,"PayerSwaption")==0))
  	   if((ptOpt->BMaturity.Val.V_DATE)<=(ptOpt->OMaturity.Val.V_DATE))
  	{
  	  Fprintf(TOSCREENANDFILE,"Option maturity greater than Bond maturity!\n");
  	  status+=1;
  	}
   
  return status;
}

extern PricingMethod MET(AP_SCHWARTZ);
PricingMethod* MOD_OPT(methods)[]={
  &MET(AP_SCHWARTZ),
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

