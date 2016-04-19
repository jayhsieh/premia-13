#include  "schwartztrolle_stdg.h" 

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);
  //TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel);
  int status=OK;

  if ((strcmp(Opt->Name,"CallFuture")==0)||(strcmp(Opt->Name,"PutFuture")==0))
	   if((ptOpt->FutureMaturity.Val.V_DATE)<(ptOpt->Maturity.Val.V_DATE))
	{
	  Fprintf(TOSCREENANDFILE,"Forward maturity greater than Option aturity!\n");
	  status+=1;
	}

  return status;
}

extern PricingMethod MET(AP_SCHWARTZTROLLE);
PricingMethod* MOD_OPT(methods)[]={
  &MET(AP_SCHWARTZTROLLE),
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

