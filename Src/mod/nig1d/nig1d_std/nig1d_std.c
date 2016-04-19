#include  "nig1d_std.h" 

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);
  TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel);
  int status=OK;
	
  if ((ptOpt->Maturity.Val.V_DATE)<=(ptMod->T.Val.V_DATE))
    {
      Fprintf(TOSCREENANDFILE,"Current date greater than maturity!\n");
      status+=1;		
    };

  return status;
}
extern PricingMethod MET(AP_fastwhamerdig_nig);
extern PricingMethod MET(AP_fastwhamer_nig);
extern PricingMethod MET(AP_CarrNIG);
extern PricingMethod MET(AP_spmNIG);
extern PricingMethod MET(CF_spmNIG);
extern PricingMethod MET(FD_ImpExp);
extern PricingMethod MET(TR_MSS_NIG);

PricingMethod* MOD_OPT(methods)[]={
  &MET(AP_fastwhamer_nig),
  &MET(AP_fastwhamerdig_nig),
  &MET(FD_ImpExp),
  &MET(AP_CarrNIG),
  &MET(AP_spmNIG),
   &MET(CF_spmNIG),
  &MET(TR_MSS_NIG),	
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

