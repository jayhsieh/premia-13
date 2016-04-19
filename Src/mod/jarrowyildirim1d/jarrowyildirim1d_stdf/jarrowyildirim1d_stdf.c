#include  "jarrowyildirim1d_stdf.h" 

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);
  TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel);
  int status=OK;
  
  if ((ptOpt->BMaturity.Val.V_DATE)<=(ptMod->T.Val.V_DATE))
    {
      Fprintf(TOSCREENANDFILE,"Current date greater than maturity!\n");
      status+=1;		
    };

 
  return status;
}

extern PricingMethod MET(CF_YI_YYIIS);
extern PricingMethod MET(CF_YI_IICAPLET);
PricingMethod* MOD_OPT(methods)[]={
  &MET(CF_YI_YYIIS),
   &MET(CF_YI_IICAPLET),
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

