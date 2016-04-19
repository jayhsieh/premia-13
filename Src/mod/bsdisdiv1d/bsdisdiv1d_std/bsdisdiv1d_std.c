#include  "bsdisdiv1d_std.h" 

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

extern PricingMethod MET(AP_EtoreGobet);
extern PricingMethod MET(TR_Vellekoop);
extern PricingMethod MET(TR_SingularPoints_Up);
extern PricingMethod MET(TR_SingularPoints_Down);

PricingMethod* MOD_OPT(methods)[]={
 
  &MET(AP_EtoreGobet),
  &MET(TR_Vellekoop),
  &MET(TR_SingularPoints_Up),
  &MET(TR_SingularPoints_Down),
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






























