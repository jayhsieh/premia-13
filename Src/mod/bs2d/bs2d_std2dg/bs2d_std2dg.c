#include  "bs2d_std2dg.h" 

static int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  TYPEOPT* ptOpt=(TYPEOPT*)(Opt->TypeOpt);
  TYPEMOD* ptMod=(TYPEMOD*)(Mod->TypeModel);
  int status=OK;
	
  if (ptOpt->Maturity.Val.V_DATE<=ptMod->T.Val.V_DATE)
    {
      Fprintf(TOSCREENANDFILE,"Current date greater than maturity!\n");
      status+=1;		
    };

	
  return status;
}
extern PricingMethod MET(AP_Spread_Carmona);
extern PricingMethod MET(AP_Spread_Bjerksund);
extern PricingMethod MET(AP_Spread_HurdZhou_bs);
PricingMethod* MOD_OPT(methods)[]={
  &MET(AP_Spread_Carmona),
   &MET(AP_Spread_Bjerksund),
 &MET(AP_Spread_HurdZhou_bs),
NULL
};

//extern DynamicTest MOD_OPT(test);
DynamicTest* MOD_OPT(tests)[]={
  //&MOD_OPT(test),
  NULL
};

Pricing MOD_OPT(pricing)={
  ID_MOD_OPT,
  MOD_OPT(methods),
  MOD_OPT(tests),
  MOD_OPT(ChkMix)
};



