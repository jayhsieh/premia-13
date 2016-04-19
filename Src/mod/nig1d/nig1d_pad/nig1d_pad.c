#include "nig1d_pad.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);
  TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel);
  int status=OK;

  if (ptOpt->Maturity.Val.V_DATE<=ptMod->T.Val.V_DATE)
    {
      Fprintf(TOSCREENANDFILE,"Current date greater than maturity!\n");
      status+=1;
    };
  if ((ptOpt->MinOrElse).Val.V_BOOL==MINIMUM)
    {
      if ((ptOpt->PathDep.Val.V_NUMFUNC_2)->Par[4].Val.V_PDOUBLE>ptMod->S0.Val.V_PDOUBLE)
	{
	  Fprintf(TOSCREENANDFILE,"Minimum greater than spot!\n");
	  status+=1;
	};
    }
  if ((ptOpt->MinOrElse).Val.V_BOOL==MAXIMUM)
    {
      if ((ptOpt->PathDep.Val.V_NUMFUNC_2)->Par[4].Val.V_PDOUBLE<ptMod->S0.Val.V_PDOUBLE)
	{
	  Fprintf(TOSCREENANDFILE,"Maximum lower than spot!\n");
	  status+=1;
	};
    }
  return status;
}

extern PricingMethod MET(AP_Asian_FMMNIG);
extern PricingMethod MET(AP_FixedAsian_FusaiMeucciNIG);
extern PricingMethod MET(MC_Nig_Fixed);
extern PricingMethod MET(MC_Nig_Floating);
extern PricingMethod MET(MC_NIG_FixedAsian);
extern PricingMethod MET(MC_NIG_FloatingAsian);
extern PricingMethod MET(AP_CernyKyriakou_NIG_FixedAsian);
extern PricingMethod MET(AP_CernyKyriakou_NIG_FloatingAsian);
extern PricingMethod MET(FFT_NIG_FloatingLookback);
extern PricingMethod MET(FFT_NIG_FixedLookback);
PricingMethod *MOD_OPT(methods)[]={
  &MET(AP_Asian_FMMNIG),
  &MET(AP_FixedAsian_FusaiMeucciNIG),
  &MET(MC_Nig_Fixed),
  &MET(MC_Nig_Floating),
  &MET(MC_NIG_FixedAsian),
  &MET(MC_NIG_FloatingAsian),
  &MET(AP_CernyKyriakou_NIG_FixedAsian),
  &MET(AP_CernyKyriakou_NIG_FloatingAsian),
  &MET(FFT_NIG_FloatingLookback),
  &MET(FFT_NIG_FixedLookback),
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


