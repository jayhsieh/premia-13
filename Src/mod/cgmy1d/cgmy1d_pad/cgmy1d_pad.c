#include "cgmy1d_pad.h"

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

extern PricingMethod MET(AP_Asian_FMMCGMY);
extern PricingMethod MET(AP_FixedAsian_FusaiMeucciCGMY);
extern PricingMethod MET(MC_CGMY_FixedLookback);
extern PricingMethod MET(MC_CGMY_FloatingLookback);
extern PricingMethod MET(MC_CGMY_FixedAsian);
extern PricingMethod MET(MC_CGMY_FloatingAsian);
extern PricingMethod MET(AP_CernyKyriakou_CGMY_FixedAsian);
extern PricingMethod MET(AP_CernyKyriakou_CGMY_FloatingAsian);
extern PricingMethod MET(FFT_CGMY_FloatingLookback);
extern PricingMethod MET(FFT_CGMY_FixedLookback);
extern PricingMethod MET(AP_WH_FloatingLookback);
extern PricingMethod MET(AP_WH_FixedLookback);
PricingMethod *MOD_OPT(methods)[]={
  &MET(AP_Asian_FMMCGMY),
  &MET(AP_FixedAsian_FusaiMeucciCGMY),
  &MET(MC_CGMY_FixedLookback),
  &MET(MC_CGMY_FloatingLookback),
  &MET(MC_CGMY_FloatingAsian),
  &MET(MC_CGMY_FixedAsian),
  &MET(AP_CernyKyriakou_CGMY_FixedAsian),
  &MET(AP_CernyKyriakou_CGMY_FloatingAsian),
  &MET(FFT_CGMY_FloatingLookback),
  &MET(FFT_CGMY_FixedLookback),
  &MET(AP_WH_FloatingLookback),
  &MET(AP_WH_FixedLookback),
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


