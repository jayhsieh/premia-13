#include "libor_affine_cir1d_stdi.h"

int MOD_OPT(ChkMix)(Option *Opt,Model *Mod)
{
  TYPEOPT* ptOpt=( TYPEOPT*)(Opt->TypeOpt);
  TYPEMOD* ptMod=( TYPEMOD*)(Mod->TypeModel);
  int status=OK;


  if ((strcmp(Opt->Name,"Floor")==0)||(strcmp(Opt->Name,"Cap")==0))
    {
      if ((ptOpt->FirstResetDate.Val.V_DATE)<=(ptMod->T.Val.V_DATE))
	{
	  Fprintf(TOSCREENANDFILE,"Current date greater than first coupon date!\n");
	  status+=1;
	}
      if ((ptOpt->FirstResetDate.Val.V_DATE)>=(ptOpt->BMaturity.Val.V_DATE))
	{
	  Fprintf(TOSCREENANDFILE,"First reset date greater than contract maturity!\n");
	  status+=1;
	}
    }

  return status;
}


extern PricingMethod MET(CF_LibAffCir1d_Direct_CapFloor);
extern PricingMethod MET(CF_LibAffCir1d_Fourier_CapFloor);
extern PricingMethod MET(CF_LibAffCir1d_Direct_Swaption);
extern PricingMethod MET(CF_LibAffCir1d_Fourier_Swaption);

PricingMethod* MOD_OPT(methods)[]={

  &MET(CF_LibAffCir1d_Direct_CapFloor),
  &MET(CF_LibAffCir1d_Fourier_CapFloor),
  &MET(CF_LibAffCir1d_Direct_Swaption),
  &MET(CF_LibAffCir1d_Fourier_Swaption),
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
