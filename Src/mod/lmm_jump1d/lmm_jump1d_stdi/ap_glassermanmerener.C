/*--------------------------------------------------------------*/
/*   CF approx. for caplet prices in one-factor LMM with jumps  */
/*   Algorithm of Glasserman/Merener                            */
/*                                                              */
/*--------------------------------------------------------------*/
/*  Sonke Blunck, Premia 2005                                   */
/*--------------------------------------------------------------*/

#include  "glassermanmerener.h"

extern "C"{
#include  "lmm_jump1d_stdi.h"


static int ap_glassermanmerenenr_caplet(NumFunc_1 *p,double l0,double t0, double sigma,double capletMat ,double strike,double tenor,double *price)
{
     
  capletMat=capletMat-t0;
  return lmm_jump_caplet_GlassMer_pricer(tenor,capletMat,strike,l0,sigma, price);

}

int CALC(AP_GM)(void *Opt,void *Mod,PricingMethod *Met)
{
  TYPEOPT* ptOpt=(TYPEOPT*)Opt;
  TYPEMOD* ptMod=(TYPEMOD*)Mod;

  return ap_glassermanmerenenr_caplet(ptOpt->PayOff.Val.V_NUMFUNC_1,ptMod->l0.Val.V_PDOUBLE,
                                      ptMod->T.Val.V_DATE,
                                      ptMod->Sigma.Val.V_PDOUBLE,
                                      ptOpt->BMaturity.Val.V_DATE,
                                      ptOpt->FixedRate.Val.V_PDOUBLE,
                                      ptOpt->ResetPeriod.Val.V_DATE,
                                      &(Met->Res[0].Val.V_DOUBLE));
}

static int CHK_OPT(AP_GM)(void *Opt, void *Mod)
{

  if ((strcmp(((Option*)Opt)->Name,"Caplet")==0))
    return OK;
  else
    return WRONG;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(AP_GM)=
{
  "AP_GlassermanMerener",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_GM),
  {{"Price",DOUBLE,{100},FORBID}/*,{"Delta",DOUBLE,{100},FORBID}*/ ,{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_GM),
  CHK_ok,
  MET(Init)
} ;
}
