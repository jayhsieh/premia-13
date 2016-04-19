#include "svj.h"
#include "grad_svj.h"
#include "hes1d_std.h"
//
//static
int Grad_FT_Call_HestMert(double St0, NumFunc_1  *p, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v,int dimx,double *grad)
{
  double K;
  SVJPARAMS *svj;
  //
  K=p->Par[0].Val.V_DOUBLE;
  //
  svj = (SVJPARAMS *)malloc(sizeof(SVJPARAMS));
  //
  svj->heston = 1;
  svj->merton = 1;
  // phi = 1 ==> Call
  svj->phi    = 1.;
  //
  svj->type_f = 1;
  //
  svj->K      = K;
  svj->St0    = St0;
  svj->T      = T;
  svj->r      = r;
  svj->divid  = divid;
  //
  svj->sigmav = sigmav;
  svj->V0     = V0;
  svj->theta  = theta;
  svj->rho    = rho;
  svj->kappa  = kappa;
  svj->lambda = lambda;
  svj->m0     = m0;
  svj->v      = v;
  //       
  calc_grad_svj(svj,dimx,grad);
  //  
  free(svj);
  //
  return OK;
}
/*
int CALC(CF_CallHeston)(void *Opt, void *Mod, PricingMethod *Met)
{
  TYPEOPT* ptOpt=(TYPEOPT*)Opt;
  TYPEMOD* ptMod=(TYPEMOD*)Mod;
  double r,divid;

   if(ptMod->Sigma.Val.V_PDOUBLE==0.0)
   {
  Fprintf(TOSCREEN,"BLACK-SHOLES MODEL\n\n\n");
   return WRONG;
   }
   else 
     {
  r=log(1.+ptMod->R.Val.V_DOUBLE/100.);
  divid=log(1.+ptMod->Divid.Val.V_DOUBLE/100.);

  return CFCallHeston(ptMod->S0.Val.V_PDOUBLE,
		  ptOpt->PayOff.Val.V_NUMFUNC_1,
		  ptOpt->Maturity.Val.V_DATE-ptMod->T.Val.V_DATE,
		  r,
		  divid, ptMod->Sigma0.Val.V_PDOUBLE
		  ,ptMod->MeanReversion.Val.V_PDOUBLE,
		  ptMod->LongRunVariance.Val.V_PDOUBLE,
		  ptMod->Sigma.Val.V_PDOUBLE,
		  ptMod->Rho.Val.V_PDOUBLE,
		  &(Met->Res[0].Val.V_DOUBLE),
		  &(Met->Res[1].Val.V_DOUBLE)
		  );
  }
		    
}



int CHK_OPT(CF_CallHeston)(void *Opt, void *Mod)
{ 
return strcmp( ((Option*)Opt)->Name,"CallEuro");   
}



static int MET(Init)(PricingMethod *Met)
{
  return OK;
}

PricingMethod MET(CF_CallHeston)=
{
  "CF_Heston",
  {{" ",END,0,FORBID}},
  CALC(CF_CallHeston),
  {{"Price",DOUBLE,100,FORBID},
   {"Delta",DOUBLE,100,FORBID} ,
   {" ",END,0,FORBID}},
  CHK_OPT(CF_CallHeston),
  CHK_ok,
  MET(Init)
};
*/
