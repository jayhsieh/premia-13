#include <stdlib.h>
#include  "kou1d_pad.h"
#include  "math/ap_fusai_levy/QDiscreteAsian.h"
#include "math/ap_fusai_levy/nrutil.h"

static int FusaiMeucciKOU_FixedAsian(double pseudo_stock,double pseudo_strike,NumFunc_2  *po,double t,double r,double divid,double sigma,double lambda,double lambdap,double lambdam,double p,int N,int M,double *ptprice,double *ptdelta)
{
  double CTtK,PTtK,Dlt,Plt;
  double lowlim=10.,uplim=10.;
  long int nfft=65536;
  double *price,*solution,delta;

 price=dvector(0,M-1);
 solution=dvector(0,M-1);
    
  /*  Call Price */
  CTtK= Asian_DE_FusaiMeucci(pseudo_stock,pseudo_strike,t,r,divid,sigma,lambda,p,lambdap,lambdam,N,lowlim,uplim,M,nfft,price,solution,&delta);

  /*  Put Price from Parity*/
  if(r==divid)
    PTtK=CTtK+pseudo_strike*exp(-r*t)-pseudo_stock*exp(-r*t);
  else
    PTtK=CTtK+pseudo_strike*exp(-r*t)-pseudo_stock*exp(-r*t)*(exp((r-divid)*t)-1.)/(t*(r-divid));

  /*Delta for call option*/
  Dlt=delta;

  /*Delta for put option*/
  if(r==divid)
    Plt=Dlt-exp(-r*t);
  else
    Plt=Dlt-exp(-r*t)*(exp((r-divid)*t)-1.0)/(t*(r-divid));

  /*Price*/
  if ((po->Compute)==&Call_OverSpot2)
    *ptprice=CTtK;
  else
    *ptprice=PTtK;
	 
  /*Delta */
  if ((po->Compute)==&Call_OverSpot2)
    *ptdelta=Dlt;
  else
    *ptdelta=Plt;
  
   free_dvector(price,0,M-1);
   free_dvector(solution,0,M-1);
   
  return OK;
}

int CALC(AP_FixedAsian_FusaiMeucciKOU)(void *Opt,void *Mod,PricingMethod *Met)
{
  TYPEOPT* ptOpt=(TYPEOPT*)Opt;
  TYPEMOD* ptMod=(TYPEMOD*)Mod;

  int return_value;
  double r,divid,time_spent,pseudo_spot,pseudo_strike;
  double t_0, T_0;
  
  r=log(1.+ptMod->R.Val.V_DOUBLE/100.);
  divid=log(1.+ptMod->Divid.Val.V_DOUBLE/100.);
  
  T_0 = ptMod->T.Val.V_DATE;
  t_0= (ptOpt->PathDep.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE;

  
  if(T_0 < t_0)
    {
      Fprintf(TOSCREEN,"T_0 < t_0, untreated case\n\n\n");
      return_value = WRONG;
    }
  /* Case t_0 <= T_0 */
  else
    {
      time_spent=(ptMod->T.Val.V_DATE-(ptOpt->PathDep.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE)/(ptOpt->Maturity.Val.V_DATE-(ptOpt->PathDep.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE);
      pseudo_spot=(1.-time_spent)*ptMod->S0.Val.V_PDOUBLE;
      pseudo_strike=(ptOpt->PayOff.Val.V_NUMFUNC_2)->Par[0].Val.V_PDOUBLE-time_spent*(ptOpt->PathDep.Val.V_NUMFUNC_2)->Par[4].Val.V_PDOUBLE;
  
	return_value= FusaiMeucciKOU_FixedAsian(pseudo_spot,pseudo_strike,ptOpt->PayOff.Val.V_NUMFUNC_2,ptOpt->Maturity.Val.V_DATE-ptMod->T.Val.V_DATE,r,divid,ptMod->Sigma.Val.V_PDOUBLE,ptMod->Lambda.Val.V_PDOUBLE,ptMod->LambdaPlus.Val.V_PDOUBLE,ptMod->LambdaMinus.Val.V_PDOUBLE,ptMod->P.Val.V_PDOUBLE,Met->Par[0].Val.V_INT2,Met->Par[1].Val.V_INT2,&(Met->Res[0].Val.V_DOUBLE),&(Met->Res[1].Val.V_DOUBLE));
    }
  
  return return_value;
}

static int CHK_OPT(AP_FixedAsian_FusaiMeucciKOU)(void *Opt, void *Mod)
{
  if ( (strcmp(((Option*)Opt)->Name,"AsianCallFixedEuro")==0) || (strcmp( ((Option*)Opt)->Name,"AsianPutFixedEuro")==0) )
    return OK;
  return WRONG;
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
   if ( Met->init == 0)
    {
      Met->init=1;
      Met->Par[0].Val.V_INT2=52;
      Met->Par[1].Val.V_INT2=5000;   
    }
  return OK;
}

PricingMethod MET(AP_FixedAsian_FusaiMeucciKOU)=
{
  "AP_FixedAsian_FusaiMeucci_Kou",
  {{"Nb.of Monitoring Dates",INT2,{2000},ALLOW },
   {"Nb.of Integration Points ",INT2,{1000},ALLOW},
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_FixedAsian_FusaiMeucciKOU),
  {{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID} ,{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_FixedAsian_FusaiMeucciKOU),
  CHK_ok,
  MET(Init)
};

  
