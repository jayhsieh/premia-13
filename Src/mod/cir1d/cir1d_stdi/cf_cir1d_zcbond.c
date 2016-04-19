#include  "cir1d_stdi.h"

/*Zero Coupon Bond*/
static int zcb_cir1d(double r0,double k, double t,double sigma,double theta,double T,double *price)
{
  double h,A,B;

  /*A,B coefficient*/
  h=sqrt(SQR(k)+2.*SQR(sigma));
  B=2.*(exp(h*(T-t))-1.)/(2.*h+(k+h)*(exp(h*(T-t))-1.));
  A=pow(h*exp(0.5*(k+h)*(T-t))/(h+0.5*(k+h)*(exp(h*(T-t))-1.)),2.*k*theta/SQR(sigma));
  
  /*Price*/
  *price=A*exp(-B*r0);

  return OK;
}

int CALC(CF_ZCBond)(void *Opt,void *Mod,PricingMethod *Met)
{
  TYPEOPT* ptOpt=(TYPEOPT*)Opt;
  TYPEMOD* ptMod=(TYPEMOD*)Mod;

  return zcb_cir1d(ptMod->r0.Val.V_PDOUBLE,ptMod->k.Val.V_DOUBLE,ptMod->T.Val.V_DATE,ptMod->Sigma.Val.V_PDOUBLE,ptMod->theta.Val.V_PDOUBLE,ptOpt->BMaturity.Val.V_DATE,&(Met->Res[0].Val.V_DOUBLE));
}

static int CHK_OPT(CF_ZCBond)(void *Opt, void *Mod)
{
  return strcmp( ((Option*)Opt)->Name,"ZeroCouponBond");
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(CF_ZCBond)=
{
  "CF_Cir1d_ZCBond",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(CF_ZCBond),
  {{"Price",DOUBLE,{100},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(CF_ZCBond),
  CHK_ok,
  MET(Init)
} ;


