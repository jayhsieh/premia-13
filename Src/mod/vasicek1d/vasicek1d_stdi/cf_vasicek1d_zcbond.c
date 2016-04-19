#include  "vasicek1d_stdi.h"

/*Zero Coupon Bond*/
static int zcb_vasicek1d(double r0,double k, double t,double sigma,double theta,double T,double *price)
{
  double A,B;

  /*A,B coefficient*/  
  B=(1./k)*(1.-exp(-k*(T-t)));
  A=exp((theta-SQR(sigma)/(2.*SQR(k)))*(B-T+t)-(SQR(sigma)/(4.*k))*SQR(B));
  
  /*Price*/
  *price=A*exp(-B*r0);

  return OK;
}

int CALC(CF_ZCBond)(void *Opt,void *Mod,PricingMethod *Met)
{
  TYPEOPT* ptOpt=(TYPEOPT*)Opt;
  TYPEMOD* ptMod=(TYPEMOD*)Mod;


  return zcb_vasicek1d(ptMod->r0.Val.V_PDOUBLE,ptMod->k.Val.V_DOUBLE,ptMod->T.Val.V_DATE,ptMod->Sigma.Val.V_PDOUBLE,ptMod->theta.Val.V_PDOUBLE,ptOpt->BMaturity.Val.V_DATE,&(Met->Res[0].Val.V_DOUBLE));
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
  "CF_Vasicek1d_ZCBond",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(CF_ZCBond),
  {{"Price",DOUBLE,{100},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(CF_ZCBond),
  CHK_ok,
  MET(Init)
} ;


