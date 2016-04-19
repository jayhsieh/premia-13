#include  "vasicek1d_stdi.h"

static double A,B;

/*Zero Coupon Bond*/
static double zcb_vasicek1d(double theta, double r,double k,double sigma,double ti,double Ti)
{
  /*A,B coefficient*/  
  B=(1./k)*(1.-exp(-k*(Ti-ti)));
  A=exp((theta-SQR(sigma)/(2.*SQR(k)))*(B-Ti+ti)-(SQR(sigma)/(4.*k))*SQR(B));

  return A*exp(-B*r);
}

/*Call Option*/
static int zbc_vasicek1d(double r, double k,double t, double sigma,double theta, double S, double T,NumFunc_1 *p,double *price,double *delta)
{
  double PtS,PtT;
  double d1,d2,sigma_p,K;

  K=p->Par[0].Val.V_DOUBLE;
  PtT=zcb_vasicek1d(theta,r,k,sigma,t,T);
  PtS=zcb_vasicek1d(theta,r,k,sigma,t,S);
  sigma_p=sigma*sqrt((1.-exp(-2.*k*(T-t)))/(2*k))*(1./k)*(1.-exp(-k*(S-T)));
  d1=1./(sigma_p)*log(PtS/(PtT*K))+0.5*sigma_p;
  d2=d1-sigma_p;

  /*Price*/
  *price=PtS*cdf_nor(d1)-K*PtT*cdf_nor(d2);

  /*Delta*/
  *delta=cdf_nor(d1);

  return OK;
}

int CALC(CF_ZCCallBondEuro)(void *Opt,void *Mod,PricingMethod *Met)
{
  TYPEOPT* ptOpt=(TYPEOPT*)Opt;
  TYPEMOD* ptMod=(TYPEMOD*)Mod;


  return zbc_vasicek1d(ptMod->r0.Val.V_PDOUBLE,ptMod->k.Val.V_DOUBLE,ptMod->T.Val.V_DATE,ptMod->Sigma.Val.V_PDOUBLE,
		       ptMod->theta.Val.V_PDOUBLE,ptOpt->BMaturity.Val.V_DATE,ptOpt->OMaturity.Val.V_DATE,ptOpt->PayOff.Val.V_NUMFUNC_1,
		       &(Met->Res[0].Val.V_DOUBLE),&(Met->Res[1].Val.V_DOUBLE));
}


static int CHK_OPT(CF_ZCCallBondEuro)(void *Opt, void *Mod)
{
  return strcmp( ((Option*)Opt)->Name,"ZeroCouponCallBondEuro");
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(CF_ZCCallBondEuro)=
{
  "CF_Vasicek1d_ZBCallEuro",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(CF_ZCCallBondEuro),
  {{"Price",DOUBLE,{100},FORBID},{"Delta",DOUBLE,{100},FORBID} ,{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(CF_ZCCallBondEuro),
  CHK_ok,
  MET(Init)
} ;
