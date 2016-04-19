#include  "vasicek1d_stdi.h"

static int nb_payement;
static double A,B;

/*Zero Coupon Bond*/
static double zcb_vasicek1d(double theta, double r,double k,double sigma,double ti,double Ti)
{
  
  B=(1./k)*(1.-exp(-k*(Ti-ti)));
  A=exp((theta-SQR(sigma)/(2.*SQR(k)))*(B-Ti+ti)-(SQR(sigma)/(4.*k))*SQR(B));

  return A*exp(-B*r);
}

/*Call Option on Zero Coupon Bond*/
static double zbc_vasicek1d(double t,double T,double S,double r, double k, double theta,double sigma,double K,double periodicity)
{
  double PtS,PtT;
  double d1,d2,sigma_p;
  double new_K;
  
  new_K=1./(1.+K*periodicity);
 
  PtT=zcb_vasicek1d(theta,r,k,sigma,t,T);
  PtS=zcb_vasicek1d(theta,r,k,sigma,t,S);
  sigma_p=sigma*sqrt((1.-exp(-2.*k*(T-t)))/(2*k))*(1./k)*(1.-exp(-k*(S-T)));
  d1=1./(sigma_p)*log(PtS/(PtT*new_K))+0.5*sigma_p;
  d2=d1-sigma_p;

  return PtS*cdf_nor(d1)-new_K*PtT*cdf_nor(d2);
}

/*Floor*/
static int floor_vasicek1d(double r,double k, double date,double sigma,double theta,double Nominal,double K,double periodicity,double first_payement,double contract_maturity,double *price)
{
  double sum,tim,tip;
  int i;

  nb_payement=(int)((contract_maturity-first_payement)/periodicity);

  /*Floor=Portfolio of zero-bond Call options*/
  sum=0.;
  for(i=0;i<nb_payement;i++)
    {
      tim=first_payement+(double)i*periodicity;
      tip=tim+periodicity;
      sum+=(1.+K*periodicity)*zbc_vasicek1d(date,tim,tip,r,k,theta,sigma,K,periodicity);
    }
  
  /*Price*/
  *price=Nominal*sum;

  return OK;
}

int CALC(CF_Floor)(void *Opt,void *Mod,PricingMethod *Met)
{
  TYPEOPT* ptOpt=(TYPEOPT*)Opt;
  TYPEMOD* ptMod=(TYPEMOD*)Mod;

  return floor_vasicek1d(ptMod->r0.Val.V_PDOUBLE,ptMod->k.Val.V_DOUBLE,ptMod->T.Val.V_DATE,ptMod->Sigma.Val.V_PDOUBLE,ptMod->theta.Val.V_PDOUBLE,ptOpt->Nominal.Val.V_PDOUBLE,ptOpt->FixedRate.Val.V_PDOUBLE,ptOpt->ResetPeriod.Val.V_DATE,ptOpt->FirstResetDate.Val.V_DATE,ptOpt->BMaturity.Val.V_DATE,&(Met->Res[0].Val.V_DOUBLE));
}


static int CHK_OPT(CF_Floor)(void *Opt, void *Mod)
{
  return strcmp( ((Option*)Opt)->Name,"Floor");
}

static int MET(Init)(PricingMethod *Met,Option *Opt)
{
  if ( Met->init == 0)
    {
      Met->init=1;
    }

  return OK;
}

PricingMethod MET(CF_Floor)=
{
  "CF_Vasicek1d_Floor",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(CF_Floor),
  {{"Price",DOUBLE,{100},FORBID},{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(CF_Floor),
  CHK_ok,
  MET(Init)
} ;


