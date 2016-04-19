
#include  "stein1d_std.h"


static double d1(double x,double t,double s,double K,double r,double T){
	double d=(log(x/K)+(r+s*s/2)*(T-t))/(s*sqrt(T-t));
	return d;
}
static double H(double t, double x, double v,double K,double r,double T){
	double a,d,HH; 
	a=d1(x,t,v,K,r,T)*d1(x,t,v,K,r,T);
	d=v*sqrt((T-t)*2*M_PI);
	HH=exp(-a/2)/d*x*(1-d1(x,t,v,K,r,T)/v/sqrt(T));
	return HH;
}

static double diffH(double v,double T, double S, double K, double r){
	return(-0.5/pow(v,3)*pow(2,0.5)/pow(M_PI,0.5)/pow(T,0.5)*(log(S/K)+(r+1./2*v*v)*T)/T/S*exp(-1./2*pow((log(S/K)+(r+1./2*v*v)*T),2)/pow(v,2)/T)*(1-1./2*(log(S/K)+(r+1./2*v*v)*T)/pow(v,2)/pow(T,0.5)*pow(2,0.5)/sqrt(M_PI)/pow(T,0.5))-1./2/pow(v,3)/M_PI/pow(T,3./2)*exp(-1./2*pow(log(S/K)+(r+1./2*v*v)*T,2)/pow(v,2)/T)/S);
}

int ApAlosStein(double S,NumFunc_1  *p, double T, double r, double divid, double v0,double kappa,double theta,double sigma,double rho,double *ptprice, double *ptdelta)
{
  int flag_call;
  double K,prix,delta,price_bs,delta_bs;
  double I1,I2,I,Ivol,vol,d,lambda;
  double sigma0;

  sigma0=sqrt(v0);
  sigma0=v0;
  K=p->Par[0].Val.V_PDOUBLE;
  
  if ((p->Compute)==&Call)
    flag_call=1;
  else
    flag_call=0;;
  
  // Quantity introduce in the paper to simplify the calculation
  lambda=sigma/sqrt(kappa);
  I2=1./4*(-9*exp(2*T*kappa)*pow(theta,2)-2*pow(theta,2)*T*kappa-4*theta*sigma0*T*kappa*exp(T*kappa)+4*pow(theta,2)*T*kappa*exp(2*T*kappa)+4*pow(theta,2)*T*kappa*exp(T*kappa)-8*theta*sigma0*exp(T*kappa)+12*pow(theta,2)*exp(T*kappa)-2*pow(sigma0,2)*T*kappa-pow(sigma0,2)+4*theta*sigma0*T*kappa+4*theta*sigma0-3*pow(theta,2)+pow(sigma0,2)*exp(2*T*kappa)+4*exp(2*T*kappa)*theta*sigma0)/pow(kappa,2)*exp(
    -2*T*kappa);
	//Calculation of the quantity denote by I1 in the paper
	I1=lambda*lambda/2*(-1./4*(2*T*exp(-2*T*kappa)*kappa+exp(-2*T*kappa)-1)/pow(kappa,2));
	I=I1+I2;
	//Calculation of the quantity denote by sigma0* in the paper
	Ivol=1./4*pow(lambda,2)*(2*T*kappa+exp(-2*T*kappa)-1)/kappa+1./2*(-pow(sigma0,2)+2*theta*sigma0-theta*theta-4*theta*sigma0*exp(T*kappa)+4*theta*theta*exp(T*kappa)+2*theta*theta*T*kappa*exp(2*T*kappa)+sigma0*sigma0*exp(2*T*kappa)+2*exp(2*T*kappa)*theta*sigma0-3*exp(2*T*kappa)*theta*theta)*exp(-2*T*kappa)/kappa;
	vol=sqrt(1/T*Ivol);
	d=diffH(vol,T,S,K,r);
 
  if(flag_call==1){
    pnl_cf_call_bs(S,K,T,r,divid,vol,&price_bs,&delta_bs);
	prix=price_bs+lambda*sqrt(kappa)*rho*H(0,S,vol,K,r,T)*I;
	delta=delta_bs+lambda*sqrt(kappa)*rho*I*d;
  }
  else{
    pnl_cf_put_bs(S,K,T,r,divid,vol,&price_bs,&delta_bs);
	prix=price_bs+lambda*sqrt(kappa)*rho*H(0,S,vol,K,r,T)*I;
	delta=delta_bs+lambda*sqrt(kappa)*rho*I*d;
  }
  
  /* Price*/
  *ptprice=prix;
    
  /* Delta */
  *ptdelta=delta;

  return OK;
}

int CALC(AP_Alos_Stein)(void *Opt, void *Mod, PricingMethod *Met)
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

      return ApAlosStein(ptMod->S0.Val.V_PDOUBLE,
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

static int CHK_OPT(AP_Alos_Stein)(void *Opt, void *Mod)
{
  if ((strcmp( ((Option*)Opt)->Name,"CallEuro")==0)
      ||(strcmp( ((Option*)Opt)->Name,"PutEuro")==0))

    return OK;
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

PricingMethod MET(AP_Alos_Stein)=
{
  "AP_Alos_Stein",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_Alos_Stein),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_Alos_Stein),
  CHK_ok,
  MET(Init)
};

