
#include  "hes1d_std.h"

/////////////////////////////////////////////////////////////////

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

int ApAlosHeston(double S,NumFunc_1  *p, double T, double r, double divid, double v0,double kappa,double theta,double sigma,double rho,double *ptprice, double *ptdelta)
{
  int flag_call;
  double K,prix,delta,price_bs,delta_bs;
  double v0et,I,d;
  
  K=p->Par[0].Val.V_PDOUBLE;
  
  if ((p->Compute)==&Call)
    flag_call=1;
  else
    flag_call=0;;

  //Calculation of the quantity denote by v0* in the paper
  v0et=sqrt(theta+ 1/(kappa*T)*(v0-theta)*(1-exp(-kappa*T)));

  // Calculation of the quantity denote by I in the paper
  d=diffH(v0et,T,S,K,r);
  I=sigma/kappa/kappa*(theta*(kappa*T-2)+v0+exp(-kappa*T)*(kappa*T*(theta-v0)+2*theta-v0));
 
  if(flag_call==1){
        pnl_cf_call_bs(S,K,T,r,divid,v0et,&price_bs,&delta_bs);
	prix=price_bs+rho/2.*H(0,S,v0et,K,r,T)*I;
	delta=delta_bs+rho/2*I*d;
  }
  else{
        pnl_cf_put_bs(S,K,T,r,divid,v0et,&price_bs,&delta_bs);
	prix=price_bs+rho/2*H(0,S,v0et,K,r,T)*I;
	delta=delta_bs+rho/2*I*d;
  }
  
  /* Price*/
  *ptprice=prix;
    
  /* Delta */
  *ptdelta=delta;

  return OK;
}

int CALC(AP_Alos_Heston)(void *Opt, void *Mod, PricingMethod *Met)
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

      return ApAlosHeston(ptMod->S0.Val.V_PDOUBLE,
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

static int CHK_OPT(AP_Alos_Heston)(void *Opt, void *Mod)
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

PricingMethod MET(AP_Alos_Heston)=
{
  "AP_Alos_Heston",
  {{" ",PREMIA_NULLTYPE,{0},FORBID}},
  CALC(AP_Alos_Heston),
  {{"Price",DOUBLE,{100},FORBID},
   {"Delta",DOUBLE,{100},FORBID} ,
   {" ",PREMIA_NULLTYPE,{0},FORBID}},
  CHK_OPT(AP_Alos_Heston),
  CHK_ok,
  MET(Init)
};

