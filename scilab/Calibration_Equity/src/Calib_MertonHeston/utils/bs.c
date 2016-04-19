#include "bs.h"
//
// ===================================================================
// Approximation de la fonction de repartition de la loi Gaussienne
// precise a 10^-7 pres.
// ===================================================================
static double MyN(double x)
{
  double p,b1,b2,b3,b4,b5,t1,t2,t3,t4,t5,e,nx;
  //
 if(x <0) {return  1. - MyN(-x);}
 //
 p  =  0.2316419;
 b1 =  0.319381530;
 b2 = -0.356563782;
 b3 =  1.781477937;
 b4 = -1.821255978;
 b5 =  1.330274429;
 //
 t1 = 1./(1.+p*x);
 t2 = t1*t1;
 t3 = t2*t1;
 t4 = t2*t2;
 t5 = t4*t1;
 e  = b1*t1+b2*t2+b3*t3+b4*t4+b5*t5;
 //
 nx = 1 - (1./sqrt(2.*PI)) * exp(-x*x/2.) * e;
 //
 return nx;
 // 
}
// ===================================================================
//
// ===================================================================
// Derivee de la fonction de repartition de la loi Gaussienne :
// e^{-x^2/2}/sqrt(2.*PI) tout simplement !
// ===================================================================
static inline  double dN(double x)
{
   return (1./sqrt(2.*PI)) * exp(-x*x/2.);
}
// ===================================================================
//
// ===================================================================
// fomule de BS pour un Call Euro
// Copie de Premia : Src/mod/bs1d/bs1d_std/cf_call.c 
// ===================================================================
int MyCall_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta){
  double sigmasqrt,d1,d2,delta;
  
  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt/2.;
  d2=d1-sigmasqrt;
  delta=exp(-divid*t)*N(d1);
  
  /*Price*/
  *ptprice=s*delta-exp(-r*t)*k*N(d2);
  
  /*Delta*/
  *ptdelta=delta;
  
  return OK;
}
// ===================================================================
// fomule de BS pour un Put Euro
// Copie de Premia : Src/mod/bs1d/bs1d_std/cf_call.c 
// ===================================================================
int MyPut_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta)
{
  double sigmasqrt,d1,d2,delta;
  
  sigmasqrt=sigma*sqrt(t);
  d1=(log(s/k)+(r-divid)*t)/sigmasqrt+sigmasqrt/2.;
  d2=d1-sigmasqrt;
  delta=-exp(-divid*t)*N(-d1);
  
  /*Price*/
  *ptprice=exp(-r*t)*k*N(-d2)+delta*s;
  
  /*Delta*/
  *ptdelta=delta;
  
  return OK;
}
//
// ===================================================================
// d(formule de BS pour un Call Euro) / d Sigma
// ===================================================================
int dCall_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice){
  double sigmasqrt,d1,d2;
  double sqrtt,A,dd1,dd2;
  //
  sqrtt     = sqrt(t);
  sigmasqrt = sigma*sqrt(t);
  A         = (log(s/k)+(r-divid)*t)/sqrtt;
  //
  d1  = A /sigma + sigmasqrt/2.;
  d2  = d1-sigmasqrt;
  dd1 = - A /(sigma*sigma) + 0.5*sqrtt;
  dd2 = dd1 - sqrtt;
  //
  *dprice = s*exp(-divid*t)*dN(d1)*dd1 -exp(-r*t)*k*dN(d2)*dd2;
  //
  return OK;
}
// ===================================================================
// d(formule de BS pour un Put Euro) / d Sigma
// ===================================================================
int dPut_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice)
{
  double sigmasqrt,d1,d2;
  double sqrtt,A,dd1,dd2;
  //
  sqrtt     = sqrt(t);
  sigmasqrt = sigma*sqrt(t);
  A         = (log(s/k)+(r-divid)*t)/sqrtt;
  //
  d1  = A /sigma + sigmasqrt/2.;
  d2  = d1-sigmasqrt;
  dd1 = - A /(sigma*sigma) + 0.5*sqrtt;
  dd2 = dd1 - sqrtt;
  //
  *dprice =-exp(-r*t)*k*dN(-d2)*dd2 + s*exp(-divid*t)*dN(-d1)*dd1;
  //
  return OK;
}
//
// ===================================================================
// Fonction pour calculer le Sigma Implicite dans la Call ou le Put
//  Europeen avec BS et de prix cbs. Methode utlisee : dichotomie.
// ===================================================================
double SigmaImplicite(double eps,double a,double b,int TypeOpt,double St,double T,double K,
			double r,double divid,double cbs,int *nbiters)
{
  //
  int n=0,ibidon;
  double fa,fb,c,fc,sol,prix,delta;
  //
  // Cas du Call
  if(TypeOpt == 1)
	{
	  ibidon = MyCall_BlackScholes_73(St,K,T,r,divid,a,&prix,&delta);
	  fa = prix - cbs;
	  ibidon = MyCall_BlackScholes_73(St,K,T,r,divid,b,&prix,&delta);
	  fb = prix - cbs;
	}
  // Cas du Put
  else {
	ibidon = MyPut_BlackScholes_73(St,K,T,r,divid,a,&prix,&delta);
	fa = prix - cbs;
	ibidon = MyPut_BlackScholes_73(St,K,T,r,divid,b,&prix,&delta);
	fb = prix - cbs;
  }
  //
  // Si pas de solution dans l'intervalle [a,b]
  if (fa*fb > 0.)
	{
	  // On indique que pas de sol trouvee
	  *nbiters = -1;
	  // On renvoit la borne qui donne la prix le plus proche.
	  if(fabs(fa-cbs) < fabs(fb-cbs) )
		{
		  sol =a;
		} else {
		  sol =b;
		}
	  return sol;
	}
  //
  while (fabs(b-a) > eps) 
	{
	  //
	  n++;  
	  c  = (a+b)/2.;
	  // Cas du Call
	  if(TypeOpt == 1)
		{
		  ibidon = MyCall_BlackScholes_73(St,K,T,r,divid,c,&prix,&delta);
		  fc = prix - cbs;
		}
	  // Cas du Put
	  else {
		ibidon = MyPut_BlackScholes_73(St,K,T,r,divid,c,&prix,&delta);
		fc = prix - cbs;
	  }
	  //
	  //
	  if(fc == 0.) {
		sol = c;
		*nbiters = n;
		return sol;
	  }
	  //
	  if (fc*fa < 0. )
		{b = c; fb = fc;}
	  else
		{a = c; fa = fc;};
	}
  //
  *nbiters = n;
  sol = (a+b)/2.;
  //
  return sol;
  //
}
// ===================================================================
