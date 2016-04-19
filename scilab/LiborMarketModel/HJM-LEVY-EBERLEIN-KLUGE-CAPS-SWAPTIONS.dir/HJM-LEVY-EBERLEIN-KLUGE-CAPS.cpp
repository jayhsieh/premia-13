/*---------------------------------------------------------------*/
/*     Exact pricing formulae for caps in Lévy                   */
/*                term strucure model                            */
/*                                                               */
/*---------------------------------------------------------------*/
/*  Audrey Drif, Premia 2006                                     */
/*---------------------------------------------------------------*/


#include <iostream>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include "HJM-LEVY-EBERLEIN-KLUGE-CAPS-SWAPTIONS.h" 
using namespace std;
typedef complex<double> dcomplex;

int main()
{
  typedef complex<double> dcomplex;
  dcomplex I(0,1);
  double sig,b,c,lamb,alpha,a,mu,var,bt,K,bT,gam,x;
  double t,T;
  int m,n,nn,A,N,mm,mmm,i;
  dcomplex w,ww,www,wwww,wr;
  //interface
   cout<<"\n";
  cout<<"Pricing Option on ZCB using HJM Levy models\n";
  cout<<"\n";
  cout<<"Choose Volatility Structure:\n";
  cout<<"\n";
  cout<<"(0) Ho-Lee volatility structure : \n";
  cout<<"(1) Vasicek volatility structure : \n";
  cout<<"(2) Moraleda-Vorst volatility structure:\n";
   cin>>mm;
  cout<<"Choose Levy measure:\n";
  cout<<"\n";
  cout<<"(0) drift + Brownian motion:\n";
  cout<<"(1) drift + Brownian motion + compound poisson :\n";
  cout<<"(2) drift + Brownian motion + Levy with Levy measure a*exp(-lambda*x)/x,x>0:\n";
  cout<<"(3) drift + Brownian motion + Levy with Levy measure a*exp(-lambda*x)/x^{2},x>0:\n";
  cout<<"(4) drift + Brownian motion + Levy with Levy measure a*exp(-lambda*x)/x^{alpha+1},x>0 :\n";
  cin>>mmm;

  cout<<"enter the strike K :";
  cin>>K;
  cout<<"enter the drift b :";
  cin>>b;
  cout<<"enter the brownian coefficient c (>0) example 0.5:";
  cin>>c;
  
  bondprice* bb;
  Ho_Lee* H;
  Vasicek* V;
  Moraleda* M;
  Browmian* B;
  Poisson* P;
  Gamma* G;
  Alpha1* A1;
  Alphag* Ag;
  IntegraleGeneric *In;
  IntegraleAnalytiqueBrownien* Inn;
  IntegraleAnalytiqueGamma* Innn;
  IntegraleAnalytiqueAlpha1* I4;
  IntegraleAnalytiqueAlphag* I5;
  cap* cc;
  bb=new bondprice();
  lvgamma* bg;
  lalpha1* ba1;
  lalphag* bag; 
  
  
  switch(mm)
    {
    case 0://case Ho-Lee: sigma(s,T)=tilde{sigma}
      cout<<"enter the parameter tilde{sigma} of Ho-Lee volatility structure (example 0.15) :";
      cin>>sig;
      H=new Ho_Lee(sig);
      switch(mmm)
	{
	case 0://case brownian motion
	  B=new Browmian(b,c);
	  In=new IntegraleGeneric(25,H,B);
	  Inn=new IntegraleAnalytiqueBrownien(b,c,sig,H);
	  cc=new cap(3000,K,-1.5,bb,Inn);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 1://compound poisson
	  cout<<"enter the parameter lambda of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter mu of the Lévy density :";
	  cin>>mu;
	  cout<<"enter the parameter var of the Lévy density :";
	  cin>>var;
	  P=new Poisson(b,c,lamb,mu,var);
	  www=P->theta(1.);
	  In=new IntegraleGeneric(25,H,P);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 2://Levy measure a*exp(-lambda*x)/x,x>0
	  cout<<"enter the parameter lambda of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  bg=new lvgamma(lamb,a);
	  b=b+bg->intlevy();
	  G=new Gamma(b,c,lamb,a);
	  In=new IntegraleGeneric(25,H,G);
	  Innn=new IntegraleAnalytiqueGamma(b,c,a,lamb,sig,H);
	  cc=new cap(3000,K,-1.5,bb,Innn);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	      break;
	case 3://Levy measure a*exp(-lambda*x)/x^{2},x>0
	  cout<<"enter the parameter lambda of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  ba1=new lalpha1(lamb,a);
	  b=b+ba1->intlevy();
	  A1=new Alpha1(b,c,lamb,a);
	  In=new IntegraleGeneric(25,H,A1);
	  I4=new IntegraleAnalytiqueAlpha1(b,c,a,lamb,sig,H);
	  cc=new cap(3000,K,-1.5,bb,I4);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 4://Levy measure a*exp(-lambda*x)/x^{alpha+1},x>0
	  cout<<"enter the parameter lambda of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  cout<<"enter the parameter alpha of the Lévy density :";
	  cin>>alpha;
	  bag=new lalphag(lamb,a,alpha);
	  b=b+bag->intlevy();
	  Ag=new Alphag(b,c,lamb,a,alpha);
	  In=new IntegraleGeneric(25,H,Ag);
	  I5=new IntegraleAnalytiqueAlphag(b,c,a,lamb,sig,alpha,H);
	  cc=new cap(3000,K,-1.5,bb,I5);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	}
      break;
    case 1://case Vasicek: sigma(s,T)=tilde{sigma}*exp(-x*(T-s))
      cout<<"enter the parameter tilde{sigma} of Vasicek volatility structure (example 0.15) :";
      cin>>sig;
      cout<<"enter the parameter x of Vasicek volatility structure(!=0 example 0.02):";
      cin>>x;
      V=new Vasicek(x,sig);
      switch(mmm)
	{
	case 0://case brownian motion
	  B=new Browmian(b,c);
	  In=new IntegraleGeneric(25,V,B);
	  cout<<w<<endl;
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;  
	  break;
	case 1://compound poisson
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter mu of the Lévy density :";
	  cin>>mu;
	  cout<<"give the parameter var of the Lévy density :";
	  cin>>var;
	  P=new Poisson(b,c,lamb,mu,var);
	  In=new IntegraleGeneric(25,V,P);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 2://Levy measure a*exp(-lambda*x)/x,x>0
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  bg=new lvgamma(lamb,a);
	  b=b+bg->intlevy();
	  G=new Gamma(b,c,lamb,a);
	  In=new IntegraleGeneric(25,V,G);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 3://Levy measure a*exp(-lambda*x)/x^{2},x>0
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  ba1=new lalpha1(lamb,a);
	  b=b+ba1->intlevy();
	  A1=new Alpha1(b,c,lamb,a);
	  In=new IntegraleGeneric(25,V,A1);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 4://Levy measure a*exp(-lambda*x)/x^{alpha+1},x>0
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  cout<<"enter the parameter alpha of the Lévy density :";
	  cin>>alpha;
	  bag=new lalphag(lamb,a,alpha);
	  b=b+bag->intlevy();
	  Ag=new Alphag(b,c,lamb,a,alpha);
	  In=new IntegraleGeneric(25,V,Ag);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	}
      break;
    case 2://case Moraleda: sigma(s,T)=tilde{sigma}*exp(-x(T-s))*(1+gamma*T)/(1+gamma*s)
      cout<<"enter the parameter tilde{sigma} of Moraleda-Vorst volatility structure (example 0.15) :";
      cin>>sig;
      cout<<"enter the parameter x of Moraleda-Vorst volatility structure :";
      cin>>x;
      cout<<"enter the parameter gamma>0 of Moraleda-Vorst volatility structure :";
      cin>>gam;
      M=new Moraleda(sig,x,gam);
      switch(mmm)
	{
	case 0://case brownian motion
	  B=new Browmian(b,c);
	  In=new IntegraleGeneric(25,M,B);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 1://compound poisson
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter mu of the Lévy density :";
	  cin>>mu;
	  cout<<"enter the parameter var of the Lévy density :";
	  cin>>var;
	  P=new Poisson(b,c,lamb,mu,var);
	  In=new IntegraleGeneric(25,M,P);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 2://Levy measure a*exp(-lambda*x)/x,x>0
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  bg=new lvgamma(lamb,a);
	  b=b+bg->intlevy();
	  G=new Gamma(b,c,lamb,a);
	  In=new IntegraleGeneric(25,M,G);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 3://Levy measure a*exp(-lambda*x)/x^{2},x>0
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  ba1=new lalpha1(lamb,a);
	  b=b+ba1->intlevy();
	  A1=new Alpha1(b,c,lamb,a);
	  In=new IntegraleGeneric(25,M,A1);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	case 4://Levy measure a*exp(-lambda*x)/x^{alpha+1},x>0
	  cout<<"enter the parameter lambda>0 of the Lévy density :";
	  cin>>lamb;
	  cout<<"enter the parameter a of the Lévy density :";
	  cin>>a;
	  cout<<"enter the parameter alpha of the Lévy density :";
	  cin>>alpha;
	  bag=new lalphag(lamb,a,alpha);
	  b=b+bag->intlevy();
	  Ag=new Alphag(b,c,lamb,a,alpha);
	  In=new IntegraleGeneric(25,M,Ag);
	  cc=new cap(3000,K,-1.5,bb,In);
	  ww=cc->pxcap();
	  cout<<real(ww)<<endl;
	  break;
	  
	    }
      break;
    }
}

