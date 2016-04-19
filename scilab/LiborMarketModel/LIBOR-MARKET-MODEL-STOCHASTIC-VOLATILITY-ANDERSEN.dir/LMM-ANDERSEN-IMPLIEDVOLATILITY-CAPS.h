/*---------------------------------------------------------------------------*/
/*----------Extended Libor Market Models with stochastic Volatility----------*/
/*---------------------------------------------------------------------------*/
/*----------------------Audrey Drif, Premia 2006-----------------------------*/
/*---------------------------------------------------------------------------*/



#include <iostream>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;
typedef complex<double> dcomplex;
const double PI = 3.141592653589793;

class Varphi{//calculate int_ {X}^{F} varphi(u)^{-1} du
 public:
  virtual double varphi(double F,double X,double alpha)=0;
};

class Varphi1 : public Varphi{//for the case varphi(u)=u^{alpha}
 public:
  Varphi1(){}
  double varphi(double F,double X,double alpha);
}; 

class Varphi2 : public Varphi{//for the case varphi(u)=u
 public:
  Varphi2(){}
  double varphi(double F,double X,double alpha);
}; 

class Sigma0{//calculate Omega0(F)
 public:
  Varphi* v;
  Sigma0(Varphi* _v) : v(_v) {}
  double sigma0(double F, double X,double alpha);
};

class Sigma1{//calculate Omega1(F)
 public:
  Sigma0 *s0;
  Varphi* v;
  Sigma1(Varphi* _v,Sigma0* _s0) : v(_v),s0(_s0){}
  double sigma1(double F, double X,double alpha);
};

class Volconst{//calculate Omega0(F)*lamb+ Omega1(F)*lamb^{3}*T=Omega/sqrt(T)
 public:
  Sigma0 *s0;
  Sigma1 *s1;
  Volconst(Sigma0* _s0,Sigma1* _s1) : s0(_s0),s1(_s1)  {}
  double volconst(double T,double lamb,double F,double X,double alpha);
};

class C1{//calculate T^{-1}*mu_{U}(0,1)
 public:
  C1(){}
  double c1(double lamb,double T);
};

class L12{//calculate l_{1,2}(0,1)
 public:
  virtual double l12(double K, double lamb,double T)=0;
};

class L12K : public L12{//for the case K>0 
 public:
  L12K(){}
  double l12(double K, double lamb,double T);
};

class L120 : public L12{//for the case K=0
 public:
  L120(){}
  double  l12(double K, double lamb,double T);
};

class L22{//calculate l_{2,2}(0,1)
 public:
  virtual double l22(double K, double lamb,double T)=0;
};

class L22K : public L22{//for the case K>0 
 public:
  L22K(){}
  double l22(double K, double lamb,double T);
};

class L220 : public L22{//in the case K=0
 public:
  L220(){}
  double l22(double K, double lamb,double T);
};

class L23{//calculate l_{2,3}(0,1)
 public:
  virtual double l23(double K, double lamb,double T)=0;
};

class L23K : public L23{//in the case K>0 
 public:
  L23K(){}
  double l23(double K, double lamb,double T);
};

class L230 : public L23{//in the case K=0
 public:
  L230(){}
  double l23(double K, double lamb,double T);
};

class Dsigma1{//the first derivative of Omega with respect to lamb
 public:
  Sigma0 *s0;
  Sigma1 *s1;
Dsigma1(Sigma0* _s0,Sigma1* _s1) : s0(_s0),s1(_s1)  {}
 double dsigma1(double T,double F, double X,double alpha,double lamb);
};

class Dsigma2{//the second derivative of Omega with respect to lamb
 public:
  Sigma0 *s0;
  Sigma1 *s1;
  Dsigma2(Sigma0* _s0,Sigma1* _s1) : s0(_s0),s1(_s1)  {}
  double dsigma2(double T,double F, double X,double alpha,double lamb);
};

class Dsigma3{//the third derivative of Omega with respect to lamb
 public:
  Sigma0 *s0;
  Sigma1 *s1;
  Dsigma3(Sigma0* _s0,Sigma1* _s1) : s0(_s0),s1(_s1)  {}
  double dsigma3(double T,double F, double X,double alpha,double lamb);
};

class Dsigma4{//the fourth derivative of Omega with respect to lamb
 public:
  Sigma0 *s0;
  Sigma1 *s1;
  Dsigma4(Sigma0* _s0,Sigma1* _s1) : s0(_s0),s1(_s1)  {}
  double dsigma4(double T,double F, double X,double alpha,double lamb);
};

class Sigma10{//the first derivative of Omega with respect to lamb divided by Omega
 public:
  Dsigma1 *Ds1;
  Volconst *vol;
  Sigma10(Volconst* _vol,Dsigma1* _Ds1) : vol(_vol),Ds1(_Ds1) {}
  double sigma10(double T,double F, double X,double alpha,double lamb);  
};

class Sigma21{//the second derivative of Omega with respect to lamb divided by the first derivative of Omega with respect to lamb 
 public:
  Dsigma1 *Ds1;
  Dsigma2 *Ds2;
  Sigma21( Dsigma1 *_Ds1,Dsigma2* _Ds2) : Ds1(_Ds1),Ds2(_Ds2) {}
  double sigma21(double T,double F, double X,double alpha,double lamb);  
};

class Sigma20{//the second derivative of Omega with respect to lamb divided by Omega
 public:
  Dsigma2 *Ds2;
  Volconst *vol;
  Sigma20(Volconst* _vol,Dsigma2* _Ds2) :  vol(_vol),Ds2(_Ds2) {}
  double sigma20(double T,double F, double X,double alpha,double lamb);  
};

class Sigma30{//the third derivative of Omega with respect to lamb divided by Omega
 public:
  Dsigma3 *Ds3;
  Volconst *vol;
  Sigma30(Volconst* _vol,Dsigma3* _Ds3) :  vol(_vol),Ds3(_Ds3) {}
  double sigma30(double T,double F, double X,double alpha,double lamb);  
};

class Sigma31{//the third derivative of Omega with respect to lamb divided by the first derivative of Omega with respect to lamb
 public:
  Dsigma1 *Ds1;
  Dsigma3 *Ds3;
  Sigma31( Dsigma1 *_Ds1,Dsigma3 *_Ds3) : Ds1(_Ds1),Ds3(_Ds3) {}
  double sigma31(double T,double F, double X,double alpha,double lamb);  
};

class Sigma41{//the fourth derivative of Omega with respect to lamb divided by the first derivative of Omega with respect to lamb 
 public:
  Dsigma1 *Ds1;
  Dsigma4 *Ds4;
  Sigma41( Dsigma1 *_Ds1,Dsigma4 *_Ds4) : Ds1(_Ds1),Ds4(_Ds4) {}
  double sigma41(double T,double F, double X,double alpha,double lamb);  
};

class Alpha0{//calculate the constant alpha0
 public:
  Volconst *vol;
  L12 *l;
  Sigma10 *s10;
  Sigma21 *s21;
Alpha0(Volconst* _vol,L12* _l,Sigma10* _s10,Sigma21* _s21) : vol(_vol),l(_l),s10(_s10),s21(_s21) {}
double alpha0(double T,double F,double X,double alpha,double lamb,double K);
};


class Alpha1{//calculate the constant alpha1
public:
  Volconst *vol;
  L12 *l;
  Sigma10 *s10;
  Alpha1(Volconst* _vol,L12* _l,Sigma10* _s10) : vol(_vol),l(_l),s10(_s10){}
double alpha1(double T,double F,double X,double alpha,double lamb,double K);
};

class Beta0{//calculate the constant beta0
 public:
  Volconst *vol;
  L12 *l;
  Sigma10 *s10;
  Sigma21 *s21;
  Sigma20 *s20;
  Sigma30 *s30;
  Sigma31 *s31;
  Sigma41 *s41;
  L22 *ll;
  L23 *lll;
  Beta0(Volconst* _vol,L12* _l,L22 *_ll,L23 *_lll,Sigma10* _s10,Sigma20* _s20,Sigma21* _s21,Sigma30* _s30,Sigma31* _s31,Sigma41* _s41) : vol(_vol),l(_l),ll(_ll),lll(_lll),s10(_s10),s20(_s20),s21(_s21),s30(_s30),s31(_s31),s41(_s41) {} 
  double beta0(double T,double F,double X,double alpha,double lamb,double K);
};

class Beta1{//calculate the constant beta1
 public:
  Volconst *vol;
  L12 *l;
  Sigma10 *s10;
  Sigma21 *s21;
  Sigma20 *s20;
  Sigma30 *s30;
  Sigma31 *s31;
  L22 *ll;
  L23 *lll;
  Beta1(Volconst* _vol,L12* _l,L22 *_ll,L23 *_lll,Sigma10* _s10,Sigma20* _s20,Sigma21* _s21,Sigma30* _s30,Sigma31* _s31) : vol(_vol),l(_l),ll(_ll),lll(_lll),s10(_s10),s20(_s20),s21(_s21),s30(_s30),s31(_s31){} 
  double beta1(double T,double F,double X,double alpha,double lamb,double K);
};

class Beta2{//calculate the constant beta2
 public:
  Volconst *vol;
  L12 *l;
  L23 *lll;
  Sigma10 *s10;
  Sigma20 *s20;
 Beta2(Volconst* _vol,L12* _l,L23 *_lll,Sigma10* _s10,Sigma20* _s20) : vol(_vol),l(_l),lll(_lll),s10(_s10),s20(_s20){} 
  double beta2(double T,double F,double X,double alpha,double lamb,double K);
};


class Cetoile{//calculate c^{*}(0,1) 
 public:
  virtual double cetoile(double T,double F,double X,double alpha,double lamb,double K)=0;
};

class C2 : public Cetoile{//to order varepsilon^{4}
 public:
  C1 *cc1;
  Alpha0 *a0;
  Alpha1 *a1;
  double eps;//constant varepsilon,
  C2(C1* _cc1,Alpha0* _a0,Alpha1* _a1,double _eps) : cc1(_cc1),a0(_a0),a1(_a1) ,eps(_eps){}
  double cetoile(double T,double F,double X,double alpha,double lamb,double K);
};

class C3 : public Cetoile{//to order varepsilon^{6}
 public:
  C1 *cc1;
  Alpha0 *a0;
  Alpha1 *a1;
  Beta0 *b0;
  Beta1 *b1;
  Beta2 *b2;
  double eps;
  C3(C1* _cc1,Alpha0* _a0,Alpha1* _a1,Beta0 *_b0,Beta1 *_b1,Beta2 *_b2,double _eps) : cc1(_cc1),a0(_a0),a1(_a1),b0(_b0),b1(_b1),b2(_b2),eps(_eps){}
  double cetoile(double T,double F,double X,double alpha,double lamb,double K);
};

class C4 : public Cetoile{//to order varepsilon^{6} and Lambda=0
 public:
  C1 *cc1;
  Alpha0 *a0;
  Alpha1 *a1;
  Beta0 *b0;
  Beta1 *b1;
  Beta2 *b2;
  double eps;
  C4(C1* _cc1,Alpha0* _a0,Alpha1* _a1,Beta0 *_b0,Beta1 *_b1,Beta2 *_b2,double _eps) : cc1(_cc1),a0(_a0),a1(_a1),b0(_b0),b1(_b1),b2(_b2),eps(_eps){}
  double cetoile(double T,double F,double X,double alpha,double lamb,double K);
};

class Volsto{//calculate Black-Scholes implied
 public:
  virtual double volsto(double T,double lamb,double F,double X,double alpha,double K)=0;
};

class Volsto1 : public Volsto{//with C2
 public:
  C2* cc2;
  Sigma0* s0;
  Sigma1 *s1;
  Volsto1(C2* _cc2,Sigma0* _s0,Sigma1* _s1) :cc2(_cc2),s0(_s0),s1(_s1) {}
  double volsto(double T,double lamb,double F,double X,double alpha,double K);
};

class Volsto2 : public Volsto{//with C3
 public:
  C3* cc3;
  Sigma0* s0;
  Sigma1 *s1;
  Volsto2(C3* _cc3,Sigma0* _s0,Sigma1* _s1) :cc3(_cc3),s0(_s0),s1(_s1) {}
  double volsto(double T,double lamb,double F,double X,double alpha,double K);
};

class Volsto3 : public Volsto{//with C4
 public:
  C4* cc4;
  Sigma0* s0;
  Sigma1 *s1;
  Volsto3(C4* _cc4,Sigma0* _s0,Sigma1* _s1) :cc4(_cc4),s0(_s0),s1(_s1) {}
  double volsto(double T,double lamb,double F,double X,double alpha,double K);
};

class bondprice{//calculate bondprice
 public:
  bondprice(){}
  double* bond1(double*,int);/*calculate bondprice with interest rate*/
  double* bond2(double*,double*,int,double*,int);/*calculate bond price by interpolation*/
};

class Caplet{//caplet pricing formulas with Volsto
 public:
  double X;//strike
  bondprice * b;
  Volsto *v;
  Caplet(double _X,bondprice * _b,Volsto* _v) : X(_X),b(_b), v(_v) {}
  double caplet(double lamb,double alpha,double K);
};

class Caplet0{//caplet pricing formulas with Volconst
 public:
  double X;//strike
  bondprice * b;
  Volconst *v;
  Caplet0(double _X,bondprice * _b,Volconst* _v) : X(_X),b(_b), v(_v) {}
  double caplet0(double lamb,double alpha);
};




double Varphi1::varphi(double F, double X, double alpha)
{
  return((pow(F,1.-alpha)-pow(X,1.-alpha))/(-alpha+1.));
}

double Varphi2::varphi(double F, double X, double alpha)
{
  return(log(F)-log(X));
}

double Sigma0::sigma0(double F, double X,double alpha)
{
  return(log(F/X)/(v->varphi(F,X,alpha)));
}

double Sigma1::sigma1(double F,double X,double alpha)
{
  double sig;
  sig=s0->sigma0(F,X,alpha);
  return((-sig/pow(v->varphi(F,X,alpha),2.))*log(sig*sqrt(F*X/pow(F*X,alpha))));
}

double Volconst::volconst(double T,double lamb,double F,double X,double alpha)
{
  return((s0->sigma0(F,X,alpha))*lamb+(s1->sigma1(F,X,alpha))*pow(lamb,3.)*T);
}

double C1::c1(double lamb,double T)
{
  return(pow(lamb,2));
}

double L12K::l12(double K,double lamb,double T)
{
  return((pow(lamb,4)/2.)*(T/pow(K,2)+(1./(2*K)-exp(-2*K*T)/(2*K))/pow(K,2)-(2./pow(K,2))*((1./K)-exp(-K*T)/K)));
}

double L120::l12(double K,double lamb,double T)
{
  return(pow(lamb,4)*pow(T,3)/6.);
}

double L22K::l22(double K,double lamb,double T)
{
  return(((3.*pow(lamb,4))/(16.*pow(K,2)))*(T/(2.*K)-(1.-exp(-2.*K*T))/(4.*pow(K,2))-exp(-2.*K*T)*T/(2.*K)+(1.-exp(-2.*K*T))/(4.*pow(K,2))+(1.-exp(-2*K*T))/pow(K,2)-2.*(1.-exp(-K*T))/pow(K,2)));
}

double L220::l22(double K,double lamb,double T)
{
  return(pow(lamb,4)*pow(T,4)/64.);
}

double L23K::l23(double K,double lamb,double T)
{
  return(((-3.*pow(lamb,6))/(4.*pow(K,3)))*(T/K-(1.-exp(-2*K*T))/(2.*pow(K,2))+2.*exp(-K*T)*T/K-2.*(1.-exp(-K*T))/pow(K,2)-(1.-exp(-K*T))/pow(K,2)+(1.-exp(-3*K*T))/(3.*pow(K,2))-exp(-2.*K*T)*T/K+(1.-exp(-2*K*T))/(2.*pow(K,2))));
}

double L230::l23(double K,double lamb,double T)
{
  return(-pow(lamb,6)*pow(T,5)/20.);
}

double Dsigma1::dsigma1(double T,double F, double X,double alpha,double lamb)
{
  return(pow(T,0.5)*(s0->sigma0(F,X,alpha))/(2.*lamb)+lamb*pow(T,1.5)*(s1->sigma1(F,X,alpha))*3./2.);
}

double Dsigma2::dsigma2(double T,double F, double X,double alpha,double lamb)
{
  return(-pow(T,0.5)*s0->sigma0(F,X,alpha)/(4.*pow(lamb,3))+3.*pow(T,1.5)*s1->sigma1(F,X,alpha)/(4.*lamb));
}

double Dsigma3::dsigma3(double T,double F, double X,double alpha,double lamb)
{
  return(pow(T,0.5)*s0->sigma0(F,X,alpha)*3./(8*pow(lamb,5))-3.*pow(T,1.5)*s1->sigma1(F,X,alpha)/(8*pow(lamb,3)));
}

double Dsigma4::dsigma4(double T,double F, double X,double alpha,double lamb)
{
  return(-15.*pow(T,0.5)*s0->sigma0(F,X,alpha)/(pow(lamb,7)*16)+9.*pow(T,1.5)*s1->sigma1(F,X,alpha)/(pow(lamb,5)*16));
}

double Sigma10::sigma10(double T,double F, double X,double alpha,double lamb)
{
  return(Ds1->dsigma1(T,F,X,alpha,lamb)/(sqrt(T)*(vol->volconst(T,lamb,F,X,alpha))));
}

double  Sigma21::sigma21(double T,double F, double X,double alpha,double lamb)
{
  return(Ds2->dsigma2(T,F,X,alpha,lamb)/Ds1->dsigma1(T,F,X,alpha,lamb));
}

double  Sigma20::sigma20(double T,double F, double X,double alpha,double lamb)
{
  return(Ds2->dsigma2(T,F,X,alpha,lamb)/(sqrt(T)*(vol->volconst(T,lamb,F,X,alpha))));
}

double  Sigma30::sigma30(double T,double F, double X,double alpha,double lamb)
{
  return(Ds3->dsigma3(T,F,X,alpha,lamb)/(sqrt(T)*(vol->volconst(T,lamb,F,X,alpha))));
}

double  Sigma31::sigma31(double T,double F, double X,double alpha,double lamb)
{
  return(Ds3->dsigma3(T,F,X,alpha,lamb)/Ds1->dsigma1(T,F,X,alpha,lamb));
}

double  Sigma41::sigma41(double T,double F, double X,double alpha,double lamb)
{
  return(Ds4->dsigma4(T,F,X,alpha,lamb)/Ds1->dsigma1(T,F,X,alpha,lamb));
}

double Alpha0::alpha0(double T,double F,double X,double alpha,double lamb,double K)
{
  return((s21->sigma21(T,F,X,alpha,lamb)-pow(sqrt(T)*vol->volconst(T,lamb,F,X,alpha),2)*s10->sigma10(T,F,X,alpha,lamb)/4.)*(l->l12(K,lamb,T))/pow(T,2));
}

double Alpha1::alpha1(double T,double F,double X,double alpha,double lamb,double K)
{
  return(l->l12(K,lamb,T)*s10->sigma10(T,F,X,alpha,lamb)/(pow(sqrt(T)*vol->volconst(T,lamb,F,X,alpha),2)*pow(T,2)));
}

double Beta0::beta0(double T,double F,double X,double alpha,double lamb,double K)
{
  double omega10,omega20,omega21,omega30,omega31,omega41,omega;
  double z12,z22,z23;
  omega10=s10->sigma10(T,F,X,alpha,lamb);
  omega20=s20->sigma20(T,F,X,alpha,lamb);
  omega21=s21->sigma21(T,F,X,alpha,lamb);
  omega30=s30->sigma30(T,F,X,alpha,lamb);
  omega31=s31->sigma31(T,F,X,alpha,lamb);
  omega41=s41->sigma41(T,F,X,alpha,lamb);
  omega=sqrt(T)*vol->volconst(T,lamb,F,X,alpha);
  z12=l->l12(K,lamb,T);
  z22=ll->l22(K,lamb,T);
  z23=lll->l23(K,lamb,T);

  return(z22*(omega21-pow(omega,2)*omega10/4.)/pow(T,2)-z23*(omega31-pow(omega21,2)-pow(omega,2)*(omega20+pow(omega10,2))/4.+pow(omega21-pow(omega,2)*omega10/4.,2))/pow(T,3)+pow(z12,2)*(omega41-3.*omega31*omega21+2*pow(omega21,3)-pow(omega,2)*omega30/4.-3.*pow(omega,2)*omega10*omega20/4.+3*(omega21-pow(omega,2)*omega10/4.)*(omega31-pow(omega21,2)-pow(omega,2)*(omega20+pow(omega10,2))/4.))/(2.*pow(T,4)));
}

double Beta1::beta1(double T,double F,double X,double alpha,double lamb,double K)
{
  double omega10,omega20,omega21,omega30,omega31,omega41,omega;
  double z12,z22,z23;
  omega10=s10->sigma10(T,F,X,alpha,lamb);
  omega20=s20->sigma20(T,F,X,alpha,lamb);
  omega21=s21->sigma21(T,F,X,alpha,lamb);
  omega30=s30->sigma30(T,F,X,alpha,lamb);
  omega31=s31->sigma31(T,F,X,alpha,lamb);
  omega=sqrt(T)*vol->volconst(T,lamb,F,X,alpha);
  z12=l->l12(K,lamb,T);
  z22=ll->l22(K,lamb,T);
  z23=lll->l23(K,lamb,T);

  return((z22*omega10/pow(T,2)-z23*(omega20-3.*pow(omega10,2)+2.*omega10*(omega21-pow(omega,2)*omega10/4.))/pow(T,3)+pow(z12,2)*(omega30-9.*omega10*omega20+12.*pow(omega10,3)+3.*omega10*(omega31-pow(omega21,2)-pow(omega,2)*(omega20+pow(omega10,2))/4.)+3.*(omega21-pow(omega,2)*omega10/4.)*(omega20-3.*pow(omega10,2)))/(2.*pow(T,4)))/pow(omega,2));
}

double Beta2::beta2(double T,double F,double X,double alpha,double lamb,double K)
{
  double omega10;
  omega10=s10->sigma10(T,F,X,alpha,lamb);
  return((-lll->l23(K,lamb,T)*pow(omega10,2)/pow(T,3)+3.*pow(l->l12(K,lamb,T),2)*omega10*(s20->sigma20(T,F,X,alpha,lamb)-3.*pow(omega10,2))/(2.*pow(T,4)))/pow(sqrt(T)*vol->volconst(T,lamb,F,X,alpha),4));
}

double C2::cetoile(double T,double F,double X,double alpha,double lamb,double K)
{
  return(cc1->c1(lamb,T)+(a0->alpha0(T,F,X,alpha,lamb,K))*pow(eps,2)+(a1->alpha1(T,F,X,alpha,lamb,K))*pow(eps,2)*pow(log(F/X),2));
}

double C3::cetoile(double T,double F,double X,double alpha,double lamb,double K)
{
  return(cc1->c1(lamb,T)+(a0->alpha0(T,F,X,alpha,lamb,K)*pow(eps,2)+b0->beta0(T,F,X,alpha,lamb,K)*pow(eps,4))+(a1->alpha1(T,F,X,alpha,lamb,K)*pow(eps,2)+b1->beta1(T,F,X,alpha,lamb,K)*pow(eps,4))*pow(log(F/X),2)+b2->beta2(T,F,X,alpha,lamb,K)*pow(eps,4)*pow(log(F/X),4)*exp(-pow(eps,2)*pow(log(F/X),2)));
}

double C4::cetoile(double T,double F,double X,double alpha,double lamb,double K)
{
  return(cc1->c1(lamb,T)+(a0->alpha0(T,F,X,alpha,lamb,K)*pow(eps,2)+b0->beta0(T,F,X,alpha,lamb,K)*pow(eps,4))+(a1->alpha1(T,F,X,alpha,lamb,K)*pow(eps,2)+b1->beta1(T,F,X,alpha,lamb,K)*pow(eps,4))*pow(log(F/X),2)+b2->beta2(T,F,X,alpha,lamb,K)*pow(eps,4)*pow(log(F/X),4));
}

double Volsto1::volsto(double T,double lamb,double F,double X,double alpha,double K)
{
  return(s0->sigma0(F,X,alpha)*sqrt(cc2->cetoile(T,F,X,alpha,lamb,K))+s1->sigma1(F,X,alpha)*pow(cc2->cetoile(T,F,X,alpha,lamb,K),1.5)*T);
}

double Volsto2::volsto(double T,double lamb,double F,double X,double alpha,double K)
{
  return(s0->sigma0(F,X,alpha)*sqrt(cc3->cetoile(T,F,X,alpha,lamb,K))+s1->sigma1(F,X,alpha)*pow(cc3->cetoile(T,F,X,alpha,lamb,K),1.5)*T);
}

double Volsto3::volsto(double T,double lamb,double F,double X,double alpha,double K)
{
  return(s0->sigma0(F,X,alpha)*sqrt(cc4->cetoile(T,F,X,alpha,lamb,K))+s1->sigma1(F,X,alpha)*pow(cc4->cetoile(T,F,X,alpha,lamb,K),1.5)*T);
}

double* bondprice::bond1(double *t,int m)
{//B(0,t)=exp(-r*t)
  double rr;
  cout<<"enter the interest rate r";
  cin>>rr;
  int i;
  double* bt;
  bt=new double[m]; 
  for(i=0;i<m;i++)
    {
      bt[i]=exp(-rr*t[i]);
    }
  return(bt);
}

double* bondprice::bond2(double *tt,double *btt,int n,double *t,int m)
{//linear interpolation
  
  double* bt;
  int *ind;
  int i,j;
  bt=new double[m];
  ind=new int[m];
  int indi;
  for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++)
	{
	  if(tt[j]<=t[i]&tt[j+1]>t[i]) 
	    {
	      ind[i]=j;
	      break;
	    }
	}
      indi=ind[i];
      bt[i]=(btt[indi+1]-btt[indi])/(tt[indi+1]-tt[indi])*t[i]+btt[indi]-tt[indi]*(btt[indi+1]-btt[indi])/(tt[indi+1]-tt[indi]);  
    }
  if(t[m-1]==tt[n-1]) bt[m-1]=btt[n-1];
  return bt;
}

double Caplet::caplet(double lamb,double alpha,double K)
{
  int m,n1,i,n,n1b,n2b;
  double F;
  double delta;
  double *T;
  double sig;
  double d1;
  double d2;
  T=new double[2];
  double *bt;
  bt=new double[2];
  double *tt;
  double *btt;
    static char init[]="initialyield.dat";
  /*Name of the file where to read P(0, T) of the market.*/
    static FILE* Entrees; /*File variable of the code*/
  int etat;
  char ligne[20];
  char* pligne;

  cout<<"enter the date T_{k} :";
  cin>>T[0];
  cout<<"enter the date T_{k+1} :";
  cin>>T[1];
  cout<<T[0]<<endl;

  cout<<"to calculate bond prices verify yields structure given in initialyealds.dat:\n";
    Entrees=fopen(init, "r");
    if(Entrees==NULL){printf("THE FILE cannot be open: verify the path \n");} else {}
    
    i=0;
    pligne=ligne;
    tt=new double[100];
    btt=new double[100];
    while(1)
      {
	pligne=fgets(ligne, sizeof(ligne), Entrees);
        if(pligne==NULL) break;
        else
          {
            sscanf(ligne, "%lf t=%lf", &(btt[i]), &(tt[i]));
            i++;
          }
      }
    n1=i;
    etat=fclose( Entrees); 

 //Linear Interpolation
  bt=b->bond2(tt,btt,n1,T,2);
   
  
  delta=T[1]-T[0];
  F=(bt[0]/bt[1]-1.)/delta;
  //cout<<F<<endl;
  //cout<<T[0]<<endl;
  //cout<<T[1]<<endl;
  //cout<<bt[0]<<endl;
  //cout<<bt[1]<<endl;
  //cout<<F<<endl;
  //F=0.06;
  sig=sqrt(T[0])*v->volsto(T[0],lamb,F,X,alpha,K);
  cout<<"The implied Black-Sholes volatility is :"; 
  cout<<v->volsto(T[0],lamb,F,X,alpha,K)<<endl;
  d1=(log(F/X)+0.5*pow(sig,2))/sig;
  d2=(log(F/X)-0.5*pow(sig,2))/sig;
  n=1000;
  double e1,e2;
  n2b=int(n*(d2+50));
  n1b=int(n*(d1+50));
 
  double *l1;
  l1=new double[n1b+1];
  double *l2;
  l2=new double[n2b+1]; 
  e2=0.;
  for(i=0;i<=n2b;i++)
    {
      l2[i]=exp(-(-50.+i/double(n))*(-50.+i/double(n))/2);
    }

  for(i=3;i<=n2b-3;i++)
    {
      e2=e2+l2[i];
    }
  
  e2=(e2+(3./8.)*l2[0]+(7./6.)*l2[1]+(23./24.)*l2[2]+(23./24.)*l2[n2b-2]+(7./6.)*l2[n2b-1]+(3./8.)*l2[n2b])*1./double(n);
 e2=e2/sqrt(2*PI);

 e1=0;
 for(i=0;i<=n1b;i++)
   {
     l1[i]=exp(-(-50.+i/double(n))*(-50.+i/double(n))/2);
   }

 for(i=3;i<=n1b-3;i++)
   {
     e1=e1+l1[i];
   }
 
 e1=(e1+(3./8.)*l1[0]+(7./6.)*l1[1]+(23./24.)*l1[2]+(23./24.)*l1[n1b-2]+(7./6.)*l1[n1b-1]+(3./8.)*l1[n1b])*1./double(n);
 e1=e1/sqrt(2*PI);
 return(bt[1]*delta*(F*e1-X*e2));
 
}

double Caplet0::caplet0(double lamb,double alpha)
{
  int m,n1,i,n,n1b,n2b;
  double F;
  double delta;
  double *T;
  double sig;
  double d1;
  double d2;
  T=new double[2];
  double *bt;
  bt=new double[2];
  double *tt;
  double *btt;
  cout<<"enter the date T_{k} :";
  cin>>T[0];
  cout<<"enter the date T_{k+1} :";
  cin>>T[1];

  n1=10;
  tt=new double[n1];
  btt=new double[n1];
  for(i=0;i<n1;i++)
    { 
      tt[i]=0.25*(double)i;
    }
  btt[0]=1.000000;
  btt[1]=0.991982;
  btt[2]=0.982995;
  btt[3]=0.973132;
  btt[4]=0.962485;
  btt[5]=0.951141;
  btt[6]=0.939181;
  btt[7]=0.926683;
  btt[8]=0.913719;
  btt[9]=0.900356;

  bt=b->bond2(tt,btt,n1,T,2);

  
  delta=T[1]-T[0];
  F=(bt[0]/bt[1]-1.)/delta;
  //cout<<F<<endl;
  //F=0.06;
  sig=sqrt(T[0])*v->volconst(T[0],lamb,F,X,alpha);
  cout<<"The implied Black-Sholes volatility is :"; 
  cout<<v->volconst(T[0],lamb,F,X,alpha)<<endl;
  d1=(log(F/X)+0.5*pow(sig,2))/sig;
  d2=(log(F/X)-0.5*pow(sig,2))/sig;
  n=1000;
  double e1,e2;
  n2b=int(n*(d2+50));
  n1b=int(n*(d1+50));
 
  double *l1;
  l1=new double[n1b+1];
  double *l2;
  l2=new double[n2b+1]; 
  e2=0.;
  for(i=0;i<=n2b;i++)
    {
      l2[i]=exp(-(-50.+i/double(n))*(-50.+i/double(n))/2);
    }

  for(i=3;i<=n2b-3;i++)
    {
      e2=e2+l2[i];
    }
  
  e2=(e2+(3./8.)*l2[0]+(7./6.)*l2[1]+(23./24.)*l2[2]+(23./24.)*l2[n2b-2]+(7./6.)*l2[n2b-1]+(3./8.)*l2[n2b])*1./double(n);
 e2=e2/sqrt(2*PI);

 e1=0;
 for(i=0;i<=n1b;i++)
   {
     l1[i]=exp(-(-50.+i/double(n))*(-50.+i/double(n))/2);
   }
 
 for(i=3;i<=n1b-3;i++)
   {
     e1=e1+l1[i];
   }
 
 e1=(e1+(3./8.)*l1[0]+(7./6.)*l1[1]+(23./24.)*l1[2]+(23./24.)*l1[n1b-2]+(7./6.)*l1[n1b-1]+(3./8.)*l1[n1b])*1./double(n);
 e1=e1/sqrt(2*PI);
 return(bt[1]*delta*(F*e1-X*e2));
}

