/*---------------------------------------------------------------*/
/*                 The Lévy Libor Model                          */
/*                                                               */
/*---------------------------------------------------------------*/
/*  Audrey Drif, Premia 2006                                     */
/*---------------------------------------------------------------*/


#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include <fstream>
#include "LMM-EBERLEIN-OZBAN-CAPSFLOOR.h"


using namespace std;
typedef complex<double> dcomplex;

int main()
{
  typedef complex<double> dcomplex;
  dcomplex I(0,1);
  dcomplex w;
  double K;
  double mu,a,var;
  double a1,a2,lamb,lamb1,lamb2,alpha1,alpha2,beta,alpha;
  int m;
  bondprice *b;
  Gamma *ml;
  Alpha1* p;
  Alphag *Ag;
  Templevyg* Tpg;
  Templevy1* Tp1;
  Templevy* Tp;
  laplace *lap;
  caplet *cpt;
  cap *cp;
   cout<<"\n";
  cout<<"Pricing Cap using LMM Levy models\n";
  cout<<"\n";
  cout<<"Choose Levy measure:\n";
  cout<<"\n";
  cout<<"(0)  case where the density Lévy is a*exp(-lambda*x)/x^2:\n";
  cout<<"(1)  case where the density Lévy is a*exp(-lambda*x)/x 1:\n";
  cout<<"(2)  case where the density Lévy is a*exp(-lambda*x)/x^{1+alpha}:\n";
  cout<<"(3)  case where the density Lévy is (a-)*exp(-(lambda-)*|x|)/|x|^{1+alpha-}*1{x<0}+(a+)*exp(-(lambda+)*x)/x^{1+alpha+}*1{x>0}: \n";
  cout<<"(4)  case where the density Lévy is (a-)*exp(-(lambda-)*|x|)/|x|)*1{x<0}+(a+)*exp(-(lambda+)*x)/x)*1{x>0}: \n";
  cout<<"(5) case where the density Lévy is (a-)*exp(-(lambda-)*|x|)/|x|^{2}*1{x<0}+(a+)*exp(-(lambda+)*x)/x^{2}*1{x>0}: \n";
  cout<<"\n";	       
  cin>>m;
  
  b=new bondprice();
  
  switch(m)
    {

    case 0://a*exp(-lambda*x)/x^2,x>0
      cout<<"enter the parameter a>0 of this Lévy density :";
      cin>>a;
      cout<<"enter the parameter lambda>0 of this Lévy density :";
      cin>>lamb;
      cout<<"enter the strike K :";
      cin>>K;
      lap=new laplace(K); 
      p=new Alpha1(a,lamb);
      cpt=new caplet(p,lap,-1.5,3000);
      cp=new cap(cpt,b);
      w=cp->pxcap();
      cout<<real(w)<<endl;
      break;
    case 1://a*exp(-lambda*x)/x,x>0
      cout<<"enter the parameter a>0 of this Lévy density :";
      cin>>a;
      cout<<"enter the parameter lambda>0 of this Lévy density :";
      cin>>lamb;
      cout<<"enter the strike K :";
      cin>>K;
      lap=new laplace(K); 
      ml=new Gamma(lamb,a);
      cpt=new caplet(ml,lap,-1.5,3000);
      cp=new cap(cpt,b);
      w=cp->pxcap();
      cout<<real(w)<<endl;
      break;
    case 2://a*exp(-lambda*x)/x^{alpha+1},x>0
      cout<<"enter the parameter a>0 of this Lévy density :";
      cin>>a;
      cout<<"enter the parameter lambda>0 of this Lévy density :";
      cin>>lamb;
      cout<<"enter the parameter alpha(!=0 and !=1) of this Lévy density :";
      cin>>alpha;
      cout<<"enter the strike K :";
      cin>>K;
      lap=new laplace(K); 
      Ag=new Alphag(a,lamb,alpha);
      cpt=new caplet(Ag,lap,-1.5,3000);
      cp=new cap(cpt,b);
      w=cp->pxcap();
      cout<<real(w)<<endl;
      break;
    case 3://(a-)*exp(-(lambda-)*|x|)/|x|^(1+(alpha-))*1{x<0}+(a+)*exp(-(lambda+)*x)/x^(1+(alpha+))*1{x>0}
      cout<<"enter the parameter a+>0 of this Lévy density :";
      cin>>a1;
      cout<<"enter the parameter a->0 of this Lévy density:";
      cin>>a2;
      cout<<"enter the parameter lambda+>0 of this Lévy density :";
      cin>>lamb1;
      cout<<"enter the parameter lambda->0 of this Lévy density :";
      cin>>lamb2;
      cout<<"enter the parameter alpha+>0 of this Lévy density :";
      cin>>alpha1;
      cout<<"enter the parameter alpha->0 of this Lévy density : ";
      cin>>alpha2;
      cout<<"enter the strike K :";
      cin>>K;
      lap=new laplace(K); 
      Tpg=new Templevyg(a1,lamb1,alpha1,a2,lamb2,alpha2);
      cpt=new caplet(Tpg,lap,-1.5,3000);
      cp=new cap(cpt,b);
      w=cp->pxcap();
      cout<<real(w)<<endl;
      break;
    case 4://(a-)*exp(-(lambda-)*|x|)/|x|)*1{x<0}+(a+)*exp(-(lambda+)*x)/x)*1{x>0}
      cout<<"enter the parameter a+>0 of this Lévy density :";
      cin>>a1;
      cout<<"enter the parameter a->0 of this Lévy density:";
      cin>>a2;
      cout<<"enter the parameter lambda+>0 of this Lévy density :";
      cin>>lamb1;
      cout<<"enter the parameter lambda->0 of this Lévy density :";
      cin>>lamb2;
     cout<<"enter the strike K :";
      cin>>K;
      lap=new laplace(K); 
      Tp=new Templevy(a1,lamb1,a2,lamb2);
      cpt=new caplet(Tp,lap,-1.5,3000);
      cp=new cap(cpt,b);
      w=cp->pxcap();
      cout<<real(w)<<endl;
      break;
 case 5://(a-)*exp(-(lambda-)*|x|)/|x|^{2})*1{x<0}+(a+)*exp(-(lambda+)*x)/x^{2})*1{x>0}
      cout<<"enter the parameter a+>0 of this Lévy density :";
      cin>>a1;
      cout<<"enter the parameter a->0 of this Lévy density:";
      cin>>a2;
      cout<<"enter the parameter lambda+>0 of this Lévy density :";
      cin>>lamb1;
      cout<<"enter the parameter lambda->0 of this Lévy density :";
      cin>>lamb2;
     cout<<"enter the strike K :";
      cin>>K;
      lap=new laplace(K); 
      Tp1=new Templevy1(a1,lamb1,a2,lamb2);
      cpt=new caplet(Tp1,lap,-1.5,3000);
      cp=new cap(cpt,b);
      w=cp->pxcap();
      cout<<real(w)<<endl;
 
    }
}




