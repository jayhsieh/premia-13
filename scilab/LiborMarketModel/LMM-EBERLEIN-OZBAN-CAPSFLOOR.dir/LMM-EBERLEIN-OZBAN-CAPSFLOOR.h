/*---------------------------------------------------------------*/
/*                 The Lévy Libor Model                          */
/*                                                               */
/*---------------------------------------------------------------*/
/*  Audrey Drif, Premia 2006                                     */
/*---------------------------------------------------------------*/


#include<iostream>
#include <cmath>
#include <complex>
#include <cstdlib>
using namespace std;
typedef complex<double> dcomplex;
const double PI = 3.141592653589793;

class bondprice{//calculate bondprice
public:
  bondprice(){}
  double* bond1(double*,int);//calculate bondprice with interest rate
  double* bond2(double*,double*,int,double*,int);//calculate bond price by interpolation
};

class CharExponent{//calculate X(u)
public:
  virtual dcomplex X(dcomplex u,int j,int n,double *lam)=0;
};

class Gamma : public CharExponent{//calculate Levy measure at time Tj for the case where the Levy density is a*exp(-lambda*x)/x, x>0
public:
  double lamb;//lambda
  double a;
  Gamma(double _lamb,double _a) : lamb(_lamb), a(_a){}
  dcomplex X(dcomplex u,int j,int n,double *lam);
};

class Alpha1 : public CharExponent{//calculate Levy measure at time Tj for the case where the Levy density is a*exp(-lambda*x)/x^2, x>0
public:
  double a;
  double lamb;//lambda
  Alpha1(double _a,double _lamb) : a(_a), lamb(_lamb){}
  dcomplex X(dcomplex u,int j,int n,double *lam);
};

class Alphag : public CharExponent{//calculate Levy measure at time Tj for the case where the Levy density is a*exp(-lambda*x)/x^(alpha+1), x>0
public:
  double a;
  double lamb;//lambda
  double alpha;
  Alphag(double _a,double _lamb,double _alpha) : a(_a), lamb(_lamb) ,alpha(_alpha){}
  dcomplex X(dcomplex u,int j,int n,double *lam);
};

class Templevyg : public CharExponent{//calculate Levy measure at time Tj for the case where the Levy density is (a-)*exp(-(lambda-)*|x|)/|x|^{1+(alpha-)}*1{x<0}+(a+)*exp(-(lambda+)*|x|)/|x|^{1+(alpha+)}*1{x>0}
public:
  double a1;//a+
  double lamb1;//lambda+
  double alpha1;//alpha+
  double a2;//a-
  double lamb2;//lambda-
  double alpha2;//alpha-
  Templevyg(double _a1,double _lamb1,double _alpha1,double _a2,double _lamb2,double _alpha2) : a1(_a1), lamb1(_lamb1) ,alpha1(_alpha1),a2(_a2),lamb2(_lamb2),alpha2(_alpha2){}
  dcomplex X(dcomplex u,int j,int n,double *lam);
};

class Templevy1 : public CharExponent{//calculate Levy measure at time Tj for the case where the Levy density is (a-)*exp(-(lambda-)*|x|)/|x|^{2}*1{x<0}+(a+)*exp(-(lambda+)*|x|)/|x|^{2}*1{x>0}
public:
  double a1;//a+
  double lamb1;//lambda+
  double a2;//a-
  double lamb2;//lambda-
  Templevy1(double _a1,double _lamb1,double _a2,double _lamb2) : a1(_a1), lamb1(_lamb1) ,a2(_a2),lamb2(_lamb2){}
  dcomplex X(dcomplex u,int j,int n,double *lam);
};

class Templevy : public CharExponent{//calculate Levy measure at time Tj for the case where the Levy density is (a-)*exp(-(lambda-)*|x|)/|x|*1{x<0}+(a+)*exp(-(lambda+)*|x|)/|x|*1{x>0}
public:
  double a1;//a+
  double lamb1;//lambda+
  double a2;//a-
  double lamb2;//lambda-
  Templevy(double _a1,double _lamb1,double _a2,double _lamb2) : a1(_a1), lamb1(_lamb1) ,a2(_a2),lamb2(_lamb2){}
  dcomplex X(dcomplex u,int j,int n,double *lam);
};

class laplace{//calculate bilateral transform of vk
public:
  double K;//strike
  laplace(double _K) : K(_K){}
  dcomplex bilaplace(double R,double u,double g);
}; 

class caplet{//calculate the time-0 price of the j-th caplet
public:
  CharExponent * Xg;
  
  laplace *la;
  int NN;
  double R; 
  caplet(CharExponent * _Xg,laplace * _la,double _R,int _NN) : Xg(_Xg),la(_la),R(_R),NN(_NN) {}
  dcomplex pxcaplet(int n,int j,double T,double g,double *lam,double *libo,double * t);
};

class cap{//calculate the time-0 price of cap
public:
  caplet * ca;
  bondprice *b;
  cap(caplet * _ca,bondprice * _b) : ca(_ca), b(_b){}
  dcomplex pxcap();
};

double* bondprice::bond1(double *t,int m)
{
  int i;
  double* bt;
  double rr;
  cout<<"enter the interest rate r:";
  cin>>rr;
  bt=new double[m]; 
  for(i=0;i<m;i++)
    {
      bt[i]=exp(-rr*t[i]);
    }
  return bt;
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
          if(tt[j]<=t[i] && tt[j+1]>t[i]) 
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

dcomplex Gamma::X(dcomplex u,int j,int n,double *lam)//calculate Levy measure at time Tj for the case where the Levy density is a*exp(-lambda*x)/x, x>0
{
  dcomplex I(0,1);
  int i;
  int k;
  dcomplex e;
  dcomplex aa;
  double lamm;
  
  lamm=lamb;
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm=lamm-lam[i];
      }
    
  
  aa=a*(u*lam[j-1]/(I*lamm)+log((I*lamm)/(u*lam[j-1]+I*lamm))-I*u*(-lam[j-1]/lamm+log(lamm/(lamm-lam[j-1]))));
  return aa;
}

dcomplex Alpha1::X(dcomplex u,int j,int n,double *lam)//calculate Levy measure at time Tj for the case where the Levy density is a*exp(-lambda*x)/x^2, x>0
{
  dcomplex I(0,1);
  int i;
  int k;
  dcomplex e;
  dcomplex aa;
  double lamm;
 
  lamm=lamb;
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm=lamm-lam[i];
      }
    
  
  aa=a*((lamm-I*u*lam[j-1])*log(1.-I*u*lam[j-1]/lamm)+I*u*lam[j-1]-I*u*((lamm-lam[j-1])*log(1.-lam[j-1]/lamm)+lam[j-1]));        
  return aa;
}

dcomplex Alphag::X(dcomplex u,int j,int n,double *lam)//calculate Levy measure at time Tj for the case where the Levy density is a*exp(-lambda*x)/x^(alpha+1), x>0
{
  dcomplex I(0,1);
  int i;
  int k;
  dcomplex e;
  dcomplex aa;
  double lamm;
 
  lamm=lamb;
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm=lamm-lam[i];
      }
    
  double h;
  double mm;
  double rr;
  h=-alpha;
  if  (alpha>0) {
    mm=int(alpha)+1.;
    
    for(i=1;i<mm;i++)
      {
        h=h*(h+i);
      }
    
    rr=tgamma(-alpha+mm)/h;
  }
 
  else{ rr=tgamma(-alpha);}
  aa=a*pow(lamm,alpha)*rr*((pow(1.-I*u*lam[j-1]/lamm,alpha)-1.+I*u*lam[j-1]*alpha/lamm)-I*u*(pow(1.-lam[j-1]/lamm,alpha)-1.+lam[j-1]*alpha/lamm));   
  return aa;
}

dcomplex Templevyg::X(dcomplex u,int j,int n,double *lam)//calculate Levy measure at time Tj for the case where the Levy density is a2*exp(-lambda2*|x|)/|x|^{1+alpha2}*1{x<0}+a1*exp(-lambda1*|x|)/|x|^{1+alpha1}*1{x>0}
{
  dcomplex I(0,1);
  int i;
  int k;
  dcomplex aa1;
  double lamm1;
  double h1;
  double mm1;
  double rr1;
  lamm1=lamb1;
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm1=lamm1-lam[i];
      }
  
  
  h1=-alpha1;
  if  (alpha1>0) {
    mm1=int(alpha1)+1.;
    
    for(i=1;i<mm1;i++)
      {
        h1=h1*(h1+i);
      }
    
    rr1=tgamma(-alpha1+mm1)/h1;
  }
  else{ rr1=tgamma(-alpha1);}
  aa1=a1*pow(lamm1,alpha1)*rr1*((pow(1.-I*u*lam[j-1]/lamm1,alpha1)-1.+I*u*lam[j-1]*alpha1/lamm1)-I*u*(pow(1.-lam[j-1]/lamm1,alpha1)-1.+lam[j-1]*alpha1/lamm1)); 
  
  dcomplex aa2;
  double lamm2;
  lamm2=lamb2;
  double h2;
  double mm2;
  double rr2;
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm2=lamm2+lam[i];
      }
  
  h2=-alpha2;
  if  (alpha2>0) {
    mm2=int(alpha2)+1.;
    
    for(i=1;i<mm2;i++)
      {
        h2=h2*(h2+i);
      }
    
    rr2=tgamma(-alpha2+mm2)/h2;
  }
  else{ rr2=tgamma(-alpha2);}
  aa2=a2*pow(lamm2,alpha2)*rr2*((pow(1.+I*u*lam[j-1]/lamm2,alpha2)-1.-I*u*lam[j-1]*alpha2/lamm2)-I*u*(pow(1.+lam[j-1]/lamm2,alpha2)-1.-lam[j-1]*alpha2/lamm2)); 
  return (aa1+aa2);
}

dcomplex Templevy::X(dcomplex u,int j,int n,double *lam)//calculate Levy measure at time Tj for the case where the Levy density is (a2)*exp(-(lambda2)*|x|)/|x|)*1{x<0}+(a1)*exp(-(lambda1)*x)/x)*1{x>0}
{
  dcomplex I(0,1);
  int i;
  int k;
  dcomplex aa1;
  double lamm1;
  lamm1=lamb1;
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm1=lamm1-lam[i];
      }
  aa1=a1*(u*lam[j-1]/(I*lamm1)+log((I*lamm1)/(u*lam[j-1]+I*lamm1))-I*u*(-lam[j-1]/lamm1+log(lamm1/(lamm1-lam[j-1]))));
  
  dcomplex aa2;
  double lamm2;
  lamm2=lamb2;
  
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm2=lamm2+lam[i];
      }
  
  aa2=a2*(-u*lam[j-1]/(I*lamm2)+log((I*lamm2)/(-u*lam[j-1]+I*lamm2))-I*u*(+lam[j-1]/lamm2+log(lamm2/(lamm2+lam[j-1]))));
    
  return (aa1+aa2);
}


dcomplex Templevy1::X(dcomplex u,int j,int n,double *lam)//calculate Levy measure at time Tj for the case where the Levy density is (a2)*exp(-(lambda2)*|x|)/|x|^{2})*1{x<0}+(a1)*exp(-(lambda1)*x)/x^{2})*1{x>0}
{
  dcomplex I(0,1);
  int i;
  int k;
  dcomplex aa1;
  double lamm1;
  lamm1=lamb1;
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm1=lamm1-lam[i];
      }
  aa1=a1*((lamm1-I*u*lam[j-1])*log(1.-I*u*lam[j-1]/lamm1)+I*u*lam[j-1]-I*u*((lamm1-lam[j-1])*log(1.-lam[j-1]/lamm1)+lam[j-1])); 
  
  dcomplex aa2;
  double lamm2;
  lamm2=lamb2;
  
  if (j<n) 
    for(i=j;i<=n-1;i++)
      {
        lamm2=lamm2+lam[i];
      }
  
  aa2=a2*((lamm2+I*u*lam[j-1])*log(1.+I*u*lam[j-1]/lamm2)-I*u*lam[j-1]-I*u*((lamm2+lam[j-1])*log(1.+lam[j-1]/lamm2)-lam[j-1])); 
  
  return (aa1+aa2);
  
}


dcomplex laplace::bilaplace(double R,double u,double g)//calculate L[vk](z) where z=R+I*u
{
  double KK;
  KK=g*K+1;
  dcomplex I(0,1);
  return(exp(log(KK)*(R+1.+I*u))/(-(R+1.+I*u))+KK*exp((R+I*u)*log(KK))/(R+I*u));
}


dcomplex caplet::pxcaplet(int n,int j,double T,double g,double *lam,double *libo,double *t)//calculate the time-0 price of the j-th caplet
{ dcomplex I(0,1);
  
 int n1;
 int i;
 dcomplex c,pri;
 dcomplex w,ww;
 double u[2*NN+1];
 dcomplex m1[2*NN+1];
 dcomplex m2[2*NN+1];
 dcomplex lp[2*NN+1];
 c=0.;
 for(i=0;i<=2*NN;i++)
   {
     u[i]=-100.+i*100./NN;
   }
 
 //calculate X(u[i]) at time Tj
 for(i=0;i<=2*NN;i++)
   {
     m2[i]=exp(t[j-1]*Xg->X(-u[i]+I*R,j,n,lam));
   }
  
 //calculate bilateral Laplace transform of vk
 for(i=0;i<=2*NN;i++)
   {
     m1[i]=la->bilaplace(R,u[i],g);
   }
  
 //calculate the time-0 price of the j-th caplet
 for(i=0;i<=2*NN;i++)
   {
     lp[i]=m1[i]*m2[i]*exp(I*u[i]*(-log(libo[j-1])));
   }

 for(i=3;i<=2*NN-3;i++)
   {
     c=c+lp[i];
   }
  
 c=(c+(3./8.)*lp[0]+lp[1]*(7./6.)+lp[2]*(23./24.)+lp[2*NN-1]*(7./6.)+lp[2*NN-2]*(23./24.)+(3./8.)*lp[2*NN])*double(100)/double(NN);
 pri=c*exp(-log(libo[j-1])*R)/(2*PI);
 return pri;
}

dcomplex cap::pxcap()//price cap
{
  int n,i;
  dcomplex ee;
  dcomplex eee;
  double T,g;
  double *t;
  double *tt;
  double *btt;
  double *bt;
  double *lam;
  int n1,m;
  double *libo;
  double KK;
  static char init[]="initialyield.dat";
  /*Name of the file where to read P(0, T) of the market.*/
  static FILE* Entrees; /*File variable of the code*/
  int etat;
  char ligne[20];
  char* pligne;
  double p, ttt;

  cout<<"enter the number of Cap dates :";
  cin>>n;
 
  bt=new double[n+1];
  t=new double[n+1];
  lam=new double[n];
  libo=new double [n];
  cout<<"enter the maturity of the Cap :";
  cin>>T;
  
  g=T/(double(n)+1);
 
  for(i=1;i<=n+1;i++)
    {
      t[i-1]=(i)*(T/(double(n)+1.));
    }
  
  for(i=0;i<n;i++)
    { 
      cout<<"enter the volatility of libor for the date:\n";
      cout<<t[i]<<endl;
      cout<<"Note that the sum of the volatility should not be greater or equal than lambda or min(lambda+,lambda-)!:";
      cin>>lam[i];
    }
 
  //calculate bond price
  cout<<"to calculate bond prices verify yields structure given in initialyealds.dat:\n";
  /* i est le nb de ligne lues */
  /* La ligne lue dans le fichier doit etre de la forme "0.943290 t=0.5" ou 0.943290 est un double pour le prix de B(0,t=0.5)*/ 
  
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
  bt=b->bond2(tt,btt,n1,t,n+1);


  //calculate 1+g*L(0,Tj)
  for(i=0;i<n;i++)
    {
      libo[i]=bt[i]/bt[i+1];
    }
  
  ee=0.;
  //calculate the price of cap
  for(i=1;i<=n;i++)
    {
      eee=bt[i]*ca->pxcaplet(n,i,T,g,lam,libo,t);
      ee=ee+eee;
    }
  
  
  delete [] t;
  delete [] bt;
  delete [] lam;
  delete [] libo;
  delete [] tt;
  delete [] btt;
  
  return(ee);
}

