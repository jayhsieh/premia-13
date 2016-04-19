/*---------------------------------------------------------------*/
/*     Exact pricing formulae for caps and swaptions in          */
/*               Lévy term strucure model                        */
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



class Volatility{//volatility stucture
 public:
  virtual double Sigma(double s, double t)=0;
};

class Ho_Lee : public Volatility{//Ho-Lee volatility structure, sigma(s,T)=tile{sigma}
 public:
  double sig;//tilde{sigma}
  Ho_Lee(double _sig) : sig(_sig){}
  double Sigma(double s,double t); 
};


class Vasicek : public Volatility{//Vasicek volatility structure, sigma(s,T)=tilde{sigma}*exp(-x(T-s))
 public:
  double x;
  double sig;//tilde{sigma}
  Vasicek(double _x,double _sig) : x(_x), sig(_sig){}
  double Sigma(double s,double t); 
};

class Moraleda : public Volatility{//Moraleda-Vorst volatility structure, sig(s,T)=tilde{sigma}*(1+gamma*T)*exp(-x(T-s))/(1+gamma*s)
 public:
  double sig;//tilde{sigma}
  double x;
  double gam;//gamma
  Moraleda(double _sig,double _x,double _gam) : sig(_sig), x(_x), gam(_gam){}
  double Sigma(double s,double t); 
};


class CharExponent{
 public:
  virtual dcomplex theta(dcomplex z)=0;
};

class Browmian : public CharExponent{//brownian motion
 public:
  double b;//drift
  double c;//coefficient brownian
  Browmian(double _b,double _c) : b(_b), c(_c){}
  dcomplex theta(dcomplex z);
};


class Gamma : public CharExponent{//Levy density a*exp(-lambda*x)/x,x>0
 public:
  double b,c;//parameters of the brownian motion
  double lamb;//lambda
  double a;
  Gamma(double _b,double _c,double _lamb,double _a) : b(_b), c(_c), lamb(_lamb), a(_a){}
  dcomplex theta( dcomplex z);
};

class Alpha1 : public CharExponent{//Levy density a*exp(-lambda*x)/x^2,x>0
 public:
  double b,c;//parameters of the brownian motion
  double lamb;//lambda
  double a;
  Alpha1(double _b,double _c,double _lamb,double _a) : b(_b), c(_c), lamb(_lamb), a(_a){}
  dcomplex theta(dcomplex z);
};

class Alphag : public CharExponent{//Levy density a*exp(-lambda*x)/x^{alpha+1},x>0}
 public:
  double b,c;//parameters of the brownian motion
  double lamb;//lambda
  double a;
  double alpha;
    Alphag(double _b,double _c,double _lamb,double _a,double _alpha):b(_b), c(_c), lamb(_lamb), a(_a), alpha(_alpha){}
  dcomplex theta( dcomplex z);
};

class Poisson : public CharExponent{//Brownian+Compound Poisson jumps, levy density lambda*exp(-(x-mu)^2)/2*var^2)/var*sqrt(2*pi)
 public:
  double b,c;//parameters of the brownian motion
  double lamb;//lambda
  double mu;
  double var;
  Poisson(double _b,double _c,double _lamb,double _mu,double _var) : b(_b), c(_c), lamb(_lamb), mu(_mu), var(_var){}
  dcomplex theta(dcomplex z);
};


class levy{//calculate levy measure
 public:
  virtual double intlevy()=0;
};

class lvgamma : public levy{//case gamma
 public:
  double lamb;
  double a;
  lvgamma(double _lamb,double _a) : lamb(_lamb),a(_a){}
  double intlevy();
};

class lalpha1 : public levy{//case alpha1
 public:
  double lamb;
  double a;
  lalpha1(double _lamb,double _a) : lamb(_lamb),a(_a){}
  double intlevy();
};

class lalphag : public levy{//case alphag
 public:
  double lamb;
  double a;
  double alpha;
  lalphag(double _lamb,double _a,double _alpha) : lamb(_lamb),a(_a),alpha(_alpha){}
  double intlevy(void);
};

class Integrale{//calculate the integrale of theta(z*sigma(s,T)+(1-z)*sigma(s,t))-theta(sigma(s,t)) from 0 to t
 public:
 virtual dcomplex integrale(double t, double T, dcomplex z)=0;
};

class IntegraleGeneric : public Integrale{ //numerical calculation
 public:
  Volatility * vol1;
  
  CharExponent * ce1;
  
  int n;
  IntegraleGeneric(int _n,Volatility *_vol1, CharExponent * _ce1) : n(_n), vol1(_vol1), ce1(_ce1){}
  dcomplex integrale(double t, double T, dcomplex z);
};

class IntegraleAnalytiqueBrownien : public Integrale{//analytical calculation for the case Brownian motion with Ho-Lee volatility
 public:
  Ho_Lee * vol;
  double b,c;//parameters of the browmian motion
  double sig;//tilde{sigma}
  IntegraleAnalytiqueBrownien (double _b,double _c,double _sig,Ho_Lee *_vol) : b(_b), c(_c),sig(_sig),vol(_vol){}
  dcomplex integrale(double t, double T, dcomplex z);
};

class IntegraleAnalytiqueGamma : public Integrale{//analytical calculation for the case Gamma with Ho-Lee volatility
 public:
  Ho_Lee * vol;
  double b,c;//parameters of the brownian motion
  double sig;//tilde{sigma}
  double a,lamb;//a,lambda
   IntegraleAnalytiqueGamma (double _b,double _c,double _a,double _lamb,double _sig, Ho_Lee *_vol) : b(_b), c(_c),a(_a),lamb(_lamb),sig(_sig),vol(_vol){}
  dcomplex integrale(double t, double T, dcomplex z);
};

class IntegraleAnalytiqueAlpha1 : public Integrale{//analytical calculation for the case Alpha1 with Ho_Lee volatility
 public:
  Ho_Lee * vol;
  double b,c;//parameters of the brownian motion
  double sig;//tilde{sigma}
  double a,lamb;//a,lambda
   IntegraleAnalytiqueAlpha1(double _b,double _c,double _a,double _lamb,double _sig, Ho_Lee *_vol) : b(_b), c(_c),a(_a),lamb(_lamb),sig(_sig),vol(_vol){}
  dcomplex integrale(double t, double T, dcomplex z);
};

class IntegraleAnalytiqueAlphag : public Integrale{//analytical calculation for the case Alphag with Ho_Lee volatility
 public:
  Ho_Lee * vol;
  double b,c;//parameters of the brownian motion
  double sig;//tilde{sigma}
  double a,lamb,alpha;//a,lambda,alpha
  IntegraleAnalytiqueAlphag(double _b,double _c,double _a,double _lamb,double _sig,double _alpha, Ho_Lee *_vol) : b(_b), c(_c),a(_a),lamb(_lamb),sig(_sig),alpha(_alpha),vol(_vol){}
  dcomplex integrale(double t, double T, dcomplex z);
};

class bondprice{//calculate bondprice
 public:
  bondprice(){}
  double* bond1(double*,int);//calculate bondprice with interest rate
  double* bond2(double*,double*,int,double*,int);//calculate bond price by interpolation
};

class cap{//the value of a call with strike K and maturity t on a bond which matures at T
 public:
  double K;
  int N;
  bondprice * b;
  double r;
  Integrale * int1;
  cap(int _N,double  _K,double _r,bondprice * _b,Integrale * _int1) : N(_N), K(_K), r(_r),b(_b),int1(_int1){}
  dcomplex pxcap();
};

class swaption{//the time-0 price of a call with strike price 1 and maturity t on a bond with maturity Tn paying to its owner an amount of C1,...,Cn at the dates T1,...,Tn
 public:
  int nn;
  bondprice* bb;
  double r;
  Integrale * int1;
  swaption(int _nn,double _r,bondprice * _bb ,Integrale * _int1) : nn(_nn),r(_r),bb(_bb), int1(_int1){}
  dcomplex pxswaption(); 
};

double sum(double*,double*,double,int);


double lvgamma::intlevy()
{
  return(a*exp(-lamb)/lamb);
}

double lalpha1::intlevy()
{//simpson's rule
  double e;
  e=0.;
  int i;
  double l[1001];
  for(i=0;i<=1000;i++)
    {
      l[i]=exp(-lamb*(1+i/100.))/(1+i/100.);
    }
  for(i=3;i<=997;i++)
    {
      e=e+l[i];
    }
   
e=a*(e+(3./8.)*l[0]+(7./6.)*l[1]+(23./24.)*l[2]+(3./8.)*l[1000]+(7./6.)*l[999]+(23./24.)*l[998])/100.;
 return(e);
}

double lalphag::intlevy()
{
  double e;
  double alphabis;
  alphabis=-alpha;
  e=0.;
  int i;
  double l[1001];
  for(i=0;i<=1000;i++)
    {
      if(alpha<0) l[i]=exp(-lamb*(1.+i/100.))*pow(1.+i/100.,alphabis);
      else  l[i]=exp(-lamb*(1+i/100.))/pow(1+i/100.,alpha);
    }
  for(i=3;i<=997;i++)
    {
      e=e+l[i];
    }
   
  e=a*(e+(3./8.)*l[0]+(7./6.)*l[1]+(23./24.)*l[2]+(3./8.)*l[1000]+(7./6.)*l[999]+(23./24.)*l[998])/100.;
 return(e);
}



double Ho_Lee::Sigma(double s,double t)
{
  return(sig*(t-s));
}

double Vasicek::Sigma(double s,double t)
{
  return(sig*(1.-exp(-x*(t-s)))/x);
}

double Moraleda::Sigma(double s,double t)
{
  return(sig*exp(x*s)*(exp(-x*s)*(1.+gam*s+gam/x)+exp(-x*t)*(-1.-gam*t-gam/x))/(x*(1+gam*s)));
}

dcomplex Browmian::theta(dcomplex z)
{
  return(b*z+0.5*c*pow(z,2));
}

dcomplex Gamma::theta(dcomplex z)
{
  return(b*z+0.5*c*pow(z,2)+a*(-z/lamb+log(lamb/(lamb-z))));
}

dcomplex Alpha1:: theta(dcomplex z)
{
  return(b*z+0.5*c*pow(z,2)+a*((lamb-z)*log(1.-z/lamb)+z));
}

dcomplex Alphag ::theta(dcomplex z)
{
  double h;
  int mm;
  int i;
  double rr;
  h=alpha;
  if  (alpha>0) {
    mm=int(alpha)+1;
    
    for(i=1;i<mm;i++)
      {
	h=h*(h+i);
      }
    
    rr=tgamma(-alpha+mm)/h;
  }
  else{ rr=tgamma(-alpha);}
  return(b*z+0.5*c*pow(z,2)+pow(lamb,alpha)*rr*a*(pow(1.-z/lamb,alpha)-1.+alpha*z/lamb));
}

dcomplex Poisson::theta(dcomplex z)
{
  return(b*z+0.5*c*pow(z,2)+lamb*(exp(mu*z+pow(var*z,2)/2.)-1.));
}

dcomplex IntegraleGeneric::integrale(double t,double T,dcomplex z)
{
  dcomplex l[n+1];
  dcomplex u;
  int i;
  double r1[n+1];
  double r2[n+1];
  dcomplex r3[n+1];
  dcomplex r4[n+1];
  dcomplex e;
  e=0.;
  
  for(i=0;i<=n;i++)//Simpson's rule
    {
      r1[i]=vol1->Sigma(i*t/n,T);
      r2[i]=vol1->Sigma(i*t/n,t);
      r3[i]=ce1->theta(z*r1[i]+(r2[i])*(1.-z));
      r4[i]=ce1->theta(r2[i]);
      l[i]=r3[i]-r4[i];
    }
 
  for(i=3;i<=n-3;i++)
    {
      e=e+l[i];
    }
  
  e=(e+(3./8.)*l[0]+(7./6.)*l[1]+(23./24.)*l[2]+(3./8.)*l[n]+(7./6.)*l[n-1]+(23./24.)*l[n-2])*t/double(n);
  
  return e;
} 

dcomplex IntegraleAnalytiqueBrownien::integrale(double t,double T,dcomplex z)
{
  return(z*vol->Sigma(t,T)*b*t+0.5*z*z*sig*vol->Sigma(t,T)*(T-t)*c*t+0.5*z*sig*vol->Sigma(t,T)*c*t*t);
}

dcomplex IntegraleAnalytiqueGamma::integrale(double t,double T,dcomplex z)
{
  return(z*vol->Sigma(t,T)*b*t+0.5*z*z*sig*vol->Sigma(t,T)*(T-t)*c*t+0.5*z*sig*vol->Sigma(t,T)*c*t*t-a*z*vol->Sigma(t,T)*t/lamb+(a/sig)*(lamb*log(lamb)-(lamb-sig*t)*log(lamb-sig*t)-sig*t-(lamb-z*vol->Sigma(t,T))*log(lamb-z*vol->Sigma(t,T))+(lamb-z*vol->Sigma(t,T)-sig*t)*log(lamb-z*vol->Sigma(t,T)-sig*t)+sig*t));
}

dcomplex IntegraleAnalytiqueAlpha1::integrale(double t,double T,dcomplex z)
{
  return(z*vol->Sigma(t,T)*b*t+0.5*z*z*sig*vol->Sigma(t,T)*(T-t)*c*t+0.5*z*sig*vol->Sigma(t,T)*c*t*t+a*z*vol->Sigma(t,T)*t+a*(pow(lamb-z*vol->Sigma(t,T),2)*log(1.-z*vol->Sigma(t,T)/lamb)/(2*sig)-pow(lamb-(z*vol->Sigma(t,T)+sig*t),2)*log(1.-(z*vol->Sigma(t,T)+sig*t)/lamb)/(2*sig)+pow(lamb-sig*t,2)*log(1-sig*t/lamb)/(2*sig)+0.5*z*vol->Sigma(t,T)*t));
}

dcomplex IntegraleAnalytiqueAlphag::integrale(double t,double T,dcomplex z)
{
  double rr,h;
  int mm,i;
  h=alpha;
  if  (alpha>0) {
    mm=int(alpha)+1;
    for(i=1;i<mm;i++)
      {
	h=h*(h+i);
      }  
    //calculate zeta
    rr=tgamma(-alpha+mm)/h;
  }
  else{ rr=tgamma(-alpha);} 
  
  return(z*vol->Sigma(t,T)*b*t+0.5*z*z*sig*vol->Sigma(t,T)*(T-t)*c*t+0.5*z*sig*vol->Sigma(t,T)*c*t*t+rr*a*pow(lamb,alpha)*(z*vol->Sigma(t,T)*t*alpha/lamb+lamb/((alpha+1)*sig)*(-1.+pow(1.-z*vol->Sigma(t,T)/lamb,alpha+1)-pow(1.-(z*vol->Sigma(t,T)+sig*t)/lamb,alpha+1)+pow(1.-sig*t/lamb,alpha+1))));
}

double* bondprice::bond1(double *t,int m)
{//B(0,t)=exp(-rrt)
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
  return (bt);
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


dcomplex cap::pxcap()
{
  dcomplex I(0,1);
  int mm,i,j;
  int n1,n2;
  dcomplex q;
  dcomplex p;
  dcomplex s;
  double *tt;
  double *btt;
  double *t;
  double *bt;
  dcomplex m1[2*N+1];
  dcomplex m2[2*N+1];
  dcomplex  z[2*N+1];
  double u[2*N+1];
  double v[2*N+1];
  int m;
   static char init[]="initialyield.dat";
  /*Name of the file where to read P(0, T) of the market.*/
  static FILE* Entrees; /*File variable of the code*/
  int etat;
  char ligne[20];
  char* pligne;

  t=new double[2];
  bt=new double[2];
  
  //calculate bond price
  cout<<" enter the option maturity t :";
  cin>>t[0];
  cout<<"enter the bond maturity T :";
  cin>>t[1]; 
  
 
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
    bt=b->bond2(tt,btt,n1,t,2);

  for(i=0;i<=2*N;i++)
    {
      u[i]=-150.+i*150./N;
      z[i]=-r-u[i]*I;
    }
  
  for(i=0;i<=2*N;i++)
    {
      m2[i]=1./((r+u[i]*I)*(1.+r+u[i]*I));
    }
  
  q=log(bt[0]/bt[1])+log(K)+(int1->integrale(t[0],t[1],1.));
  
  //calculate the random generating function
  for(i=0;i<=2*N;i++)
    {
      m1[i]=exp(int1->integrale(t[0],t[1],z[i]));

    }
  
  //Simpson's rule
  dcomplex l[2*N+1];
  dcomplex e;
  e=0.;
  
  for(i=0;i<=2*N;i++)
    {
      l[i]=exp((r+u[i]*I)*q)*m1[i]*m2[i];
    }
  
  for(i=3;i<=2*N-3;i++)
    {
      e=e+l[i];
    }

  e=(e+(3./8.)*l[0]+l[1]*(7./6.)+l[2]*(23./24.)+l[2*N-1]*(7./6.)+l[2*N-2]*(23./24.)+(3./8.)*l[2*N])*double(150)/double(N);
  s=e*K*bt[0]/(2.*PI);
  return(s); 
}

dcomplex swaption::pxswaption()
{
  dcomplex prix;
  dcomplex p;
  dcomplex I(0,1);
  int i,j;
  double xx,y,zz;
  int nb;
  dcomplex m1[2*nn+1];
  dcomplex z[2*nn+1];
  double u[2*nn+1];
  double *tt;
  double *btt;
  double *t;
  double *bt;
  double *b;
  double *s;
  double *d;
  dcomplex L[2*nn+1];
  int m,n1;
  cout<<"enter number of coupons :";
  cin>>nb;
  t=new double[nb+1];
  bt=new double[nb+1];
  d=new double[nb]; 
  s=new double[nb]; 
  b=new double[nb];
   static char init[]="initialyield.dat";
  /*Name of the file where to read P(0, T) of the market.*/
  static FILE* Entrees; /*File variable of the code*/
  int etat;
  char ligne[20];
  char* pligne;

  //calculate bond price 
  cout<<"enter the option date t:";  
  cin>>t[0];
  cout<<"enter the bond maturity t:";  
  cin>>t[nb];

  for(i=1;i<nb+1;i++)
   t[i]=t[0]+(double)i*(t[nb]-t[1])/(double)nb;

    for(i=0;i<nb;i++)
    {
      cout<<"enter coupon C[i] at time :";
      cout<<t[i+1]<<endl;
      cin>>s[i]; 
    }
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
    bt=bb->bond2(tt,btt,n1,t,nb+1);

  for(i=0;i<=2*nn;i++)
    {
      u[i]=-150.+i*150./nn;
      z[i]=-r-u[i]*I;
    }
  
  //calculate Di
  for(i=0;i<nb;i++)
    {
      d[i]=(bt[i+1]*s[i])*real(exp(-int1->integrale(t[0],t[i+1],1.)))/(bt[0]);
    }
  
  //calculate the random generating function
  for(i=0;i<=2*nn;i++)
    {
      m1[i]=exp(int1->integrale(t[0],t[nb],z[i]));
    }
  
  //calculate bi
  for(i=0;i<nb;i++)
    {
      b[i]=(t[i+1]-t[0])/(t[nb]-t[0]);
    }
  
  //calculate zero
  xx=0;
  y=-1;
  if (abs(sum(b,d,xx,nb)-1.)< 0.2)  xx=xx;
  if (sum(b,d,xx,nb)-1.>0.) 
    {
      while(sum(b,d,y,nb)-1.>0.)
	{
	  
	  xx=y;
	  y=y-1.;
	}
      do
	{
	  zz=(xx+y)/2.;
	  if  (sum(b,d,zz,nb)-1.>0.)  xx=zz;
	  else y=zz;
	}
      while(abs(sum(b,d,zz,nb)-1.)>0.2);
    }
  if(sum(b,d,xx,nb)-1.<0.)  
    { 
    while (sum(b,d,y,nb)-1.<0.)  
      {
	xx=y;
	y=y+1;
      }
    do
      {
	zz=(xx+y)/2;
       if  (sum(b,d,zz,nb)-1.>0.) y=zz;
       else xx=zz;
      }
    while(abs(sum(b,d,zz,nb)-1.)>0.2);
    }
 
  
  //calculate bilateral laplace transform
  for(i=0;i<=2*nn;i++)
    {   
      for(j=0;j<nb;j++)
	{
	  L[i]=L[i]-d[j]*exp(b[j]*zz)*(1./(b[j]+r+u[i]*I));
	}
      L[i]=exp(r*zz)*(L[i]+1./(-z[i]));
    }
  //Simpson's rule
  dcomplex ll[2*nn+1];
  dcomplex e;
  e=0.;
  
  for(i=0;i<=2*nn;i++)
    {
      ll[i]=exp(u[i]*zz*I)*m1[i]*L[i];
    }
  
  for(i=3;i<=2*nn-3;i++)
    {
      e=e+ll[i];
    }
  
  e=(e+(3./8.)*ll[0]+(7./6.)*ll[1]+(23./24.)*ll[2]+(3./8.)*ll[2*nn]+(7./6.)*ll[2*nn-1]+(23./24.)*ll[2*nn-2])*double(150)/double(nn);  
  
  prix=e*double(bt[0]/(2*PI));
  return(prix);
}

double sum(double* b,double* d,double x,int n)
{
  int i;
  double y;
  y=0.;
  for(i=0;i<n;i++)
    {
      y=y+d[i]*exp(b[i]*x);
    }
  return y;
}

