#include "costFunction.h"
#include "pnl/pnl_complex.h"
#include <stdlib.h>
#include "utils/intg.h"
#include "bsvanillas.h"
#include "utils1.h"

// Variables global 
static int initOk=0,TypeNorme,TypeModel;
static int LogNorme = 0;
static int MelangeNormes = 0;
static int sigmaimpok=0,sigmaimpko=0;
static DataMarket *DM;


//=============================== SYNTHETIC DATAMARKET ====================

int CreateSyntheticDM(DataMarket *DMS, int nt, double *_T, double *ri, double *di, double St0, int nk, double *K, int *TypeOpt, int dimx, double *x)
{
	int i, j, nc, j1;
	int nbdata, ibidon;
	double price, delta;
  //
  nbdata = nk*nt;

  if(DMS == NULL) 
    DMS = CreateDataMarket(nbdata);
  //
  DMS->nbdata     = nbdata;
  for(i=0;i<nt;i++)
      for(j=0; j<nk; j++)
	{
	  nc = i * nk + j;
	  DMS->TypeOpt[nc] = TypeOpt[j];
	  DMS->St0[nc]     = St0;
	  DMS->K[nc]       = K[j];
	  DMS->T[nc]       = _T[i];
	  DMS->r[nc]       = ri[i];
	  DMS->d[nc]       = di[i];
	  DMS->wi[nc]      = 1.0;
	  if(TypeOpt[j] == 1) {
		ibidon=FT_Call_Heston(St0, K[j], _T[i], ri[i], di[i], x[0], x[1], x[2], x[3], x[4], &price, &delta);
	      } 
	  else {
		ibidon=FT_Put_Heston(St0, K[j], _T[i], ri[i], di[i], x[0], x[1], x[2], x[3], x[4], &price, &delta);
	      }
	  DMS->PrixObs[nc] = price;
	  DMS->SigmaImpObs[nc] = SigmaImplicite( TypeOpt[j], St0, _T[i], K[j], ri[i], di[i], price, &j1);
	}
  //
  
  return 1;
}

//===================================== COST FUNCTION ====================
// 
double costFunction(int dimx, double *x)
{
  //
  int i;
  double fx=0,fxi;
  int ibidon;
  double V0,kappa,theta,sigmav,rho,prix,delta;
  double St0,K,T,prixobs,sigma_imp_obs;
  int nbok,nbko;
  
  double ri,di;
  nbok = 0;
  nbko = 0;

	  if(initOk != 1) {
		printf("Error in costFunction, paramaetres are not initialized.\n ");	
		exit(-1);
	  }
	 
	  V0     = x[0];
	  kappa  = x[1];
	  theta  = x[2];
	  sigmav = x[3];
	  rho    = x[4];

	  for(i=0;i<DM->nbdata;i++)
	    {
	      //
	      St0 = DM->St0[i];
	      T   = DM->T[i];
	      K   = DM->K[i];
	      ri = DM->r[i];
	      di = DM->d[i];
	      //
	      if(DM->TypeOpt[i] == 1) {
		ibidon=FT_Call_Heston(St0,K,T,ri,di,V0,kappa,theta,sigmav,rho,&prix,&delta);
	      } else {
		ibidon=FT_Put_Heston(St0,K,T,ri,di,V0,kappa,theta,sigmav,rho,&prix,&delta);
	      }
	      //
		if(prix<0.0) {prix = 1e-14;}
	      prixobs = DM->PrixObs[i];
	      sigma_imp_obs = DM->SigmaImpObs[i];
	      fxi = NormeMarket(TypeNorme, DM->TypeOpt[i], St0, T, K, ri, di, prix, prixobs, sigma_imp_obs, &DM->SigmaImp[i]);
	      //
	      DM->PrixMod[i] = prix;

	      fx = fx + fxi*DM->wi[i];
	    }
  return fx;
}

//
void initParams(int nbdata,DataMarket *_DM,int _TypeNorme,int _TypeModel,int _LogNorme,int _MelangeNormes)
{
  int i;
  //
  if(DM == NULL) 
    DM = CreateDataMarket(nbdata);
  //
  DM->nbdata     = _DM->nbdata;

  for(i=0;i<nbdata;i++)
	{
	  DM->TypeOpt[i] = _DM->TypeOpt[i];
	  DM->St0[i]     = _DM->St0[i];
	  DM->K[i]       = _DM->K[i];
	  DM->T[i]       = _DM->T[i];
	  DM->r[i]       = _DM->r[i];
	  DM->d[i]       = _DM->d[i];
	  DM->wi[i]      = _DM->wi[i];
	  DM->PrixObs[i] = _DM->PrixObs[i];
	  DM->SigmaImpObs[i] = _DM->SigmaImpObs[i];
	}
  //
  TypeNorme     = _TypeNorme;
  TypeModel     = _TypeModel;
  LogNorme      = _LogNorme;
  MelangeNormes = _MelangeNormes;
  //
  initOk        = 1;
  // 
}

//=================================== NORME MARKET ====================

double NormeMarket(int TypeNorme,int TypeOpt,double St0,double T,double K,double r,double divid,double prix,double prixobs,double sigma_imp_obs,double *sigmaimp)
{
  
  int j1;
  double fx,sigma_imp,sigma_imp0;
  
	  sigma_imp     = SigmaImplicite(TypeOpt,St0,T,K,r,divid,prix,&j1);
	  *sigmaimp = sigma_imp;

  // Case 1 : (absolute error)^2
  if(TypeNorme == 1)
	{
	  fx     = (prix - prixobs);
	  fx     = fx * fx;
	}
  // Cas 2 : (relative error)^2
  else if(TypeNorme == 2)
	{
	  fx     = (prix - prixobs)/prixobs;
	  fx     = fx * fx;
	}
  // Cas 3 : (absolute error in SigmaImplicites)^2
  else if(TypeNorme == 3)
	{
	  if(j1==-1) {
		sigmaimpko++;
		fx = NormeMarket(2,TypeOpt,St0,T,K,r,divid,prix,prixobs,sigma_imp_obs,&sigma_imp);
		} 
	else {
		sigmaimpok++;
		sigma_imp0 = 1.e-2;
		
		fx     = (sigma_imp - sigma_imp_obs)/sigma_imp0;
		fx     = fx * fx;
	  }
	  
	} 
  // Case 4: logarithmic norm
  else if(TypeNorme == 4)
	{
		fx = log(prix/prixobs)*(prix - prixobs);
	}
  else {
	  printf("Error in NormeMarket, TypeNorme = %d is unknown \n",TypeNorme);  
	  exit(-1);
	}
  //  
  return fx;
}
//====================================FINALS============================
void checksigmaimp()
{
  printf("SigmImp OK = %d  Sigma Imp KO = %d \n",sigmaimpok,sigmaimpko);
}
//
void FreeCostFunction()
{
  FreeDataMarket(DM);
  DM = NULL;  
}
//=============================== PRINT RESULTS ========================
void AffichePrixDataMarket()
{
  PrintDataMarket(DM);
}

//
void AfficheDataMarketFile()
{
  PrintDataMarketFile(DM);
}

// ===============================VANILLAS HESTON======================

static double T,sigma,rho,k, v, r, divid, teta, lambda, S, K;

static double charact_funct1(double uu)
{
  double a,b,rs,rsp,sig,tau,tpf1,tpf2, f10, c0, d0;
  dcomplex g,z,w,tp1,tp2,DD,CN, ans,d,expo;


  tau=T;
  a=k*teta;
  rs=rho*sigma;
  rsp=rs*uu;
  sig=sigma*sigma;

  b=k+lambda-rs;
  if(uu==0)
    {
      if(b==0)
 {
   c0=a*T*T/4.0;
          d0=T/2.0;
 }
      else
        {
          c0=0.5*a*(exp(-b*T)+b*T - 1.0)/b/b;
          d0=0.5*(1.0-exp(-b*T))/b;
 }
      f10=log(S/K)+(r-divid)*T+c0+d0*v;

      return f10;
    }
  z=Complex(-b,rsp);
  z=Cmul(z,z);
  w=RCmul(sig,Complex(-uu*uu,uu));
  d=Csqrt(Csub(z,w));
  tp1=Complex(d.r+b,d.i-rsp);
  tp2=Complex(-d.r+b,-d.i-rsp);
  g=Cdiv(tp2,tp1);

 expo=Cexp(RCmul(-tau,d));
     DD=Csub(Complex(1,0), expo);
     DD=Cdiv(DD,Csub(Complex(1,0),Cmul(g,expo)));
     DD=Cmul(DD,RCmul(1.0/sig,tp2));

  CN=Csub(Cmul(g,expo), Complex(1,0));
  CN=Cdiv(CN,Csub(g, Complex(1,0) ));
  tpf1=a*(tau*tp2.r-2.0*Clog(CN).r)/sig;
  tpf2=a*(tau*tp2.i-2.0*Clog(CN).i)/sig;

  tpf2+=(r-divid)*uu*tau;
  ans=Complex(tpf1+v*DD.r,tpf2+v*DD.i+uu*log(S));
  ans=Cmul(Cexp(ans),Cexp(Complex(0,-uu*log(K))));
  ans=Cdiv(ans,Complex(0,uu));
  
  return ans.r;
}

static double charact_funct2(double uu)
{
  double a,b,rsp,sig,tau,tpf1,tpf2, f20, c0, d0;
  dcomplex g,z,w,tp1,tp2,DD,CN, ans,d,expo;
  
  tau=T;
  a=k*teta;
  rsp=rho*sigma*uu;
  sig=sigma*sigma;

  b=k+lambda;
  if(uu==0)
    {
      c0=0.5*a*(exp(-b*T)+b*T - 1.0)/b/b;
      d0=0.5*(1.0-exp(-b*T))/b;
      f20=log(S/K)+(r-divid)*T-c0-d0*v;

      return f20;
    }
  z=Complex(-b,rsp);
  z=Cmul(z,z);
  w=RCmul(sig,Complex(-uu*uu,-uu));
  d=Csqrt(Csub(z,w));
  tp1=Complex(d.r+b,d.i-rsp);
  tp2=Complex(-d.r+b,-d.i-rsp);
  g=Cdiv(tp2,tp1);
  expo=Cexp(RCmul(-tau,d));
  DD=Csub(Complex(1,0),expo);
  DD=Cdiv(DD,Csub(Complex(1,0),Cmul(g,expo)));
  DD=Cmul(DD,RCmul(1.0/sig,tp2));

  CN=Csub(Cmul(g,expo), Complex(1,0));
  CN=Cdiv(CN,Csub(g, Complex(1,0) ));
  tpf1=a*(tau*tp2.r-2.0*Clog(CN).r)/sig;
  tpf2=a*(tau*tp2.i-2.0*Clog(CN).i)/sig;
  
  tpf2+=(r-divid)*uu*tau;
  ans=Complex(tpf1+v*DD.r,tpf2+v*DD.i+uu*log(S));
  ans=Cmul(Cexp(ans),Cexp(Complex(0,-uu*log(K))));
  ans=Cdiv(ans,Complex(0,uu));

  return ans.r;
}

static double  probabilities(int n)
{
  double tp, cinf, f0, Lamb, abserr, temp;
  int i;
  
  cinf=sqrt(1.0-rho*rho)/sigma*(v+k*teta*T);
  
  if(n==1)
    {
      f0=charact_funct1(0.0);

      Lamb=2.0*(log(fabs(f0))+12.0*log(10.0))/cinf;

	temp=0.0;
	Lamb=Lamb/100.0;
	intg(0.0, Lamb, charact_funct1, 1e-14, 1e-10, &tp, &abserr);
		temp += tp;
	for(i=1; i<101;i++)
	{
		intg(i*Lamb, (i+1)*Lamb, charact_funct1, 1e-14, 1e-10, &tp, &abserr);
		temp += tp;

	}
	
      intg(0.0, Lamb, charact_funct1, 1e-14, 1e-10, &tp, &abserr);

      tp=0.5+temp/M_PI;
      return tp;

    }
  else
    {
      f0=charact_funct2(0.0);

      Lamb=2.0*(log(fabs(f0))+12.0*log(10.0))/cinf;
	temp=0.0;
	Lamb=Lamb/100.0;
	intg(0.0, Lamb, charact_funct2, 1e-14, 1e-10, &tp, &abserr);
		temp += tp;
	for(i=1; i<101;i++)
	{
		intg(i*Lamb, (i+1)*Lamb, charact_funct2, 1e-14, 1e-10, &tp, &abserr);
		temp += tp;

	}

      intg(0.0, Lamb, charact_funct2, 1e-14, 1e-10, &tp, &abserr);

      tp=0.5+temp/M_PI;
      return tp;
    }
}

int FT_Call_Heston(double s, double strike, double t, double ri, double dividi, double sigma0,double ka,double theta,double sigma2,double rhow,double *ptprice, double *ptdelta)
{
  double proba1,proba2,temp;

  K=strike;
  S=s;
  T=t;
  sigma=sigma2;
  v=sigma0;
  teta=theta;
  lambda=0.;
  r=ri;
  divid=dividi;
  rho=rhow;
  k=ka;

  proba1=probabilities(1);
  proba2=probabilities(2);
  
  temp=s*proba1*exp(-divid*t);
  temp-=K*exp(-r*t)*proba2;

  /* Price*/
  *ptprice=temp;
    
  /* Delta */
  *ptdelta=proba1*exp(-divid*t);

  return 0;
}

int FT_Put_Heston(double s, double strike, double t, double ri, double dividi, double sigma0,double ka,double theta,double sigma2,double rhow,double *ptprice, double *ptdelta)
{
  double proba1,proba2,temp;

  K=strike;
  S=s;
  T=t;
  sigma=sigma2;
  v=sigma0;
  teta=theta;
  lambda=0.;
  r=ri;
  divid=dividi;
  rho=rhow;
  k=ka;

  proba1=probabilities(1);
  proba2=probabilities(2);
  
 temp=K*exp(-r*t)*(1.-proba2);
  temp-=s*(1.-proba1)*exp(-divid*t);

   /* Price*/
  *ptprice=temp;
    
  /* Delta */
*ptdelta=-(1.-proba1)*exp(-divid*t);

  return 0;
}

