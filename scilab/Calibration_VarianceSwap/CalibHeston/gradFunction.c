
#include <stdlib.h>
#include "pnl/pnl_complex.h"
#include "pnl/pnl_mathtools.h"
#include "gradFunction.h"
#include "utils/my_integral.h"
#include "bsvanillas.h"


static int initOk=0,TypeNorme,TypeModel;
static int LogNorme=0;
static int MelangeNormes=0;
static int sigmaimpok=0,sigmaimpko=0;
static DataMarket *DM;

static double T,sigma,rho,kappa, V0, r, divid, teta, S, K;
static int func_type=0;
static double phi=0.;
static int bk,initlog=0;
static dcomplex lk_1;

extern int FT_Call_Heston(double St0, double strike, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,double *ptprice, double *ptdelta);
extern int FT_Put_Heston(double St0, double strike, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,double *ptprice, double *ptdelta);

int calc_grad_heston(int ifcall, double K, double St0, double T, double r, double divid, double sigmav, double V0, double theta, double rho, double kappa, int dimx, double *grad);
int Grad_FT_Call_Heston(double St0, double strike, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,int dimx,double *grad);
int Grad_FT_Put_Heston(double St0, double strike, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,int dimx,double *grad);

//===================================== GRAD FUNCTION ==========================
// Computes gradient of the cost function

void  gradFunction(int dimx,double *x, double *grad)
{
  
  int i,n;
  double fx=0,fxi;
  int ibidon;
  double V0,kappa,theta,sigmav,rho,prix,delta;
  double St0,K,T,prixobs,sigma_imp_obs;
  double *gradi;
  double ri,di;
  
  gradi = (double*)malloc(dimx*sizeof(double));
  
  for(n=0;n<dimx;n++){gradi[n] =0.;  grad[n] =0.;}
  //
  	  if(initOk != 1) {
	    printf("Erreur, avant d'appeler gradFunction, il faut initialiser les paramaetres.\n ");	
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
	      ri =DM->r[i];
	      di =DM->d[i];
	      //
	      if(DM->TypeOpt[i] == 1) {
		ibidon=FT_Call_Heston(St0, K, T, ri, di, V0, kappa, theta, sigmav, rho, &prix, &delta);
		ibidon=Grad_FT_Call_Heston(St0, K, T, ri, di, V0, kappa, theta, sigmav, rho, dimx, gradi);
	      } else {
		ibidon=FT_Put_Heston(St0, K, T, ri, di, V0, kappa, theta, sigmav, rho, &prix, &delta);
		ibidon=Grad_FT_Put_Heston(St0, K, T, ri, di, V0, kappa, theta, sigmav, rho, dimx, gradi);
	      }
	      //
		if(prix<0.0) {prix = 1e-14;}
	      prixobs = DM->PrixObs[i];
	      sigma_imp_obs = DM->SigmaImpObs[i];
	      fxi = dNormeMarket(TypeNorme, DM->TypeOpt[i], St0, T, K, ri, di, prix, prixobs, sigma_imp_obs, dimx, gradi);
	      //
	      fx = fx + fxi*DM->wi[i];
	      for(n=0;n<dimx;n++) grad[n] = grad[n] + gradi[n]*DM->wi[i];
	    }
}

void initParamsGrad(int nbdata,DataMarket *_DM,int _TypeNorme,int _TypeModel,int _LogNorme,int _MelangeNormes)
{
  
  int i;
  
  if(DM == NULL)
    DM = CreateDataMarket(nbdata);
  
  for(i=0;i<nbdata;i++)
    {
      DM->nbdata     = _DM->nbdata;
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
  
  TypeNorme     = _TypeNorme;
  TypeModel     = _TypeModel;
  LogNorme      = _LogNorme;
  MelangeNormes = _MelangeNormes;
  
  initOk        = 1;
  
}

//=============================== D NORME MARKET ==========================
// 
double dNormeMarket(int TypeNorme,int TypeOpt,double St0,double T,double K,double r,double divid,double prix,double prixobs,double sigma_imp_obs,int dimx,double *gradi)
{
  //
  int j1,i,ibidon;
  double fx,sigma_imp,dprice,sigma_imp0,dsigma_imp;
  //
  sigma_imp     = SigmaImplicite(TypeOpt,St0,T,K,r,divid,prix,&j1);
      
  // Case 1 : (absolute error)^2
  if(TypeNorme == 1)
    {
      fx     = (prix - prixobs);
	  for(i=0;i<dimx;i++)
	    gradi[i] = gradi[i] * 2.* fx;
	  fx     = fx * fx;
    }
  // Case 2 : (relative error)^2
  else if(TypeNorme == 2)
    {
      fx     = (prix - prixobs)/prixobs;
      for(i=0;i<dimx;i++)
	gradi[i] = gradi[i] * 2.* fx/prixobs;
      fx     = fx * fx;
    }
  // Case 3 : (absolute error in SigmaImplicites)^2
  else if(TypeNorme == 3)
    {
      //
      if(j1==-1) {
	sigmaimpko++;
	  fx = dNormeMarket(2,TypeOpt,St0,T,K,r,divid,prix,prixobs,sigma_imp_obs,dimx,gradi);
      } 
	else {
	sigmaimpok++;
	if(TypeOpt==1) {
	  ibidon =  dCall_BlackScholes_73(St0,K,T,r,divid,sigma_imp,&dprice);
	} else {
	  ibidon =  dPut_BlackScholes_73(St0,K,T,r,divid,sigma_imp,&dprice);
	}
	
	if(dprice==0.) {
	  dprice = 1.e-3;
	}
	// pour empecher dprice d'exploser.
	if(dprice<-1.e8) { dprice=-1.e8;}
	else if(dprice>1.e8) { dprice=1.e8;}		
	//
	dsigma_imp = 1. / dprice ;
	//
	sigma_imp0 = 1.e-2;
	//
	fx     = (sigma_imp - sigma_imp_obs)/sigma_imp0;
	for(i=0;i<dimx;i++)
	  gradi[i] = gradi[i] * 2.* fx * dsigma_imp/sigma_imp0;
	fx     = fx * fx;
      }
    } 
  // Case 4: logarithmic norm
  else if(TypeNorme == 4) 
  {
	fx     = (prix - prixobs)/prix + log(prix/prixobs);
	  for(i=0;i<dimx;i++)
	    gradi[i] = gradi[i] * fx;
  }
  else {
      printf("Error in dNormeMarket, TypeNorme = %d is unknown\n",TypeNorme);  
      exit(-1);
    }
  //
  return fx;
}
//
void FreeGradFunction()
{
  FreeDataMarket(DM);
  DM = NULL;
}
//
//================================VANILLAS GRAD=================

int Grad_FT_Call_Heston(double St0, double strike, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,int dimx,double *grad)
{
  double K;
  int ifcall=1;
   //
  K=strike;
 
  calc_grad_heston(ifcall, K, St0, T, r, divid, sigmav, V0, theta, rho, kappa,dimx,grad);
  //  
  return 0;
}

int Grad_FT_Put_Heston(double St0, double strike, double T, double r, double divid, double V0,double kappa,double theta,double sigmav,double rho,int dimx,double *grad)
{
  double K;
int ifcall=-1;

 
  K=strike;
  
  calc_grad_heston(ifcall, K, St0, T, r, divid, sigmav, V0, theta, rho, kappa,dimx,grad);
  
  return 0;
}

//=======calc grad===================

static void init_log(void)
{
  initlog = 0;
  bk = 0;
}
//
static void grad_charact_func(double k,int dimx,double *grad)
{
  double X,tau,roeps,u,b,I,eps,eps2;
  dcomplex Ak,Bk,z0,z1,z2,z3,z4,z5,zeta,psi_moins,psi_plus,expo,expo1, ans;
  dcomplex dlk;
  dcomplex gradV0;
  dcomplex gradkappa,gradteta,gradsigma,gradrho,dAkt,dAkk,dAks,dAkr,dBkt,dBkk,dBks,dBkr;
  dcomplex dzetak,dzetas,dzetar,dpsipk,dpsips,dpsipr,dpsimk,dpsims,dpsimr,dz3k,dz3s,dz3r;
  dcomplex dz0k,dz0s,dz0r;
  //
  tau   = T;
  eps   = sigma;
  roeps = rho*eps;
  X     = log(S/K) + (r - divid)*tau; 
  eps2  = eps*eps;
  //
  if(func_type==1)
	{
	  u = 1.;
	  b = kappa - roeps;
	  I = 1.;
	}
  else if(func_type==2)
	{
	  u = -1.;
	  b = kappa;
	  I =  0.;
	}
  else
	{
	  printf("error in charact_func: func_type is not initialized.\n");
	  exit(-1);
	}
	  //------------------------------------
	  z1   = Complex(k*k,-u*k);
	  z0   = Complex(b,-roeps*k);
	  z2   = Cmul(z0,z0);
	  //
	  dz0k = Complex(1.,0.);
	  dz0s = Complex(-rho*I,-rho*k);
	  dz0r = Complex(-eps*I,-eps*k);
	  //------------------------------------
	  zeta = Cadd(z2,RCmul(eps2,z1));
	  zeta = Csqrt(zeta);
	  //
	  dzetak = Cdiv( z0 , zeta );
	  dzetas = Cmul( dz0s , z0 );
	  dzetas = Cadd( dzetas , RCmul(eps,z1) );
	  dzetas = Cdiv( dzetas , zeta );
	  dzetar = Cmul( dz0r , z0 );
	  dzetar = Cdiv( dzetar , zeta );
	  //------------------------------------
	  psi_moins = Complex(b,-roeps*k);
	  psi_plus  = RCmul(-1.,psi_moins);
	  psi_moins = Cadd(psi_moins,zeta);
	  psi_plus  = Cadd(psi_plus,zeta);
	  //
	  dpsimk = dz0k;
	  dpsims = dz0s;
	  dpsimr = dz0r;
	  dpsipk = RCmul(-1.,dpsimk);
	  dpsips = RCmul(-1.,dpsims);
	  dpsipr = RCmul(-1.,dpsimr);
	  //
	  dpsimk = Cadd(dpsimk,dzetak);
	  dpsims = Cadd(dpsims,dzetas);
	  dpsimr = Cadd(dpsimr,dzetar);
	  //
	  dpsipk = Cadd(dpsipk,dzetak);
	  dpsips = Cadd(dpsips,dzetas);
	  dpsipr = Cadd(dpsipr,dzetar);
	  //------------------------------------
	  expo = Cexp( RCmul(-tau,zeta) );
	  expo1= Cmul(psi_plus,expo);
	  z3   = Cadd( psi_moins , expo1 );
	  expo1= RCmul(-tau,expo1);
	  //
	  dz3k = Cadd( dpsimk , Cmul(dpsipk,expo) );
	  dz3k = Cadd( dz3k , Cmul(expo1,dzetak) );
	  dz3s = Cadd( dpsims , Cmul(dpsips,expo) );
	  dz3s = Cadd( dz3s , Cmul(expo1,dzetas) );
	  dz3r = Cadd( dpsimr , Cmul(dpsipr,expo) );
	  dz3r = Cadd( dz3r , Cmul(expo1,dzetar) );
	  //------------------------------------
	  z1 = RCmul(-1.,z1);
	  z4 = Csub(Complex(1.,0),expo);
	  Bk = Cmul( z1 , z4 );
	  Bk = Cdiv(Bk,z3);
	  expo1 = RCmul(tau,expo);
	  expo1 = Cmul(expo1,z3);
	  //
	  dBkk  = Csub( Cmul(expo1,dzetak) , Cmul(z4,dz3k) );
	  dBkk  = Cdiv(dBkk,Cmul(z3,z3));
	  dBkk  = Cmul(dBkk,z1);
	  dBks  = Csub( Cmul(expo1,dzetas) , Cmul(z4,dz3s) );
	  dBks  = Cdiv(dBks,Cmul(z3,z3));
	  dBks  = Cmul(dBks,z1);
	  dBkr  = Csub( Cmul(expo1,dzetar) , Cmul(z4,dz3r) );
	  dBkr  = Cdiv(dBkr,Cmul(z3,z3));
	  dBkr  = Cmul(dBkr,z1);
	  dBkt  = Complex(0.,0.);
	  //------------------------------------
	  Ak = Cdiv( z3 , RCmul(2.,zeta) );
	  Ak = Clog(Ak);
	  // Test log
	  if(initlog>0)
		{
		  dlk = Csub(Ak,lk_1);
		  if(dlk.i < -M_PI)
			{
			  bk = bk + 1;
			}
		  else if(dlk.i > M_PI)
			{
			  bk = bk - 1;
			}
		  initlog++;
		  lk_1 = Ak;
		} else {
		  initlog++;
		  lk_1 = Ak;
		}
	  //
	  Ak = Cadd(Ak, Complex(0.,2*M_PI*bk)); 
	  //
	  Ak = RCmul( 2. , Ak );
	  Ak = Cadd( RCmul(tau,psi_plus) , Ak);
	  z4 = Ak;
	  Ak = RCmul( -kappa*teta/eps2 , Ak);
	  //
	  z5 = Cmul(zeta,z3);
	  z5 = Cdiv(Complex(2.,0),z5);
	  //
	  dAkk = Csub( Cmul(dz3k,zeta) , Cmul(z3,dzetak) );
      dAks = Csub( Cmul(dz3s,zeta) , Cmul(z3,dzetas) );
	  dAkr = Csub( Cmul(dz3r,zeta) , Cmul(z3,dzetar) );
	  //
	  dAkk = Cmul(dAkk,z5);
	  dAks = Cmul(dAks,z5);
	  dAkr = Cmul(dAkr,z5);
	  //
	  dAkk = Cadd(dAkk,RCmul(tau,dpsipk));
	  dAks = Cadd(dAks,RCmul(tau,dpsips));
	  dAkr = Cadd(dAkr,RCmul(tau,dpsipr));
	  //
	  dAkk = RCmul( -kappa*teta/eps2 , dAkk );
	  dAks = RCmul( -kappa*teta/eps2 , dAks );
	  dAkr = RCmul( -kappa*teta/eps2 , dAkr );
	  dAkk = Cadd( dAkk , RCmul( -teta/eps2 , z4) );
	  dAks = Cadd( dAks , RCmul( 2.*kappa*teta/(eps2*eps) , z4) );
	  dAkt = RCmul( -kappa/eps2 , z4);
	  //------------------------------------
	  // derivee / V0
	  gradV0 = Bk;
	  // derivee / kappa
	  gradkappa = Cadd( dAkk , RCmul(V0,dBkk) );
	  // derivee / teta
	  gradteta = Cadd( dAkt , RCmul(V0,dBkt) );
	  // derivee / sigmav
	  gradsigma = Cadd( dAks , RCmul(V0,dBks) );
	  // derivee / rho
	  gradrho = Cadd( dAkr , RCmul(V0,dBkr) );
	  //
	
  
  //
  ans = Cadd( Ak , RCmul(V0,Bk) );
  ans = Cadd( ans , Complex(0.,k*X) );
  ans = Cexp(ans);
  ans = Cdiv(ans,Complex(0.,k));
  //
 
	  gradV0     = Cmul( ans , gradV0 );
	  gradkappa  = Cmul( ans , gradkappa );
	  gradteta   = Cmul( ans , gradteta );
	  gradsigma  = Cmul( ans , gradsigma );
	  gradrho    = Cmul( ans , gradrho );
	  //
	  grad[0] = gradV0.r;
	  grad[1] = gradkappa.r;
	  grad[2] = gradteta.r;
	  grad[3] = gradsigma.r;
	  grad[4] = gradrho.r;

}

static void  grad_probabilities(int n,int dimx,double *grad)
{
  double *tp,*tpi,tpimax;
  int i,m,ngauss=10,N=100;
  double a,b,h=10.,tol=1.e-12;

  init_gauss(ngauss);
  //
  if(n==1)
    {
	  func_type = 1;
	}
  else
	{
	  func_type = 2;
	}
  //
  //================
  tp  =  (double*)malloc(sizeof(double)*dimx);
  tpi =  (double*)malloc(sizeof(double)*dimx);
  //================
  for(m=0;m<dimx;m++)
	tp[m]  = 0.;
  tpimax = 1.;
  i   = 0;
  while( (tpimax>tol) && (i<N) )
	{
	  a    = i*h;
	  b    = a + h;
      integrale_gauss_vect(grad_charact_func,a,b,dimx,tpi);
	  for(m=0;m<dimx;m++)
	{
		tp[m] = tp[m] + tpi[m];
	}
	  // on cherche le max des increments
	  tpimax = fabs(tpi[0]);
	  for(m=1;m<dimx;m++)
		tpimax = MAX(tpimax,fabs(tpi[m]));
	  //
	  i++;
	}
  //
  for(m=0;m<dimx;m++)
	grad[m] =  phi*tp[m]/M_PI;
  //================
  //
  free(tp);
  free(tpi);
  free_gauss();
}

//========================================== grad heston ================

int calc_grad_heston(int ifcall, double strike, double St0, double matur, double r0, double divid0, double sigmav0, double vv0, double theta, double rho0, double kappa0, int dimx, double *grad)
{
  int i;
  double *grad1,*grad2;
  //
  grad1 =  (double*)malloc(sizeof(double)*dimx);
  grad2 =  (double*)malloc(sizeof(double)*dimx);
  //
  K=strike;
  S=St0;
  T=matur;
  sigma=sigmav0;
  V0=vv0;
  teta=theta;
  r=r0;
  divid=divid0;
  rho=rho0;
  kappa=kappa0;
 
  phi=ifcall;
   //
 
	  init_log();
	  grad_probabilities(1,dimx,grad1);
	  init_log();
	  grad_probabilities(2,dimx,grad2);
	  //
	  for(i=0;i<dimx;i++)
		grad[i] = phi*(S*grad1[i]*exp(-divid*T) - K*exp(-r*T)*grad2[i]);

  free(grad1);
  free(grad2);
  //
  return 0;
}
