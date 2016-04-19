#include "svj.h"
#include "hes1d_std.h"
//
static double T,sigma,rho,kappa, V0, r, divid, teta, lambda0,m0,v, S, K;
static int func_type=0;
static int heston=0,merton=0;
static double phi=0.;
static int probadelta = 0;
static int bk,initlog=0;
static fcomplex lk_1;
//
static void init_log(void)
{
  initlog = 0;
  bk = 0;
}
//
static void grad_charact_func(double k,int dimx,double *grad)
{
  double X,tau,roeps,u,b,I,eps,eps2;
  fcomplex Ak,Bk,Ck,Dk,Lambdak,z0,z1,z2,z3,z4,z5,zeta,psi_moins,psi_plus,expo,expo1,ans;
  fcomplex dlk;
  fcomplex dLambdakm,dLambdakv,dz1m,dz1v,dz2m,dz2v,gradV0,gradlambda,gradm,gradv;
  fcomplex gradkappa,gradteta,gradsigma,gradrho,dAkt,dAkk,dAks,dAkr,dBkt,dBkk,dBks,dBkr;
  fcomplex dzetak,dzetas,dzetar,dpsipk,dpsips,dpsipr,dpsimk,dpsims,dpsimr,dz3k,dz3s,dz3r;
  fcomplex dz0k,dz0s,dz0r;
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
	  printf("erreur : dans charact_func il faut initialiser func_type a 1 ou 2.\n");
	  exit(-1);
	}
  //
  if(heston==1)
	{
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
	  // Test du log
	  if(initlog>0)
		{
		  dlk = Csub(Ak,lk_1);
		  if(dlk.i < -PI)
			{
			  bk = bk + 1;
			}
		  else if(dlk.i > PI)
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
	  Ak = Cadd(Ak, Complex(0.,2*PI*bk)); 
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
	}
  else
	{
	  Ak = Complex(0.,0.);
	  Bk = Complex( -0.5*tau*k*k , 0.5*tau*u*k );
	}
  //
  if(merton==1)
	{
	  z1   = Complex( -0.5*v*v*k*k + I*(m0+0.5*v*v) , (m0+I*v*v)*k );
	  z1   = Cexp(z1); 
	  dz1m = Cmul( z1 , Complex(I,k) );
	  dz1v = Cmul( z1 , Complex(-k*k*v+I*v,2.*k*I*v) );
	  
	  z2   = Complex(I,k);
	  dz2m = RCmul( exp(m0+0.5*v*v) , z2);
	  dz2v = RCmul( v*exp(m0+0.5*v*v) , z2);
	  z2   = RCmul( exp(m0+0.5*v*v)  -1, z2);
	  z2   = Cadd( Complex(1.,0.) , z2 );


	  Lambdak   = Csub(z1,z2);
      dLambdakm = Csub(dz1m,dz2m);
	  dLambdakv = Csub(dz1v,dz2v);
	  //
	  // Ck pour le cas d'une intensite non contante. Non encore prevu. 
	  Ck = Complex(0.,0.);
	  Dk = RCmul(tau,Lambdak);
	  //
	  // derivee / V0
	  gradV0 = Bk;
	  // derivee / Lambda
	  gradlambda = Dk;
	  // derivee / m0
	  gradm = RCmul( tau*lambda0 , dLambdakm );
	  // derivee / v
	  gradv = RCmul( tau*lambda0 , dLambdakv );
	  //
	}
  else
	{
	  Ck = Complex(0.,0.);
	  Dk = Complex(0.,0.);
	}
  //
  ans = Cadd( Ak , RCmul(V0,Bk) );
  ans = Cadd( ans , Ck );
  ans = Cadd( ans , RCmul(lambda0,Dk) );  
  ans = Cadd( ans , Complex(0.,k*X) );
  ans = Cexp(ans);
  ans = Cdiv(ans,Complex(0.,k));
  //
  if(merton==0)
	// cas heston seul 
	{
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
  else if(heston==0)
	// cas merton seul
	{
	  gradV0     = Cmul( ans , gradV0 );
	  gradlambda = Cmul( ans , gradlambda );
	  gradm      = Cmul( ans , gradm );
	  gradv      = Cmul( ans , gradv );
	  //
	  grad[0] = gradV0.r;
	  grad[1] = gradlambda.r;
	  grad[2] = gradm.r;
	  grad[3] = gradv.r;
	}
  else
	// cas heston + merton
	{
	  gradV0     = Cmul( ans , gradV0 );
	  gradkappa  = Cmul( ans , gradkappa );
	  gradteta   = Cmul( ans , gradteta );
	  gradsigma  = Cmul( ans , gradsigma );
	  gradrho    = Cmul( ans , gradrho );
	  gradlambda = Cmul( ans , gradlambda );
	  gradm      = Cmul( ans , gradm );
	  gradv      = Cmul( ans , gradv );
	  //
	  grad[0] = gradV0.r;
	  grad[1] = gradkappa.r;
	  grad[2] = gradteta.r;
	  grad[3] = gradsigma.r;
	  grad[4] = gradrho.r;
	  grad[5] = gradlambda.r;
	  grad[6] = gradm.r;
	  grad[7] = gradv.r;	  
	}
  //  return ans.r;
  
}
//
static void grad_charact_func0(double k,int dimx,double *grad)
{
  double X,tau,roeps,u,eps,eps2;
  fcomplex Ak,Bk,Ck,Dk,Lambdak,z0,z1,z2,z3,z4,z5,zeta,psi_moins,psi_plus,expo,expo1,ans;
  fcomplex dlk;
  fcomplex dLambdakm,dLambdakv,dz1m,dz1v,dz2m,dz2v,gradV0,gradlambda,gradm,gradv;
  fcomplex gradkappa,gradteta,gradsigma,gradrho,dAkt,dAkk,dAks,dAkr,dBkt,dBkk,dBks,dBkr;
  fcomplex dzetak,dzetas,dzetar,dpsipk,dpsips,dpsipr,dpsimk,dpsims,dpsimr,dz3k,dz3s,dz3r;
  fcomplex dz0k,dz0s,dz0r;
  //  
  tau   = T;
  eps   = sigma;
  roeps = rho*eps;
  X     = log(S/K) + (r - divid)*tau;
  //
  u = kappa - roeps/2.;
  //
  eps2 = eps*eps;
  //
  if(heston==1)
	{
	  //------------------------------------
	  z0   = Complex(u,roeps*k);
	  dz0k = Complex(1.,0.);
	  dz0s = Complex(-0.5*rho,rho*k);
	  dz0r = Complex(-0.5*eps,eps*k);
	  //------------------------------------
	  zeta.r = k*k*eps2*(1.-rho*rho) + u*u + eps2*0.25;
	  zeta.i = 2.*k*roeps*u;
	  zeta   = Csqrt(zeta);
	  //
	  dzetak = z0;
	  dzetak = Cdiv( z0 , zeta );
	  dzetas = Complex( k*k*eps*(1.-rho*rho) - 0.5*rho*u +eps*0.25 , k*rho*u  );
	  dzetas = Cdiv( dzetas , zeta );
	  dzetar = Complex( -rho*k*k*eps2 - 0.5*eps*u , k*eps*u  -0.5*eps*k*rho*eps );
	  dzetar = Cdiv( dzetar , zeta );
	  //------------------------------------
	  psi_moins = z0;
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
	  z1 = Complex( -(k*k+0.25) , 0. );
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
	  //
	  Ak = Cdiv( z3 , RCmul(2.,zeta) );
	  Ak = Clog(Ak);
	  // Test du log
	  if(initlog>0)
		{
		  dlk = Csub(Ak,lk_1);
		  if(dlk.i < -PI)
			{
			  bk = bk + 1;
			}
		  else if(dlk.i > PI)
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
	  Ak = Cadd(Ak, Complex(0.,2*PI*bk));
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
	}
  else
	{
	  Ak = Complex(0.,0.);
	  Bk = Complex( -0.5*tau*(k*k+0.25) ,0.);
	}
  //
  if(merton==1)
	{
	  z1 = Complex( 0.5*m0-0.5*v*v*(k*k-0.25) , -k*(m0+0.5*v*v) );
	  z1 = Cexp(z1);
	  dz1m = Cmul( z1 , Complex(0.5,-k) );
	  dz1v = Cmul( z1 , Complex(-(k*k-0.25)*v,-k*v) );

	  z2 = Complex(0.5,-k);
	  dz2m = RCmul( exp(m0+0.5*v*v) , z2);
	  dz2v = RCmul( v*exp(m0+0.5*v*v) , z2);
	  z2 = RCmul( exp(m0+0.5*v*v) - 1. , z2);
	  z2 = Cadd( Complex(1.,0.) , z2 );
	  	  
	  Lambdak   = Csub(z1,z2);
	  dLambdakm = Csub(dz1m,dz2m);
	  dLambdakv = Csub(dz1v,dz2v);
	  //
	  // Ck pour le cas d'une intensite non contante. Non encore prevu.
	  Ck = Complex(0.,0.);
	  Dk = RCmul(tau,Lambdak);
	  //
	  // derivee / V0
	  gradV0 = Bk;
	  // derivee / Lambda
	  gradlambda = Dk;
	  // derivee / m0
	  gradm = RCmul( tau*lambda0 , dLambdakm );
	  // derivee / v
	  gradv = RCmul( tau*lambda0 , dLambdakv );
	  //
	}
  else
	{
	  Ck = Complex(0.,0.);
	  Dk = Complex(0.,0.);
	}
  //  
  ans = Cadd( Ak , RCmul(V0,Bk) );
  ans = Cadd( ans , Ck );
  ans = Cadd( ans , RCmul(lambda0,Dk) );
  ans = Cadd( ans , RCmul(X,Complex(0.5,-k) ) );
  ans = Cexp(ans);
  ans = Cdiv(ans,Complex(k*k+0.25,0.));
  //
  /*
  if(probadelta == 1)
	{
	  ans = Cmul( ans , Complex(0.5,-k) );
	  ans = RCmul( 1./S , ans );
	}
  */
  //
  if(merton==0)
	// cas heston seul
	{
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
  else if(heston==0)
	// cas merton seul
	{
	  //
	  gradV0     = Cmul( ans , gradV0 );
	  gradlambda = Cmul( ans , gradlambda );
	  gradm      = Cmul( ans , gradm );
	  gradv      = Cmul( ans , gradv );
	  //
	  grad[0] = gradV0.r;
	  grad[1] = gradlambda.r;
	  grad[2] = gradm.r;
	  grad[3] = gradv.r;
	}
  else
	// cas heston + merton
	{
	  gradV0     = Cmul( ans , gradV0 );
	  gradkappa  = Cmul( ans , gradkappa );
	  gradteta   = Cmul( ans , gradteta );
	  gradsigma  = Cmul( ans , gradsigma );
	  gradrho    = Cmul( ans , gradrho );
	  gradlambda = Cmul( ans , gradlambda );
	  gradm      = Cmul( ans , gradm );
	  gradv      = Cmul( ans , gradv );
	  //
	  grad[0] = gradV0.r;
	  grad[1] = gradkappa.r;
	  grad[2] = gradteta.r;
	  grad[3] = gradsigma.r;
	  grad[4] = gradrho.r;
	  grad[5] = gradlambda.r;
	  grad[6] = gradm.r;
	  grad[7] = gradv.r;	  
	}
  //  return ans.r;
}
//
static void  grad_probabilities(int n,int dimx,double *grad)
{
  double *tp,*tpi,tpimax;
  int i,m,ngauss=10,N=100;
  double a,b,h=10.,tol=1.e-12;
  //  
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
  //
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
		tp[m] = tp[m] + tpi[m];
	  // on cherche le max des increments
	  tpimax = fabs(tpi[0]);
	  for(m=1;m<dimx;m++)
		tpimax = MAX(tpimax,fabs(tpi[m]));
	  //
	  i++;
	}
  //
  for(m=0;m<dimx;m++)
	grad[m] =  phi*tp[m]/PI;
  //================
  //
  free(tp);
  free(tpi);
  free_gauss();
  //
  //================
  //
  
}
//
static void  grad_probabilities2(int dimx,double *grad)
{
  double *tp,*tpi,tpimax;
  int i,m,ngauss=10,N=1000;
  double a,b,h=1.,tol=1.e-12;
  //
  init_gauss(ngauss);
  //
  tp  =  (double*)malloc(sizeof(double)*dimx);
  tpi =  (double*)malloc(sizeof(double)*dimx);
  //
  //================
  for(m=0;m<dimx;m++) 
	tp[m]  = 0.;
  tpimax = 1.;
  i   = 0;
  while( (tpimax>tol) && (i<N) )
	{
	  a    = i*h;
	  b    = a + h;
	  integrale_gauss_vect(grad_charact_func0,a,b,dimx,tpi);
	  for(m=0;m<dimx;m++)  
		tp[m] = tp[m] + tpi[m];
	  // on cherche le max des increments
	  tpimax = fabs(tpi[0]);
	  for(m=1;m<dimx;m++)
		tpimax = MAX(tpimax,fabs(tpi[m]));
	  //
	  i++;
	}
  //
  for(m=0;m<dimx;m++)
	grad[m] = K*tp[m]/PI;
  //================
  //
  free(tp);
  free(tpi);
  free_gauss();
  //
  //================
  //
}
//
int calc_grad_svj(SVJPARAMS *svj,int dimx, double *grad)
{
  int i;
  double *grad1,*grad2;
  //
  grad1 =  (double*)malloc(sizeof(double)*dimx);
  grad2 =  (double*)malloc(sizeof(double)*dimx);
  //
  K=svj->K;
  S=svj->St0;
  T=svj->T;
  sigma=svj->sigmav;
  V0=svj->V0;
  teta=svj->theta;
  r=svj->r;
  divid=svj->divid;
  rho=svj->rho;
  kappa=svj->kappa;
  lambda0 = svj->lambda;
  m0 = svj->m0;
  v  = svj->v;
  phi=svj->phi;
  heston = svj->heston;
  merton = svj->merton;
  //
  if(svj->type_f==1)
	{
	  init_log();
	  grad_probabilities(1,dimx,grad1);
	  init_log();
	  grad_probabilities(2,dimx,grad2);
	  //
	  for(i=0;i<dimx;i++)
		grad[i] = phi*(S*grad1[i]*exp(-divid*T) - K*exp(-r*T)*grad2[i]);
	}
  else if(svj->type_f==2)
	{	
	  init_log();
	  probadelta = 0;
	  grad_probabilities2(dimx,grad1);
	  //
	  for(i=0;i<dimx;i++)
		grad[i] = -grad1[i]*exp(-r*T);
	}
  else
	{
	  printf("Erreur dans svj.c : parametre svj->type_f inconnu.\n");
	  exit(-1);
	}
  //
  free(grad1);
  free(grad2);
  //
  return OK;
}
/*
int CALC(CF_CallHeston)(void *Opt, void *Mod, PricingMethod *Met)
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

  return CFCallHeston(ptMod->S0.Val.V_PDOUBLE,
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



int CHK_OPT(CF_CallHeston)(void *Opt, void *Mod)
{ 
return strcmp( ((Option*)Opt)->Name,"CallEuro");   
}



static int MET(Init)(PricingMethod *Met)
{
  return OK;
}

PricingMethod MET(CF_CallHeston)=
{
  "CF_Heston",
  {{" ",END,0,FORBID}},
  CALC(CF_CallHeston),
  {{"Price",DOUBLE,100,FORBID},
   {"Delta",DOUBLE,100,FORBID} ,
   {" ",END,0,FORBID}},
  CHK_OPT(CF_CallHeston),
  CHK_ok,
  MET(Init)
};
*/
