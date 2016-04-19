#include "lmm_cev_pricer.h"

#define PI 3.14159
#define K 41/* MAX NUMBER OF MATURITIES + 1 */

const double maturity_structure[K]={0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20};
const double tenor=0.5;
const double epsilon=20;/* factor of LCEV model*/


//int_pts= MC or number of discretiz fore auxiliary functions
void cev_price(size_t type_model,size_t type_product,size_t type_pricing, 
		 size_t type_scheme, double alpha, double fzero, double H,
		 double expiry,double swp_lenght,int int_pts,double* out){
  
  double prix,implied,err; 
  int maturity=(int)expiry*2;
  int swptn_num_of_paym=(int)swp_lenght*2;
  
 
  //CEV Model---------------------------
  if(type_model==0){
	if(type_product==0){
	  if(type_pricing==0){
	    out[0]=10000*prix_caplet(0,maturity,alpha,H,fzero,int_pts); 
	    out[1]=0.0;
	  }
	  else{
	    srand(time(NULL));
	    
	    double delta_i=0.25;/* discretisation for the simulation */
	    int L=(K-1)*(int)(tenor/delta_i);/* dimension of the subdivision */
	    double t[L+1];/* vector of discretisation */
	    double v=0;
	    double u=0/* price */,g=0/* standard deviation */,a=0,b=0;
	    double c[K];
	    int i=0,j=0,r=0;
	    for(i=0;i<K;++i){c[i]=0.0;}
	    for(i=0;i<=L;i++) t[i]=i*delta_i;
	    for(r=0;r<K;r++) c[r]=0;
	    
	    
	    for(j=0;j<int_pts;j++){
	      /*time evolution, c contains F_{0,1}(T_0),F_{1,2}(T_1)...*/
	      function_F(maturity+1,t,c,fzero,delta_i,alpha,type_model,type_scheme);
	      
	      if(c[maturity]-H>0){
		b=function_delta(maturity)*(c[maturity]-H);
		a=b/function_B_en_maturite(maturity_structure[maturity+1],c);
	      }
	      else{
		b=0;
		a=0;
	      }
	      
	      v=v+a; 
	      g=g+a*a;
	     
	    }
	    
	    u=v/int_pts;
	    out[0]=u*10000;
	    out[1]=sqrt(g/int_pts-u*u)/sqrt(int_pts)*10000;
		
	  }
	}
   
	else{
	  if(type_pricing==0){
	    double vs,black;
	    double t=0;
	     vs=somme_integral(0,maturity_structure[maturity],maturity,swptn_num_of_paym,0,alpha,fzero,int_pts);
	     out[0]=10000*prix_swaption(t,0.06,vs,alpha,maturity,swptn_num_of_paym,fzero);
	     out[1]=0;
	  }
	  else{
	  
	    srand(time(NULL));
	    double delta_i=0.125;// discretisation for the simulation 
	    int L=(K-1)*(int)(tenor/delta_i);// dimension of the subdivision 
	    double t[L+1];
	    double v=0,u=0,g=0,a,b;
	    double fts[K];
	    int i,j,l,m,h;
	    double coupon[K];
	    double fktk[K];
	    double teta=0;

	    for(i=0;i<K;++i){fts[i]=0.0;fktk[i]=0.0;coupon[i]=0.;}

	    for(i=0;i<L+1;i++) t[i]=i*delta_i;
	    
	    for(j=0;j<int_pts;j++){
	      b=0;
	      a=0;
	      function_F_swap(maturity,swptn_num_of_paym,t,fts,fzero,delta_i,alpha,type_model,type_scheme,fktk);
	      
	      coupon[maturity]=1;
	      for(l=maturity+1;l<swptn_num_of_paym+1;l++){
		coupon[l]=coupon[l-1]/(1+function_delta(l-1)*fts[l-1]);
	      }
	      for(h=maturity;h<swptn_num_of_paym;h++) b+=function_delta(h)*coupon[h+1];
	      b=b*teta;
	      b=1-coupon[swptn_num_of_paym]-b;
	      if(b<=0){a=0;}
	      if (function_B_en_maturite(maturity_structure[maturity],fktk)==0) {printf("ERROR division by zero");return;}
	      else {a=b/function_B_en_maturite(maturity_structure[maturity],fktk);}
	      v=v+a;
	      g=g+a*a;
	    }
	    u=v/int_pts;
	    out[0]=u*10000;
	    out[1]=sqrt(g/int_pts-u*u)/sqrt(int_pts)*10000;	  
	  }
	}
  }
  else{
    if(type_product==0){ 
      if(type_pricing==0) {
	printf("Sorry, no closed form for LCEV caplets\n going on with MC\n");
	type_pricing=1;
      }
      else{}	
      srand(time(NULL));
      
      double delta_i=0.25;/* discretisation for the simulation */
      int L=(K-1)*(int)(tenor/delta_i);/* dimension of the subdivision */
      
      double t[L+1];/* vector of discretisation */
      double v=0;
      double u=0/* price */,g=0/* standard deviation */,a=0,b=0;
      double c[K];
      int i=0,j=0,r=0;
      
      for(i=0;i<K;++i){c[i]=0.0;}
      for(i=0;i<=L;i++) t[i]=i*delta_i;
      for(r=0;r<K;r++) c[r]=0;
      for(j=0;j<int_pts;j++){
	/*time evolution, c contains F_{0,1}(T_0),F_{1,2}(T_1)...*/
	function_F(maturity+1,t,c,fzero,delta_i,alpha,type_model,type_scheme);
	
	if(c[maturity]-H>0){
	  b=function_delta(maturity)*(c[maturity]-H);
	  a=b/function_B_en_maturite(maturity_structure[maturity+1],c);
	}
	else{b=0;a=0;}
	v=v+a;
	g=g+a*a;
      }
      u=v/int_pts;
      out[0]=u*10000;
      out[1]=sqrt(g/int_pts-u*u)/sqrt(int_pts)*10000;
      
      
    }
    
    
    else{
      if(type_pricing==0){
	double vs,black;
	double t;
	t=0;
	vs=somme_integral(0,maturity_structure[maturity],maturity,swptn_num_of_paym,type_model,alpha,fzero,int_pts);
	out[0]=10000*prix_swaption(t,0.06,vs,alpha,maturity,swptn_num_of_paym,fzero);
	out[1]=0;
      }
      else{
	srand(time(NULL));
	double delta_i=0.125;/* discretisation for the simulation */
	int L=(K-1)*(int)(tenor/delta_i);/* dimension of the subdivision */
	
	double t[L+1];
	double v;
	double u,g,a,b;
	double fts[K];
	int i,j,l,m,h;
	double coupon[K];
	double fktk[K];
	double teta=0.06;
	for(i=0;i<K;++i){fts[i]=0.0;fktk[i]=0.0;}
	for(i=0;i<L+1;i++)  {t[i]=i*delta_i;}
	v=0;
	g=0;
	
	for(j=0;j<int_pts;j++){
	  b=0;
	  a=0;
	  function_F_swap(maturity,swptn_num_of_paym,t,fts,fzero,delta_i,alpha,type_model,type_scheme,fktk);
	  coupon[maturity]=1;
	  for(l=maturity+1;l<swptn_num_of_paym+1;l++){coupon[l]=coupon[l-1]/(1+function_delta(l-1)*fts[l-1]);}
	  
	  for(h=maturity;h<swptn_num_of_paym;h++)	b+=function_delta(h)*coupon[h+1];
	  b=b*teta;
	  b=1-coupon[swptn_num_of_paym]-b;
	  if(b<=0){a=0;}
	  else {
	    a=b/function_B_en_maturite(maturity_structure[maturity],fktk);
	  }
	  v=v+a;
	  g=g+a*a;
	}
	u=v/int_pts;
	out[0]=u*10000;
	out[1]=sqrt(g/int_pts-u*u)/sqrt(int_pts)*10000;
      }
    }
  }
  
  return;
}


//-----------------------body of auxiliary functions --------------------------
int function_n(double t)
{
  int i,h;

  for(i=1;i<K;i++)
	{

	  if((t<maturity_structure[i]) && (t>=maturity_structure[i-1]))
		h=i;
	}

  return h;
}


/* simulation of uniform variable */
double u_random()
{
  double r=0;
  r=(double)rand() / ((double)RAND_MAX + 1);
  return r;
}



/* simulation of independant Gaussian variable */

double simul_normale()
{
  double u1,u2,n;
  u1=u_random();
  u2=u_random();
  n=sqrt(-2*log(u1))*cos(2*PI*u2);
  return n;
}



/* function_fi  programmed for LCEV model */

double function_fi(double x,int i,double alpha)
{
  double a;
  if(i==0)/* CEV model */{
	a=pow(x,alpha);
  }


  if(i==1)
	{
	  if(epsilon<pow(x,alpha-1))
		a=x*epsilon;
	  else{
		a=pow(x,alpha);
	  }
	}
  return a;
  
}


double function_delta( int k)
{
  double delta;
  delta=maturity_structure[k+1]-maturity_structure[k];

  return delta;
}

double function_lambda_j(double t,int k)
{
  double a;

  a=0.05;

  return a;
}



double function_lambda_k_carre(double t){
  double res=0;
  res=0.0025;
  
  return res;
}

double price_zero_coupon_bond(double t,int k,double fzero)
{
  int i;
  double c;
  c=1/(1+fzero*maturity_structure[1]);
  if(k==0)
	return c;
  else{
    for(i=1;i<k;i++)
	  c=c/(1+function_delta(i)*fzero);
    return c;
  }
  
 
  
}

double function_B_s(double t,int debut,int fin,double fzero)
{
  double b=0;
  int i;
  
  for(i=debut;i<fin;i++){
	 b=b+function_delta(i)*price_zero_coupon_bond(t,i+1,fzero);
  }
  
  return b;
}

double function_R(double t,int debut,int fin,double fzero)
{
  double r=0;

  r=(price_zero_coupon_bond(t,debut,fzero)-price_zero_coupon_bond(t,fin,fzero))/function_B_s(t,debut,fin,fzero);

  return r;
  
}


double function_dr_df(double t,int j,int debut,int fin,double fzero)
{
  double a=0,b=0,c=0,d=0;
  int i;
  
  a=function_delta(j)*function_R(t,debut,fin,fzero)/(1+function_delta(j)*fzero);
  b=price_zero_coupon_bond(t,fin,fzero)/(price_zero_coupon_bond(t,debut,fzero)-price_zero_coupon_bond(t,fin,fzero));
  for(i=j;i<fin;i++)
	c=c+function_delta(i)*price_zero_coupon_bond(t,i+1,fzero);
  d=a*(b+(c/function_B_s(t,debut,fin,fzero)));

  return d;
}

double function_wj(double t, int j,double alpha,int debut,int fin,int i,double fzero)
{
  double w=0;

  w=function_dr_df(t,j,debut,fin,fzero)*(function_fi(fzero,i,alpha)/function_fi(function_R(t,debut,fin,fzero),i,alpha));
  return w;
}

double function_somme(double t,int debut,int fin,int i,double alpha,double fzero)
{
  double somme=0,res=0;
  int j;

  for(j=debut;j<fin;j++)
	somme=somme+function_wj(0,j,alpha,debut,fin,i,fzero)*function_lambda_j(t,j);
  res=somme*somme;
  
  return res;
  
}

double somme_integral(double a,double b,int debut,int fin,int i,double alpha,double fzero,int n)
{
  double R=0.0,res,delta;
  int j;

  delta=(b-a)/(double)n;
  for(j=0;j<n;j++){
	R+=function_somme(a+(double)j*delta,debut,fin,i,alpha,fzero);
	
  }
  res=R*delta;

  return res;
}

double calcul_integral(double (*f) (double),double a, double b, int n)
{
  double R=0.0,res,delta;
  int i;

  delta=(b-a)/(double)n;
  for(i=0;i<n;i++){
	R+=f(a+(double)i*delta);
  }
  res=R*delta;

  return res;
}

int function_a(int k, double t,double alpha,double H,int n)
{
  double vk,p,q,a;
  
  vk=calcul_integral(function_lambda_k_carre,t,maturity_structure[k],n);
  p=pow(H,2*(1-alpha));
  q=pow(1-alpha,2);

  a=p/(q*vk);
  return a;
}

double function_b(double alpha)
{
  double b;

  b=1/(1-alpha);
  return b;
  
}

double function_d(double principal_teta,double t,double vs,double alpha)
{
  double p,q,a;

  p=pow(principal_teta,2*(1-alpha));
  q=(1-alpha)*(1-alpha);

  a=p/(q*vs);
  return a;
}


double function_petitf(double t,double vs,double alpha,int debut,int fin,double fzero)
{
  double p,q,c;
  p=pow(function_R(t,debut,fin,fzero),2*(1-alpha));
  q=(1-alpha)*(1-alpha);

  c=p/(q*vs);
  return c;
}

double function_gplus(double t,double principal_teta,double vs,int debut,int fin,double fzero)
{
  double p,gplus;
  p=log(function_R(t,debut,fin,fzero)/principal_teta);
  
  gplus=(p+0.5*vs)/sqrt(vs);

  return gplus;
  
}

double function_gmoins(double t,double principal_teta,double vs,int debut,int fin,double fzero)
{
  double p,gmoins;
  p=log( function_R(t,debut,fin,fzero)/principal_teta);

  gmoins=(p-0.5*vs)/sqrt(vs);

  return gmoins;
  
}

double function_c(int k,double t,double alpha,double fzero,int n)
{
  double vk,p,q,c;
  vk=calcul_integral(function_lambda_k_carre,t,maturity_structure[k],n);
  p=pow(fzero,2*(1-alpha));
  q=pow((1-alpha),2);

  c=p/(q*vk);
  return c;
}

double function_xplus(double t,int k,double alpha,double H,double fzero,int n)
{
  double vk,p,xplus;
  vk=calcul_integral(function_lambda_k_carre,t,maturity_structure[k],n);
  p=log(fzero/H);
  
  xplus=(p+0.5*vk)/pow(vk,0.5);

  return xplus;
  
}

double function_xmoins(double t, int k,double alpha,double H,double fzero,int n)
{
  double vk,p,xmoins;
  vk=calcul_integral(function_lambda_k_carre,t,maturity_structure[k],n);
  p=log( fzero/H);

  xmoins=(p-0.5*vk)/sqrt(vk);

  return xmoins;
  
}


/* // SIMULATION LOI CHI2 */



  
double distrib_normale_positive(double x)
{
  double N,p,b1,b2,b3,b4,b5,t;

  p=0.2316419;
  b1=0.319381530;
  b2=-0.356563782;
  b3=1.781477937;
  b4=-1.821255978;
  b5=1.330274429;
  t=1/(1+p*x);

  N=1-exp(-x*x/2)*(b1*t+b2*pow(t,2)+b3*pow(t,3)+b4*pow(t,4)+b5*pow(t,5))/pow(2*PI,0.5);
  return N;
  
}


double distrib_normale_negative(double x)
{
  double N,N1;

  N1=distrib_normale_positive(-x);
  N=1-N1;

  return N;
  
}

double distrib_normale(double x)/*marche ok*/
{
  double N;
  
  if(x==0)
	N=0.5;

  if(x>0)
	N=distrib_normale_positive(x);

  if(x<0)
	N=distrib_normale_negative(x);

  return N;
}




double alngam(double xvalue)
{
  long xlgst;
  long xlge;
  double alr2pi,x,x1,x2,y,alngam;
  long r1[9]={-2.66685511495, -2.44387534237*10, -2.19698958928*10, 1.11667541262*10,  3.13060547623,  6.07771387771/10, 1.19400905721*10,  3.14690115749*10, 1.52346874070*10};
  long r2[9]={-7.83359299449*10, -1.42046296688*100,1.37519416416*100,  7.86994924154*10, 4.16438922228,  4.70668766060*10, 3.13399215894*100,  2.63505074721*100,4.33400022514*10};
  long r3[9]={-2.12159572323*100000,  2.30661510616*100000, 2.74647644705*10000, -4.02621119975*10000,-2.29660729780*1000, -1.16328495004*100000,-1.46025937511*100000, -2.42357409629*10000,-5.70691009324*100};
  
  long r4[5]={2.79195317918525/10, 4.917317610505968/10,6.92910599291889/100, 3.350343815022304,6.012459259764103};
  
  alr2pi=9.18938533204673/10;
  xlge=5.10*1000000;
  xlgst=pow(10,30);

  alngam=0;
  x=xvalue;

  if(x<1.5){
	if(x<0.5){
	  alngam=-log(x);
	  y=x+1;
	  if(y==1)
		exit(3);
	  alngam=0;
	  y=x;
	  x=(x-0.5)-0.5;
	}
	alngam=alngam+x*((((r1[4]*y+r1[3])*y+r1[2])*y+r1[1])*y+r1[0])/((((y+r1[8])*y+r1[7])*y+r1[6])*y+r1[5]);
  }

  if(x<4 && x>=1.5){
	y=(x-1)-1;
	alngam=y*((((r2[4]*x+r2[3])*x+r2[2])*x+r2[1])*x+r2[0])/((((x+r2[8])*x+r2[7])*x+r2[6])*x+r2[5]);
  }

  if(x<12 && x>=4){
	alngam=((((r3[5]*x + r3[4])*x + r3[3])*x + r3[2])*x + r3[1]) /((((x + r3[9])*x + r3[8])*x + r3[7])*x + r3[6]);
  }

  if(x>=12){
	y=log(x);
	alngam=x*(y-1)-0.5*y*alr2pi;
	x1=1/x;

	x2=x1*x1;
	alngam=alngam+x1*((r4[2]*x2+r4[1])*x2+r4[0])/((x2+r4[4])*x2+r4[3]);
  }

  return alngam;
  
}

  
double distrib_chi2(double x,double theta/* degrés */,double f/* decentrage */)
{
  double errmax,lam,n,u,v,x2,theta2,t,term,alogam,chi2;
  int itrmax;

  errmax=0.0000000001;
  itrmax=50000;

  chi2=x;
  
  if(theta<=0 || f<0){
	return 0;
	}
  
  if(x<0){
	return 0;
	
  }
  
  lam=f/2;

  n=1;
  u=exp(-lam);
  v=u;
  x2=x/2;
  theta2=theta/2;

  alogam=alngam(theta2+1);
  t=pow(x2,theta2)*exp(-x2)/exp(alogam);

  term=v*t;
  chi2=term;

  
  while(theta+2*n-x<=0){
	u=u*lam/n;
	v=v+u;
	t=t*x/(theta+2*n);
	term=v*t;
	chi2=chi2+term;
	n=n+1;
  }
  
  while(t*x/(theta+2*n-x)>errmax && n<=itrmax){
	u=u*lam/n;
	v=v+u;
	t=t*x/(theta+2*n);
	term=v*t;
	chi2=chi2+term;
	n=n+1;
  }

  return chi2;
}
	  
  
/* FIN simul loi chi2 */


/* CALCUL PRIX CAPLET */




double prix_caplet(double t,int k,double alpha,double H,double fzero,int n)
{
  
  double C_k_t,a,b,c,xplus,xmoins;
  
  if(alpha<=0)
	return 0;

  if(alpha<1 && alpha>0){
    a=function_a(k,t,alpha,H,n);
	b=function_b(alpha);
    c=function_c(k,t,alpha,fzero,n);
	 
	C_k_t=function_delta(k)*price_zero_coupon_bond(t,k+1,fzero)*(fzero*(1-distrib_chi2(a,b+2,c))-H*distrib_chi2(c,b,a));
	 
  }

  if(alpha==1){
	xplus=function_xplus(t,k,alpha,H,fzero,n);
	xmoins=function_xmoins(t,k,alpha,H,fzero,n);
	C_k_t=function_delta(k)*price_zero_coupon_bond(t,k+1,fzero)*(fzero*distrib_normale(xplus)-H*distrib_normale(xmoins));
  }

  if(alpha>1){
    a=function_a(k,t,alpha,H,n);
	b=function_b(alpha);
    c=function_c(k,t,alpha,fzero,n);
	C_k_t=function_delta(k)*price_zero_coupon_bond(t,k+1,fzero)*(fzero*(1-distrib_chi2(c,-b,a))-H*distrib_chi2(a,2-b,c));
  }

  
  return C_k_t;
}

/* CALCUL PRIX swaption */




double prix_swaption(double t,double principal_teta,double vs,double alpha,int debut,int fin,double fzero)
{
  double C_k_t,d,b,f,gplus,gmoins;
  
  if(alpha<=0)
	return 0;

  if(alpha<1 && alpha>0){
    d=function_d(principal_teta,t,vs,alpha);
	b=function_b(alpha);
    f=function_petitf(t,vs,alpha,debut,fin,fzero);
	 
	 C_k_t=function_B_s(t,debut,fin,fzero)*(function_R(t,debut,fin,fzero)*(1-distrib_chi2(d,b+2,f))-principal_teta*distrib_chi2(f,b,d));
	 
  }

  if(alpha==1){
	gplus=function_gplus(t,principal_teta,vs,debut,fin,fzero);
	gmoins=function_gmoins(t,principal_teta,vs,debut,fin,fzero);
	C_k_t=function_B_s(t,debut,fin,fzero)*(function_R(t,debut,fin,fzero)*distrib_normale(gplus)-principal_teta*distrib_normale(gmoins));
  }

  if(alpha>1){
	d=function_d(principal_teta,t,vs,alpha);
	b=function_b(alpha);
    f=function_petitf(t,vs,alpha,debut,fin,fzero);
	 
	 C_k_t=function_B_s(t,debut,fin,fzero)*(function_R(t,debut,fin,fzero)*(1-distrib_chi2(f,-b,d))-principal_teta*distrib_chi2(d,2-b,f));
  }

  
  return C_k_t;
}


double Black(int k,double sigma,double H,double fzero){
  double d1,d2,black;
  d1=(log(fzero/H)+sigma*sigma*maturity_structure[k]/2)/(sigma*pow(maturity_structure[k],0.5));
  d2=(log(fzero/H)-sigma*sigma*maturity_structure[k]/2)/(sigma*pow(maturity_structure[k],0.5));
  black=price_zero_coupon_bond(0,k+1,fzero)*function_delta(k)*(fzero*distrib_normale(d1)-H*distrib_normale(d2));
  
  
  return black;
}

double Black_prime(int k,double sigma,double H,double fzero)
{
  double black_prime,d1,n_d1;
  d1=(log(fzero/H)+sigma*sigma*maturity_structure[k]/2)/(sigma*pow(maturity_structure[k],0.5));
  /*  fprintf(stdout,"%f d1 \n",d1); */
  
  n_d1=exp(-d1*d1/2)/pow(2*PI,0.5);
  black_prime=fzero*n_d1*pow(maturity_structure[k],0.5)*function_delta(k)*price_zero_coupon_bond(0,k+1,fzero);
  
  return black_prime;
}


double implied_vol(double sigma_n,int k,double prix_caplet,double H,double fzero)
{
  while(fabs((prix_caplet-Black(k,sigma_n,H,fzero)))>0.00000000000001){ 

	sigma_n=sigma_n+(prix_caplet-Black(k,sigma_n,H,fzero))/Black_prime(k,sigma_n,H,fzero);
  }

  
  return sigma_n;
  
  
}

void function_F(int e,double t[],double final[],double fzero,double delta_i,double alpha,int type_model,int type_scheme)
{
  
  double mu=0;
  double simul=0;
  double lambda1=0;
  double lambda2=0;
  
  double tableau_fj_ti[150],ancientableau[181];
  double produit=0;
  
  int q,j,k,i,z,o,p,c,a;
  
  for(q=0;q<e+1;q++) tableau_fj_ti[q]=fzero;

  final[0]=fzero;
  for(i=0;i<e*function_delta(e)/delta_i;i++){
    
    simul=simul_normale();
    for(k=function_n(t[i]);k<e+1;k++){
      mu=0;
      ancientableau[k]=tableau_fj_ti[k];
      for(j=function_n(t[i]);j<k+1;j++){
	lambda1=function_lambda_j(j,i);
	mu+=lambda1*(function_delta(j)*function_fi(ancientableau[j],type_model,alpha))/(1+function_delta(j)*ancientableau[j]);
      }
      lambda2=function_lambda_j(k,i);
      /*  Schema log-Eulerien */
      if(type_scheme==0){
	produit=lambda2*((mu-0.5*function_fi(ancientableau[k],type_model,alpha)*lambda2/ancientableau[k])*delta_i+simul*sqrt(delta_i));
	tableau_fj_ti[k]=ancientableau[k]*exp(function_fi(ancientableau[k],type_model,alpha)*produit/ancientableau[k]);
      }
      /*  Schema Eulerien */
      if(type_scheme==1){
	produit=lambda2*(mu*delta_i+simul*sqrt(delta_i));
	tableau_fj_ti[k]=ancientableau[k]+function_fi(ancientableau[k],type_model,alpha)*produit;
      }
      if(k==(i+1)*delta_i/function_delta(e)){final[k]=tableau_fj_ti[k];}
      
    }
  }
  
}


void function_F_swap(int debut,int fin,double t[],double final[],double fzero,double delta_i,double alpha,int type_model,int type_scheme,double finalcor[])
{
  double mu=0;
  double simul=0;
  double lambda1=0;
  double lambda2=0;
  
  double tableau_fj_ti[150],ancientableau[181];
  double produit=0;
  
  int q,j,k,i,z,o,p,c,a;
  
  for(q=0;q<fin+1;q++) tableau_fj_ti[q]=fzero;

  final[0]=fzero;
  
  for(i=0;i<debut*function_delta(fin)/delta_i;i++){
    
    simul=simul_normale();
    for(k=function_n(t[i]);k<fin+1;k++){
      mu=0;
      ancientableau[k]=tableau_fj_ti[k];
      for(j=function_n(t[i]);j<k+1;j++){
	lambda1=function_lambda_j(j,i);
	mu+=lambda1*(function_delta(j)*function_fi(ancientableau[j],type_model,alpha))/(1+function_delta(j)*ancientableau[j]);
      }
      lambda2=function_lambda_j(k,i);
      /*  Schema log-Eulerien */
      if(type_scheme==0){
	produit=lambda2*((mu-0.5*function_fi(ancientableau[k],type_model,alpha)*lambda2/ancientableau[k])*delta_i+simul*sqrt(delta_i));
	
	tableau_fj_ti[k]=ancientableau[k]*exp(function_fi(ancientableau[k],type_model,alpha)*produit/ancientableau[k]);
      }
      

      /*  Schema Eulerien */
      if(type_scheme==1){
	produit=lambda2*(mu*delta_i+simul*sqrt(delta_i));
	tableau_fj_ti[k]=ancientableau[k]+function_fi(ancientableau[k],type_model,alpha)*produit;
      }
      /*  to get the Fk in T[debut] */
      if((i+1)*delta_i==debut*function_delta(fin)) final[k]=tableau_fj_ti[k];
      /*   to get the Fk in Tk */
      if(k==(i+1)*delta_i/function_delta(fin)){finalcor[k]=tableau_fj_ti[k];}
      
    }
  }
  
}


double function_B_en_maturite(double o,double F[])
{

  int i;
  double b=1;
  

  for(i=0;i<function_n(o)-1;i++){
	b *= (1+function_delta(i)*F[i]);
  }
  
  return b;
}

