#include "my_integral.h"
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#define FUNC(x) ((*func)(x))
#define FUNK(x) (2.0*(x)*(*func2)(aa+(x)*(x)))
#define FUNKY(x) ((*func1)(1.0/(x))/((x)*(x)))
// pour le changement de variable qui ramene [a,b] a [0,1]
#define FUNCG(x) ((*funcg)(a + (b-a)*(x)))
//
#define FUNCG_VECT(x,n,fx) ((*funcg_vect)( ((a + (b-a)*(x))) , n , fx) )                                                                                                                            
#define NR_END 1
#define FREE_ARG char*


static int ngauss=-1;
static double *xi,*wi;



double *vector(long nl,long nh){ 
  double *v;

  v=(double *)malloc((size_t)((nh-nl+1+NR_END)*sizeof(double)));
  /*if(!v) callerror("allocation failure in vector()");*/
  return v-nl+NR_END;
}


void free_vector(double *v,long nl,long nh){
  free((FREE_ARG)(v+nl-NR_END));
}


double midpnt(double (*func)(double), double a, double b, int n){
  double x,tnm,sum,del,ddel;
  static double s;
  int it,j;

  if(n==1){
    s=(b-a)*FUNC(0.5*(a+b));
	return s;
  } else {
	for(it=1,j=1;j<n-1;j++) it*=3;
	tnm=it;
	del=(b-a)/(3.0*tnm); 
	ddel=del+del;        
	x=a+0.5*del;        
	sum=0.0;
	for(j=1;j<=it;j++) {
	  sum+=FUNC(x);
	  x+=ddel;
	  sum+=FUNC(x);
	  x+=del;
	  
	}
	s=( midpnt(func,a,b,n-1) + (b-a)*sum/tnm )/3.0;
	return s;                
  }
}

double midpntbis(double (*func)(double), double a, double b, int n){
  double x,tnm,sum,del,ddel;
  static double s;
  int it,j;

  if(n==1){
      s=(b-a)*FUNC(0.5*(a+b));
	return s;
	
  } else {
	for(it=1,j=1;j<n-1;j++) it*=3;
	tnm=it;
	del=(b-a)/(3.0*tnm); 
	ddel=del+del;        
	x=a+0.5*del;         
	sum=0.0;
	for(j=1;j<=it;j++){
	  sum+=FUNC(x);
	  x+=ddel;
	  sum+=FUNC(x);
	  x+=del;
	}
	s=( midpntbis(func,a,b,n-1) + (b-a)*sum/tnm )/3.0;
	return s;                
  }
}

double midsql(double (*func2)(double), double aa, double bb, int n){
  double x,tnm,sum,del,ddel,b,a;
  static double s;
  int it,j;

  b=sqrt(bb-aa);
  a=0.0;
  if(n==1) {
         s=(b-a)*FUNK(0.5*(a+b));
	return s;
  } else {
	for(it=1,j=1;j<n-1;j++) it*=3;
	tnm=it;
	del=(b-a)/(3.0*tnm); 
	ddel=del+del;       
	x=a+0.5*del;         
	sum=0.0;
	for(j=1;j<=it;j++){
	  sum+=FUNK(x);
	  x+=ddel;
	  sum+=FUNK(x);
	  x+=del;
	}
	s=( midsql(func2,aa,bb,n-1) + (b-a)*sum/tnm )/3.0;
	return s;               
  }
}


double midsqlbis(double (*func2)(double), double aa, double bb, int n){
  double x,tnm,sum,del,ddel,b,a;
  static double s;
  int it,j;

  b=sqrt(bb-aa);
  a=0.0;
  if(n==1) {
    s=(b-a)*FUNK(0.5*(a+b));
	return s;
  } else {
	for(it=1,j=1;j<n-1;j++) it*=3;
	tnm=it;
	del=(b-a)/(3.0*tnm); 
	ddel=del+del;        
	x=a+0.5*del;         
	sum=0.0;
	for(j=1;j<=it;j++) {
	  sum+=FUNK(x);
	  x+=ddel;
	  sum+=FUNK(x);
	  x+=del;
	}
	s=( midsqlbis(func2,aa,bb,n-1) + (b-a)*sum/tnm )/3.0;
	return s;               
  }
}

double midinf(double (*func1)(double), double aa, double bb, int n){
  double x,tnm,sum,del,ddel,b,a;
  static double s;
  int it,j;

  b=1.0/aa; 
  a=1.0/bb; 
  if(n==1){
    s=(b-a)*FUNKY(0.5*(a+b));
    return s;
  } else {
	for(it=1,j=1;j<n-1;j++) it*=3;
	tnm=it;
	del=(b-a)/(3.0*tnm); 
	ddel=del+del;        
	x=a+0.5*del;         
	sum=0.0;
	for(j=1;j<=it;j++){
	  sum+=FUNKY(x);
	  x+=ddel;
	  sum+=FUNKY(x);
	  x+=del;
	}
	s=( midinf(func1,aa,bb,n-1) + (b-a)*sum/tnm )/3.0;
	return s;                
  }
}


/*--------------------------------------------------------------------*/
double integrale_gauss(double (*funcg)(double), double a, double b){
  int i;
  double sum=0.;

  if(ngauss<0) {
	 printf("Erreur : vous devez initialiser les points de les poids de Gauss.\n");
	 exit(-1);
  }
  
  for(i=1;i<=ngauss;i++)
	{
	  sum+=wi[i] * FUNCG(xi[i]);
	}
// pour le changement de variable qui ramene [a,b] a [0,1]
  sum*=(b-a);
  return sum;
}

/*--------------------------------------------------------------------*/
void integrale_gauss_vect(void (*funcg_vect)(double,int,double *), double a, double b, int dimx, double *sum)
{
  int i,n;
  double *fx;
  double x;
  // initialisation
  for(n=0;n<dimx;n++)  sum[n]=0.;
  //
  fx=(double*)malloc(sizeof(double)*dimx);
  //
  if(ngauss<0) {
	printf("Erreur : vous devez initialiser les points de les poids de Gauss.\n");
	exit(-1);
  }
  for(i=1;i<=ngauss;i++)
	{
	  x = xi[i];
	  FUNCG_VECT(x,dimx,fx);
	  for(n=0;n<dimx;n++)    
		sum[n] = sum[n] + wi[i] * fx[n];
	}
  // pour le changement de variable qui ramene [a,b] a [0,1]
  for(n=0;n<dimx;n++)
	sum[n]=sum[n]*(b-a);
  //
  free(fx);
  //
}
/*--------------------------------------------------------------------*/

/* Given the lower und upper limits of integration x1 and x2, 
   and given n, this routine returns arrays x[1..n] and 
   w[1..n], containing the abscissas and weights of the 
   Gauss-Legendre n-point quadrature formula.
   The integral of f will be computed in the following way:

   double integral = 0.0;
   for(i = 1; i < (n+1); i++)
     integral += w[i] * f(x[i]);
*/

void gauleg(double x1, double x2, double x[], double w[], int n)
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  double EPS1=3.0e-11;
  
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) {
	z=cos(3.141592654*(i-0.25)/(n+0.5));
	do {
	  p1=1.0;
	  p2=0.0;
	  for (j=1;j<=n;j++) {
		p3=p2;
		p2=p1;
		p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	  }
	  pp=n*(z*p1-p2)/(z*z-1.0);
	  z1=z;
	  z=z1-p1/pp;
	} while (fabs(z-z1) > EPS1);
	x[i]=xm-xl*z;
	x[n+1-i]=xm+xl*z;
	w[i]=2.0*xl/((1.0-z*z)*pp*pp);
	w[n+1-i]=w[i];
  }
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/

void init_gauss(int nbpts)
{
  ngauss = nbpts;
  xi = (double*)malloc(sizeof(double)*(ngauss+1));
  wi = (double*)malloc(sizeof(double)*(ngauss+1));
  gauleg(0.,1.,xi,wi,ngauss);
}
/*--------------------------------------------------------------------*/
void free_gauss()
{
  //
  free(xi);
  free(wi);
  //  
}
