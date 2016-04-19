#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lmm_mathtools.h"

// pour le changement de variable qui ramene [a,b] a [0,1]
#define FUNCG(x) ((*funcg)(a + (b-a)*(x)))
//


static int ngauss=-1;
static double *xi,*wi;


/*One-Dimensional Normal Law. Cumulative distribution function. */
/*Abramowitz, Milton et Stegun, Handbook of MathematicalFunctions, 1968, Dover Publication, New York, page 932 (26.2.18).Precision 10-7*/
double N(double x)
{
  const double p= 0.2316419;
  const double b1= 0.319381530;
  const double b2= -0.356563782;
  const double b3= 1.781477937;
  const double b4= -1.821255978;
  const double b5= 1.330274429;
  const double one_over_twopi= 0.39894228;
  
  double t;
  
  if(x >= 0.0)
    {
      t = 1.0 / ( 1.0 + p * x );
      return (1.0 - one_over_twopi * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    } 
  else 
    {/* x < 0 */
      t = 1.0 / ( 1.0 - p * x );
      return ( one_over_twopi * exp( -x * x / 2.0 ) * t * ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
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

/* Given the lower und upper limits of integration x1 and x2, 
   and given n, this routine returns arrays x[1..n] and 
   w[1..n], containing the abscissas and weights of the 
   Gauss-Legendre n-point quadrature formula.
   The integral of f will be computed in the following way:

   double integral = 0.0;
   for(i = 1; i < (n+1); i++)
     integral += w[i] * f(x[i]);
*/

static void gauleg(double x1, double x2, double x[], double w[], int n)
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
/*-------------------------------------------------------------------*/
/* Useful Routines to solve Least Squares Regression problems        */
/*-------------------------------------------------------------------*/ 


//Cholesky square root of a dim*dim symmetric matrix M; 
//Input:upper triangle of M; Output:lower triangle of sqrt(M)
//(square root is lower triangular matrix)

int Cholesky(double *M, int dim)
{
  int i,j,k, error=0;
  double somme;
 
  i=0;
  while (!(error)&&(i<dim)){
    for (j=i;j<dim;j++){
      for (somme=M[i*dim+j],k=i-1;k>=0;k--)
        somme-=M[i*dim+k]*M[j*dim+k];
      if (i==j) {
        if (somme<=0) error=1;
        else M[j*dim+i]=sqrt(somme);
      }
          else M[j*dim+i]=somme/M[i*dim+i];
    }
    i++;
  }
  //setting upper triangle of M to zero
  for (i=0;i<dim;i++){
    for (j=i+1;j<dim;j++) M[i*dim+j]=0.0;
  }
  return error;
}

//  Resolution of the equation 
//  M*Res=AuxR,  
//  M=lower triangular dim*dim matrix, AuxR=1*dim array
//  Res=output.

void Resolution(double *AuxR, double *Res, double *M, int dim)
{
  int i,k;
  double somme;

  for (i=0;i<dim;i++){
	for (somme=AuxR[i],k=i-1;k>=0;k--)
	  somme-=M[i*dim+k]*Res[k];
	Res[i]=somme/M[i*dim+i];
  }
  for (i=dim-1;i>=0;i--){
	for (somme=Res[i],k=i+1;k<dim;k++)
	  somme-=M[k*dim+i]*Res[k];
	Res[i]=somme/M[i*dim+i];
  }
}


