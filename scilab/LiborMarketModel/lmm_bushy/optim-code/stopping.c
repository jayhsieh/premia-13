//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//


/* 
** FILE
** *****************************************************************
** 
** stopping.c
** 
** PURPOSE
** 
** stopping criteria for iterative optimizers.
**  
*****************************************************************
*/

#include <stdio.h>
#include <math.h>
#include "malloc.h"
#include "stopping.h"




/* 
** FUNCTION
** *****************************************************************
** 
** max(double* t, int n)
** 
** DESCRIPTION
** 
** returns the maximum of a vector of numbers. The size of the vector is n.
** 
*****************************************************************
*/

double max(double* t, int n)
{
  double m=t[0];
  double p;
  int i;
  for (i=1;i<n;i++)
    {
      p=t[i];
      if (p>m) {m=p;}
    }
  return m;
}


/* 
** FUNCTION
** *****************************************************************
** 
** relgrad
** 
** DESCRIPTION
** 
** computes the current relative gradient: max ( | grad[i]*x[i] / f(x) | )
**                                          i
** To avoid problems when x[i] or f(x) are near zero, one replaces
** x[i] by max(|x[i]|,typx[i]) and f(x) by max(|f(x)|,typf)
** where typf ( resp. typx[i]) are typical magnitude of f (resp. x[i])
**
*****************************************************************
*/

double relgrad(int n, 
	       double* x0, 
	       double* grad0, 
	       double fx0, 
	       double* typx,
	       double typf)
{
  int i;
  double xi, typi, num, den, res;
  double* t = (double *) malloc( n * sizeof(double));
  for(i=0;i<n;i++)
    {
      xi  = fabs(x0[i]);
      typi= fabs(typx[i]);
      num = xi>typi ? xi:typi;
      num = num*fabs(grad0[i]);
      den = fabs(fx0)>fabs(typf) ? fabs(fx0):fabs(typf);
      t[i]= num/den;
    }
  res=max(t,n);
  free(t);
  return res;
}


/* 
** FUNCTION
** *****************************************************************
** 
** relx
** 
** DESCRIPTION
** 
** computes the current relative change in x : max ( |x'[i]-x[i]| / |x'[i]| )
**                                              i
** To avoid problems when x'[i] is near zero, one replaces x'[i] by
** max(|x'[i]|,typx[i]) where typx[i] is typical magnitude of x[i]. 
** 
*****************************************************************
*/

double relx(int n,double* x1,double* x0,double* typx)
{
  int i;
  double xi, typi, num, den;
  double* t = (double *) malloc( n * sizeof(double));
  for(i=0;i<n;i++)
    {
      xi  = fabs(x1[i]);
      typi= fabs(typx[i]);
      num = fabs(x1[i]-x0[i]);
      den = xi>typi ? xi:typi;
      t[i]= num/den;
    }
  return max(t,n);
}
