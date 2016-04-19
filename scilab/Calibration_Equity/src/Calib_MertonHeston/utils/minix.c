#include "minix.h"
//
static int nbinit=0;
static int dimx0, *indicx0;
static double *x0,*x1,*grad1;
//
double costFunction1(int dimx,double *x)
{
  double fx;
  //
  minix2x(x1,x);
  fx =  costFunction(dimx0,x1);
  //
  return fx;
  //
}
//
void gradFunction1(int dimx,double *x,double *grad)
{
  //
  minix2x(x1,x);
  gradFunction(dimx0,x1,grad1);
  x2minix(grad1,grad);
  //
}
//
void x2minix(double *x,double *minix)
{
  int i,j=0;
  
  for(i=0;i<dimx0;i++)
	{
	  if(indicx0[i] == 1)
		j++;
	  else 
		minix[i-j] = x[i];
	}
}
//
void minix2x(double *x,double *minix)
{
  int i,j=0;
  for(i=0;i<dimx0;i++)
	{
	  if(indicx0[i] == 1)
		{
		  j++;
		  x[i] = x0[i];
		} else {
		  x[i] = minix[i-j];
		}
	}
}
//
int initminix(int dimx,int *nbd,double *xmin,double *xmax,int *mininbd,double *minixmin,double *minixmax)
{
  int i,j=0;
  int dimminix;
  if(nbinit==0)
	{
	  indicx0 = (int *) malloc(dimx*sizeof(int));
	  x0 = (double *) malloc(dimx*sizeof(double));
	  x1 = (double *) malloc(dimx*sizeof(double));
	  grad1  = (double *) malloc(dimx*sizeof(double));
	}
  nbinit++;
  
  dimx0    = dimx;
  dimminix = dimx;
  //
  for(i=0;i<dimx;i++)   indicx0[i] = 0;
  
  for(i=0;i<dimx;i++)
	{ if(nbd[i]==2 && xmin[i] == xmax[i])
	  {
		indicx0[i] = 1;
		x0[i] = xmin[i];
		j++;
	  } else {
		mininbd[i-j]  = nbd[i];
		minixmin[i-j] = xmin[i];
		minixmax[i-j] = xmax[i];
	  }
	
	}
  //
  dimminix = dimx - j;
  //
  //
  return dimminix;
  //
}
//
void FreeMinix()
{
  free(indicx0);
  free(x0);
  free(x1);
  free(grad1);
  nbinit = 0;
  
}

  
