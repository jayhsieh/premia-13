//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//


/* 
** FILE
** *****************************************************************
** 
** cholesky.c
** 
** 
** PURPOSE
** 
** Some methods for the cholesky factorization LDL' of a matrix.
** 
** 
** 
** 
*****************************************************************
*/

#include <stdio.h>
#include <math.h>

#include "cholesky.h"


/* 
** FUNCTION
** *****************************************************************
** 
** init
**  
** DESCRIPTION
** 
** initialize with the identity matrix.
** 
*****************************************************************
*/

void init(struct chol *C)
{
  int i;
  int n=C->size;

  for(i=0;i<n;i++)
    {
      C->d[i]=1.;
    }
  for(i=0;i<(n*(n-1))/2;i++)
    {
      C->l[i]=0.;
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** elementLij
**  
** DESCRIPTION
** 
** Returns the element of the vector l correspoding to the element
** (i,j) of the matrix L. 
** 
*****************************************************************
*/

double elementLij(int i,int j,struct chol *C)    // i>j !!!
{
  int k= i -2 +( (j-1)*( 2*C->size -j -2 ) )/2;
  return C->l[k];
}


/* 
** FUNCTION
** *****************************************************************
** 
** elementDi
** 
** DESCRIPTION
** 
** Returns the ith element of the vector d.
** 
*****************************************************************
*/

double elementDi(int i,struct chol *C)
{
  return C->d[i-1];
}


/* 
** FUNCTION
** *****************************************************************
** 
** elementMij
** 
** DESCRIPTION
** 
** returns the element (i,j) of the product LDL' (for i>=j, remember Mij=Mji)
**  
*****************************************************************
*/

double elementMij(int i,int j,struct chol *C)// i>=j, remember Mij=Mji
{
  double m=0;
  int k;
  for (k=1;k<j;k++)
    {
      m=m + elementLij(i,k,C)*elementDi(k,C)*elementLij(j,k,C);
    }
  if (i!=j) {m=m + elementDi(j,C)*elementLij(i,j,C);} else {m=m + elementDi(j,C);}
  return m;
}


/* 
** FUNCTION
** *****************************************************************
** 
** modify_l
** 
** DESCRIPTION
** 
** Changes the element (i,j) of the matrix L into the value ll
** 
*****************************************************************
*/

void modify_l(int i,int j,double ll,struct chol *C)
{
  C->l[i -2 +( (j-1)*( 2*C->size -j -2 ) )/2]=ll;
}

/* 
** FUNCTION
** *****************************************************************
** 
** modify_d
** 
** DESCRIPTION
** 
** Changes the ith element of the matrix D into the value dd
** 
*****************************************************************
*/

void modify_d(int i,double dd,struct chol *C)
{
  C->d[i-1]=dd;
}

/* 
** FUNCTION
** *****************************************************************
** 
** scalarMult
** 
** DESCRIPTION
** 
** Multiplies the matrix by the scalar a
** 
*****************************************************************
*/

void scalarMult(struct chol *C, double a)
{
  int i;
  int n=C->size;

  for(i=0;i<n;i++)
    {
      C->d[i]*=a;
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** prodMatVect(struct chol *C,double* v, double* res)
** 
** DESCRIPTION
** 
** Returns in res the product Cv
** 
*****************************************************************
*/

void prodMatVect(struct chol *C,double* v, double* res)
{
  int i,j;
  int n=C->size;
  for (i=1;i<=n;i++)
    {
      res[i-1]=elementMij(i,i,C)*v[i-1];
      for (j=1;j<i;j++)
	{
	  double m=elementMij(i,j,C);
	  res[i-1]+= m*v[j-1];
	  res[j-1]+= m*v[i-1];
	}
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** Lsolve(struct chol *C, double* x, double* b)
** 
** DESCRIPTION
** 
** Changes x into the solution of the linear system Lx=b
** 
*****************************************************************
*/

void Lsolve(struct chol *C, double* x, double* b)
{
  int n=C->size;
  int i,j;
  double sum;

  x[0]=b[0];
  for (i=2;i<n+1;i++)
    {
      sum=0;
      for (j=1;j<i;j++)
	{
	  sum+= elementLij(i,j,C) * x[j-1];
	}
      x[i-1]= b[i-1] - sum;
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** LTsolve(struct chol *C, double* x, double* b)
** 
** DESCRIPTION
** 
** Changes x into the solution of the linear system L'x=b
** 
*****************************************************************
*/

void LTsolve(struct chol *C, double* x, double* b)
{
  int n=C->size;
  int i,j;
  double sum;

  x[n-1]=b[n-1];
  for (i=n-1;i>0;i--)
    {
      sum=0;
      for (j=i+1;j<n+1;j++)
	{
	  sum+= elementLij(j,i,C) * x[j-1];
	}
      x[i-1]= b[i-1] - sum;
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** cholSolve(struct chol *C, double* x, double* b)
** 
** DESCRIPTION
** 
** Changes x into the solution of the linear system Ax=b
** 
*****************************************************************
*/

void cholSolve(struct chol *C, double* x, double* b)// solve Ax=b
{
  int i;

  Lsolve(C, x, b);
  for (i=0; i<C->size; i++)
    {
      x[i]=x[i]/ (C->d[i]);
    }
  LTsolve(C, x, x);
}


/* 
** FUNCTION
** *****************************************************************
** 
** cholSolveOpp(struct chol *C, double* x, double* b)
** 
** DESCRIPTION
** 
** Changes x into the solution of the linear system Ax=-b
** 
*****************************************************************
*/

void cholSolveOpp(struct chol *C, double* x, double* b)// solve Ax=-b
{
  int i;

  Lsolve(C, x, b);
  for (i=0; i<C->size; i++)
    {
      x[i]=-x[i]/ (C->d[i]);
    }
  LTsolve(C, x, x);
}


/* 
** FUNCTION
** *****************************************************************
** 
** positiveRankOneModification(struct chol *C, double s, double* z)
** 
** DESCRIPTION
** 
** C is changed into the cholesky factorization of C + s*z'z (s>0)
** 
*****************************************************************
*/

void positiveRankOneModification(struct chol *C, double s, double* z)
{
  double t=1/s;
  double tt, v, beta;
  int i,j;
  
  for (j=1;j<C->size+1;j++)
    {
      v=z[j-1];
      tt=t+v*v/C->d[j-1];
      beta=(v/C->d[j-1])/tt;
      modify_d(j,C->d[j-1]*tt/t,C);
      for (i=j+1;i<C->size+1;i++)
	{
	  z[i-1]+= -v*elementLij(i,j,C);
	  modify_l(i,j,elementLij(i,j,C)+ beta*z[i-1],C);
	}
      t=tt;
      
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** negativeRankOneModification
** 
** DESCRIPTION
** 
** C is changed into the cholesky factorization of C - s*z'z (s>0)
** The algorithm is different from the previous one:
**    - a vector v s.t. Lv=z is needed.
**    - the machine precision "epsilon" is needed.
** 
*****************************************************************
*/

void negativeRankOneModification(struct chol *C, double s, double* z, double* v, double epsilon)  /* v is such that Lv=z, s<0 */
{
  double tt=1/s;
  double t, beta, vv, zz;
  int i,j;
  
  for (j=1;j<C->size+1;j++)
    {
      tt=tt+v[j-1]*v[j-1]/C->d[j-1];
    }
  if (s*tt<=epsilon) {tt=epsilon/s;}
  for (j=C->size;j>=1;j--)
    {
      vv=v[j-1];
      t=tt-vv*vv/C->d[j-1];
      beta=(vv/C->d[j-1])/tt;
      modify_d(j,C->d[j-1]*tt/t,C);
      z[j-1]=vv;
      for (i=j+1;i<C->size+1;i++)
	{
	  zz=z[i-1];
	  z[i-1]+= vv*elementLij(i,j,C);
	  modify_l(i,j,elementLij(i,j,C)+ beta*zz,C);
	}
      tt=t;
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** printChol, printCholM
** 
** DESCRIPTION
** 
** printChol prints the matrix L and D in the terminal.
** printCholM  prints the matrix M=LDL' in the terminal.
**  
*****************************************************************
*/

void printChol(struct chol *C)
{
  int i,j;
  
  printf("matrix L: \n");
  for (i=1;i<C->size+1;i++)
    {
      for (j=1;j<i;j++)
	{
	  printf("%f",C->l[i -2 +( (j-1)*( 2*C->size -j -2 ) )/2] );
	}
      printf("%f\n",1. );
    }
  printf("\n");
  printf("matrix D: \n");
  for (j=1;j<i;j++) {printf("%f",C->d[j-1]);}
  printf("\n");
}

void printCholM(struct chol *C)
{
  int i,j;
  
  printf("matrix M: \n");
  for (i=1;i<=C->size;i++)
    {
      for (j=1;j<=i;j++) {printf(" %f",elementMij(i,j,C));}
      for (j=i+1;j<=C->size;j++) {printf(" %f",elementMij(j,i,C));}
      printf("\n");
    }  
}
