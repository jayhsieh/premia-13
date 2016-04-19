//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//


/* 
** FILE
** *****************************************************************
** 
** BFGSupdate.c 
** 
** 
** PURPOSE
** 
** Implementation of the BFGS update of the Cholesky factorization
** of the matrix which approximates the hessian. (The modification is 
** considered as two succesive rank-one modifications)
** 
** 
*****************************************************************
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "BFGSupdate.h"
#include "cholesky.h"


/* 
** FUNCTION
** *****************************************************************
** 
** BFGSupdate
** 
** 
** DESCRIPTION
** 
** updates B - Bs(Bs)'/(s'Bs) + yy'/y's in a new cholesky factorization
** 
** where s = x1 - x0
**       y = g1 - g0
** 
** 
** 
** 
*****************************************************************
*/

void BFGSupdate(struct chol *C,//the cholesky factorization to update 
		double* x1,//the current point "x(n+1)" 
		double* x0,//the current point "x(n)"
		double* g1,//the current gradient "g(n+1)" at x(n+1)
		double* g0,//the current gradient "g(n)" at x(n)
		int n,//the size of the problem 
		double epsilon,//the machine precision 
		int *endBFGS,//test wether the denominators don't equal zero 
		int verbosity)//level of information choosen by the user
{
  int i;
  
  double* s;//s=x1-x0
  double* y;//y=g1-g0
  double* t;//"t=Bs"
  double* v;//v is s.t.  Lv=t (B=LDL')
  double* Bs;//the product "new BS"*s suppose to equal y (verification
             //of the secant equation)
  
  
  double alpha=0;//alpha will be 1/y's
  double beta=0;//beta will be 1/(s'Bs)
  
  s=(double *) malloc( n * sizeof(double));
  y=(double *) malloc( n * sizeof(double));
  t=(double *) malloc( n * sizeof(double));
  v=(double *) malloc( n * sizeof(double));
  Bs=(double *) malloc(n*sizeof(double));
  
  if (verbosity>=3) {printf("BFGS update:");}

  for (i=1;i<=n;i++)
    {
      s[i-1]= x1[i-1]-x0[i-1];
      y[i-1]= g1[i-1]-g0[i-1];
      alpha+= y[i-1]*s[i-1];//saclar product y's
    }
  
  prodMatVect(C,s,t);//"t=Bs"
  
  for (i=1;i<=n;i++) { beta+= s[i-1]*t[i-1]; }//saclar product s'Bs
  
  if (beta==0) { *endBFGS=1;printf("\tbeta infinite in BFGS update"); } 
  else 
    {
      beta = -1/beta;//beta=1/(s'Bs)
      if (verbosity>=3) {printf("\t-1/sBs=%e", beta);}
      if (alpha==0) { *endBFGS=1;printf("\talpha infinite in BFGS update"); }
      else {
	alpha= 1/alpha;
	if (verbosity>=3) {printf("  1/ys=%e\n", alpha);}

	
	/************ The two successive rank-one modifications ********/
	
	positiveRankOneModification(C, alpha, y);
	Lsolve(C,v,t);//v is s.t.  Lv=t
	negativeRankOneModification(C, beta, t, v, epsilon); 
	
	/******* End of the two successive rank-one modifications ******/
	
	
	if (verbosity>=3)
	  {
	    double N2=0;
	    printf("\t\tVerification of secant equation: ");//We verify (new B)s=y
	    prodMatVect(C,s,Bs);
	    for (i=0;i<n;i++) { N2+= pow((g1[i]-g0[i])-Bs[i],2);}
	    printf("N2( B's-y )=%e\n ");
	  }
      }
    }
  
  //memory disallocation
  free(s);
  free(y);
  free(t);
  free(v);
  free(Bs);
}


