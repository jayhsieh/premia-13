//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//


/* 
** FILE
** *****************************************************************
** 
** BFGSupdate.h 
** 
** 
** PURPOSE
** 
** Implementation of the BFGS update of the Cholesky factorization
** of the matrix which approximates the hessian. (The modification is 
** considered as two succesive rank-one modifications):
** updates B - Bs(Bs)'/(s'Bs) + yy'/y's in a new cholesky factorization,
** where s = x1 - x0
**       y = g1 - g0. 
** 
*****************************************************************
*/

#ifndef BFGSUPDATE
#define BFGSUPDATE

#include "cholesky.h"

void BFGSupdate(struct chol *C,//the cholesky factorization to update 
		double* x1,//the current point "x(n+1)" 
		double* x0,//the current point "x(n)"
		double* g1,//the current gradient "g(n+1)" at x(n+1)
		double* g0,//the current gradient "g(n)" at x(n)
		int n,//the size of the problem 
		double epsilon,//the machine precision 
		int *endBFGS,//test wether the denominators don't equal zero 
		int verbosity);//level of information choosen by the user


#endif
