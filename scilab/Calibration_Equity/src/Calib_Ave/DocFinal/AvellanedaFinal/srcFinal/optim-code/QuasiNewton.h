//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//

/* 
** FILE
** *****************************************************************
** 
** QuasiNewton.h
** 
** 
** PURPOSE
** 
** Implements a Quasi-Newton unconstrained optimizer. Uses a BFGS update
** of the Cholesky factors of the approximation of the hessian. An initial
** scaling is made in the first iteration.
** 
** 
*****************************************************************
*/

#ifndef QUASINEWTON
#define QUASINEWTON

double macheps();//Calculates the machine epsilon

void QuasiNewton(int n, 
		 double* x0, 
		 double (*costFunction)(double*), 
		 void (*gradCostFunction)(double*,double*),
		 double* sol,
		 double gradtol,
		 double steptol,
		 int maxCounter,
		 int verbosity,
		 int saveSuccessiveXinFile);

void testQuasiNewton();

#endif
