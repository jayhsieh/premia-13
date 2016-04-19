//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//

/* 
** FILE
** *****************************************************************
** 
** stopping.h
** 
** PURPOSE
** 
** stopping criteria for iterative optimizers.
**  
*****************************************************************
*/

#ifndef STOPPING
#define STOPPING


double relgrad(int n, 
	       double* x0, 
	       double* grad0, 
	       double fx0, 
	       double* typx,
	       double typf);
// Computes the current relative gradient: max ( | grad[i]*x[i] / f(x) | )
//                                          i
// To avoid problems when x[i] or f(x) are near zero, one replaces
// x[i] by max(|x[i]|,typx[i]) and f(x) by max(|f(x)|,typf)
// where typf ( resp. typx[i]) are typical magnitude of f (resp. x[i])


double relx(int n,double* x1,double* x0,double* typx);
// Computes the current relative change in x : max ( |x'[i]-x[i]| / |x'[i]| )
//                                              i
// To avoid problems when x'[i] is near zero, one replaces x'[i] by
// max(|x'[i]|,typx[i]) where typx[i] is typical magnitude of x[i]. 

#endif
