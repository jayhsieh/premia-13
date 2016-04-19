//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//

/* 
** FILE
** *****************************************************************
** 
** linesearch.h
** 
** 
** PURPOSE
** 
** implementation of the Wolfe's line-search with cubic interpolation.
** 
** 
** 
** 
*****************************************************************
*/

#ifndef LINESEARCH
#define LINESEARCH 


double meriteFunction( int n,//dimension of the problem
		       double (*costFunction)(double*),//cost function F
		       double* point,//previous point of the optimisation algorithm
		       double* direction,//direction of the line search 
		       double t);//step t                               
//computes F(p+td)


double difMeriteFunction( int n,//dimension of the problem
			  void (*gradCostFunction)(double*, double*),//gradient of F
			  double* point,//previous point of the optimisation algorithm
			  double* direction, //direction of the line search 
			  double t);//step t
//computes (d/dt)(F(p+td))


double cubicMin(double t1, double t2, double q1, double q2, double dq1, double dq2);
//Given the doubles t1, t2, q1, q2, dq1, dq2, computes the argmin of a cubic function f such that: f(t1)=q1, f'(t1)=dq1, f(t2)=q2 and f'(t2)=dq2.


void testcubicMin();


void WolfeRec_cubic(int n,
		    double (*costFunction)(double*),
		    void (*gradCostFunction)(double*,double*),
		    double* point, 
		    double* direction,
		    double h0, 
		    double dh0, 
		    double *tg, 
		    double *tf, 
		    double *td,
		    double m1,
		    double m2,
		    int *nI,
		    int *nE,
		    double *t_,
		    double *h_,
		    double *dh_,
		    int verbosity);
//the recursive Wolfe's method used in Wolfe_cubic.


double Wolfe_cubic( int n,
		    double (*costFunction)(double*),
		    void (*gradCostFunction)(double*, double*),
		    double* point, 
		    double* direction,
		    double *fx0,
		    int *totalCallsJ,
		    int *totalCallsGradJ,
		    int verbosity);
//returns the step given by the line-search

#endif
