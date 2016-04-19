//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//


/* 
** FILE
** *****************************************************************
** 
** linesearch.c
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lineSearch.h"


/* 
** FUNCTION
** *****************************************************************
** 
** meriteFunction
** 
** DESCRIPTION
** 
** Given a cost-function F, a point p and a direction d, 
** computes F(p+td)
** 
*****************************************************************
*/

double meriteFunction( int n, //dimension of the problem
		       double (*costFunction)(double*),    //cost function F
		       double* point,  //previous point of the optimisation algorithm
		       double* direction, //direction of the line search 
		       double t) //step t
{
  double result=0;
  double* pp= (double*) malloc(n*sizeof(double));
  int i;

  for (i=0; i<n; i++)
    {
      pp[i]=point[i]+t*direction[i];
    }
  result=costFunction(pp);
  free(pp);
  return result;
}


/* 
** FUNCTION
** *****************************************************************
** 
** difMeriteFunction
** 
** DESCRIPTION
** 
** Given the gradient of a cost-function F, a point p and a direction d,
** computes the derivative of F(p+td) with respect to the variable t.
** 
** 
** 
*****************************************************************
*/

double difMeriteFunction( int n,
			  void (*gradCostFunction)(double*, double*),
			  double* point, 
			  double* direction, 
			  double t)
{
  double* pp= (double*) malloc(n*sizeof(double));
  double* grad= (double*) malloc(n*sizeof(double));
  
  double dif=0;
  int i;
  for (i=0; i<n; i++)
    {
      pp[i]=point[i]+t*direction[i];
    }
  gradCostFunction(pp,grad);
  for (i=0; i<n; i++)
    {
      dif+= direction[i]*grad[i];
    }
  free(pp);
  free(grad);
  return dif;
}


/* 
** FUNCTION
** *****************************************************************
** 
** cubicMin
** 
** DESCRIPTION
** 
** Given the doubles t1, t2, q1, q2, dq1, dq2,
** computes the argmin of a cubic function f such that:
**    - f(t1)=q1, f'(t1)=dq1
**    - f(t2)=q2, f'(t2)=dq2
** 
*****************************************************************
*/

double cubicMin(double t1, 
		double t2, 
		double q1, 
		double q2, 
		double dq1, 
		double dq2)
{
  double e     = t1  - t2;
  double qq    = (q1-q2)/e;
  double p     = dq1 + dq2 - 3*qq;
  double a     = dq1 + dq2 + 2*p;
  double b     = dq1 + p;
  double discr = p*p - dq1*dq2;
  double d     = 0;

  if (discr>=0.)
    {
      d= sqrt(discr);
    }
  else {return (t1+t2)/2;}
  if (e>0)
    {
      if (b<=0)
	{
	  if (a!=0) {return t1 + e*((-b)+d)/a;} else {return (t1+t2)/2;}
	}
      else
	{
	  return t1 + e*(-dq1)/(b+d);
	}
    }
  else
    {
      if (b>=0)
	{
	  if (a!=0) {return t1 + (-e)*(b+d)/a;} else {return (t1+t2)/2;}
	}
      else
	{
	  return t1 + e*(dq1)/((-b)+d);
	}
    }
}

void testcubicMin()
{
  printf("cubicMin(0,1,1,1,-1,2)=%f\n",cubicMin(0,1,1,1,-1,2));
  printf("cubicMin(0,1,1,1,-1,1)=%f\n",cubicMin(0,1,1,1,-1,1));
}


/* 
** FUNCTION
** *****************************************************************
** 
** WolfeRec_cubic
** 
** DESCRIPTION
** 
** Recursive Wolfe's method with cubic interpolation.
** 
** 
** 
** 
*****************************************************************
*/

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
		    int *callsJ,
		    int *callsGradJ,
		    double *t_,
		    double *h_,
		    double *dh_,
		    int verbosity)
{
  double htf = meriteFunction( n,costFunction,point,direction,*tf);
  double dhtf= difMeriteFunction( n,gradCostFunction,point,direction,*tf);
  (*callsJ)++;(*callsGradJ)++;
      
  if (verbosity>=3) {
    if (verbosity>=3) {printf("\t* t values:\n");}
    printf("\t\ttg=%e\n",*tg);
    if (*td==-1) {printf("\t\ttd not yet defined\n");}
    else {printf("\t\ttd=%e\n",*td);}
    printf("\t\ttf=%e : ",*tf);
    printf("h(tf)=%e, h'(tf)=%e\n",htf,dhtf);
  }
  
  if ( htf > ( h0+m1*(*tf)*dh0 ) && (*callsJ<100))
    {
      *td=*tf;
      *tf=cubicMin(*td,*t_,htf,*h_,dhtf,*dh_);
      printf("\t\tt*=%e\n",*tf);
      if (verbosity>=3) {printf("\t     ->interpolation(td decreased).");}
    
      if ( *tf > (*td-0.5*(*td-*tg)) )
	{
	  *tf= *td-0.5*(*td-*tg);
	  if (verbosity>=3) {printf(" (tf= td-0.5*(td-tg) is used)");}
	}
      if ( *tf < (*tg+0.01*(*td-*tg)) ) 
	{ 
	  *tf= *tg+0.01*(*td-*tg);
	  if (verbosity>=3) {printf(" (tf= tg+0.01*(td-tg) is used)");}
	}
      if (verbosity>=3) {printf("\n");}
     
      WolfeRec_cubic(n,costFunction,gradCostFunction,point,direction,h0,dh0,tg,tf,td,m1,m2,callsJ,callsGradJ,t_,h_,dh_,verbosity);
    }
  else
    {
      if ( dhtf < m2*dh0 && (*callsJ<100) && (*callsGradJ<100))
	{
	  *tg=*tf;
	  if (*td==-1.)
	    {
	      *tf=cubicMin(*tg,*t_,htf,*h_,dhtf,*dh_);
	      *t_  =  *tg;
	      *h_  = htf;
	      *dh_  = dhtf;
	      printf("\t\tt*=%e\n",*tf);
	      if ( *tf < (2.*(*t_)) ) { *tf=2.*(*t_); } 
	      if (verbosity>=3) {printf("\t     ->extrapolation.\n");}
	      
	      WolfeRec_cubic(n,costFunction,gradCostFunction,point,direction,h0,dh0,tg,tf,td,m1,m2,callsJ,callsGradJ,t_,h_,dh_,verbosity);
	    }
	  else
	    {
	      *tf=cubicMin(*tg,*t_,htf,*h_,dhtf,*dh_);
	      *t_  =  *tg;
	      *h_  = htf;
	      *dh_  = dhtf;
	      printf("\t\tt*=%e\n",*tf);
	      if (verbosity>=3) {printf("\t     ->interpolation(tg increased).");}
	      if ( *tf > (*td-0.01*(*td-*tg)) ) 
		{ 
		  *tf= *td-0.01*(*td-*tg); 
		  if (verbosity>=3) {printf(" (tf= td-0.01*(td-tg) is used)");}
		}
	      if ( *tf < (*tg+0.5*(*td-*tg)) ) 
		{ 
		  *tf= *tg+0.5*(*td-*tg); 
		  if (verbosity>=3) {printf(" (tf= tg+0.5*(td-tg) is used)");}
		}
	      if (verbosity>=3) {printf("\n");}
	      
	      WolfeRec_cubic(n,costFunction,gradCostFunction,point,direction,h0,dh0,tg,tf,td,m1,m2,callsJ,callsGradJ,t_,h_,dh_,verbosity);
	    }
	}
      else return;
    }
}


/* 
** FUNCTION
** *****************************************************************
** 
** Wolfe_cubic
** 
** 
** DESCRIPTION
** 
** Wolfe's cubic line-search.
** 
** 
** 
** 
*****************************************************************
*/


double Wolfe_cubic( int n,
		    double (*costFunction)(double*),
		    void (*gradCostFunction)(double*,double*),
		    double* point, 
		    double* direction,
		    double *fx0,
		    int *totalCallsJ,
		    int *totalCallsGradJ,
		    int verbosity )
{
  //parameters of the Wolfe's algorithm
  double m1=0.0001;
  double m2=0.99;
  
  double h0 = meriteFunction( n,costFunction,point,direction,0 );
  double dh0= difMeriteFunction( n,gradCostFunction,point,direction,0 );
  
  //initialisations
  double tg=0.;
  double tf=1.;
  double td=-1.;
  double t_=0;
  double h_=h0;
  double dh_=dh0;
  int callsJ=1;
  int callsGradJ=1;

  *fx0=h0;//value of the Cost-Function used in QuasiNewton
  
  if (verbosity>=1) {printf("J = %e\n", h0);}
  if (verbosity>=2) {printf("Line Search:\n ", h0, dh0);}
  if (verbosity>=3) {printf("\th(0)=%e, h'(0)=%e\n", h0, dh0);}


  //Calls the recursive function WolfeRec_cubic
  WolfeRec_cubic( n,costFunction,gradCostFunction,point,direction,h0,dh0,&tg,&tf,&td,m1,m2,&callsJ,&callsGradJ,&t_,&h_,&dh_,verbosity);
  
  (*totalCallsJ)+=callsJ;
  (*totalCallsGradJ)+=callsGradJ;
  
  if (verbosity>=3) {printf("\tfinal step=%e\n",tf);}
  if (verbosity>=2) {printf("\t");}
  if (verbosity>=3) {printf("End of the Line Search:");} 
  if (verbosity>=2) {printf(" #calls J:%d  gradJ:%d\n",callsJ,callsGradJ);}

  //the final result is the step given by the line-search
  return tf; 
}
