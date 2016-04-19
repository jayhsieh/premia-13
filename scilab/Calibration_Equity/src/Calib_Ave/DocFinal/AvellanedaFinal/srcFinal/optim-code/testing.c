//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//


/* 
** FILE
** *****************************************************************
** 
** testing.c
** 
** 
** PURPOSE
** 
** Some test-problems in order to test optimizer.
** See:
**      - Schittkowski, "More Test Examples for Nonlinear Programming Codes"
** (Springer-Verlag, 1987).
**      - J. E. Dennis and R.B. Schnabel, "Numericals Methods for Unconstrained
** Optimization and Nonlinear Equations" (Prentice-Hall, 1983)
** 
*****************************************************************
*/

#include <stdio.h>
#include <math.h>

#include "testing.h"

int i;



/* 
** TEST-PROBLEM
** *****************************************************************
** 
** Rosenbrock
** 
**
*****************************************************************
*/

double RosenStart[2]={-1.2,1};
double RosenSol[2]={1,1};

double Rosen(double* p)
{
  double x=p[0];
  double y=p[1];
  return pow((x-1),2) + 10.*pow((x*x-y),2);
}

void RosenGrad(double* p, double* grad)
{
  double x=p[0];
  double y=p[1];
  grad[0]=2*(x-1)+ 40*x*(x*x-y);
  grad[1]=-20*(x*x-y);
}

const testProblem Rosenbrock={2,RosenStart,Rosen,RosenGrad,RosenSol};



/* 
** TEST-PROBLEM
** *****************************************************************
** 
** Rosenbrock100
** 
**
*****************************************************************
*/

double Rosen100Start[2]={-1.2,1};
double Rosen100Sol[2]={1,1};

double Rosen100(double* p)
{
  double x=p[0];
  double y=p[1];
  return pow((x-1),2) + 100.*pow((x*x-y),2);
}

void Rosen100Grad(double* p, double* grad)
{
  double x=p[0];
  double y=p[1];
  grad[0]=2*(x-1)+ 400*x*(x*x-y);
  grad[1]=-200*(x*x-y);
}

const testProblem Rosenbrock100={2,Rosen100Start,Rosen100,Rosen100Grad,Rosen100Sol};



/* 
** TEST-PROBLEM
** *****************************************************************
** 
** Quadratic
** 
*****************************************************************
*/

double QuadStart[2]={8,9};
double QuadSol[2]={1,1};

double Quad(double* p)
{
  double x=p[0];
  double y=p[1];
  return 4*pow((x-5),2) + pow((y-6),2);
}

void QuadGrad(double* p, double* grad)
{
  double x=p[0];
  double y=p[1];
  grad[0]=8*(x-5);
  grad[1]=2*(y-6);
}

const testProblem Quadratic={2,QuadStart,Quad,QuadGrad,QuadSol};


/* 
** TEST-PROBLEM
** *****************************************************************
** 
** Quadratic2
** 
*****************************************************************
*/

double Quad2Start[3]={100,-1,2.5};
double Quad2Sol[3]={0,0,0};

double Quad2(double* p)
{
  double x=p[0];
  double y=p[1];
  double z=p[2];
  return pow((x-y+z),2) + pow((-x+y+z),2) + pow((x+y-z),2);
}

void Quad2Grad(double* p, double* grad)
{
  double x= p[0];
  double y= p[1];
  double z= p[2];
  grad[0] = 2*(x-y+z)+2*(x+y-z)-2*(-x+y+z);
  grad[1] = 2*(-x+y+z)+2*(x+y-z)-2*(x-y+z);
  grad[2] = 2*(x-y+z)+2*(-x+y+z)-2*(x+y-z);
}

const testProblem Quadratic2={3,Quad2Start,Quad2,Quad2Grad,Quad2Sol};


/* 
** TEST-PROBLEM
** *****************************************************************
** 
** WoodProblem
** 
*****************************************************************
*/

double WoodStart[4]={-3,-1,-3,-1};
double WoodSol[4]={1,1,1,1};

double Wood(double* p)
{
  double x=p[0];
  double y=p[1];
  double z=p[2];
  double t=p[3];
  return 
    100*pow(x*x-y,2)+
    pow(1-x,2)+
    90* pow(z*z-t,2)+
    pow(1-z,2)+
    10*pow(y+t-2,2)+
    0.1*pow(y-t,2);
}

void WoodGrad(double* p, double* grad)
{
  double x= p[0];
  double y= p[1];
  double z= p[2];
  double t= p[3];
  grad[0] = 2 * ( 200*x*(x*x-y) + (1-x) );
  grad[1] = 2 * ( 100*(y-x*x) + 10*(y+t-2) + 0.1*(y-t) );
  grad[2] = 2 * ( 180* z*(z*z-t) + (z-1) ); 
  grad[3] = 2 * ( 90* (t-z*z) + 10* (y+t-2) + 0.1*(t-y) ); 
}

const testProblem WoodProblem={4,WoodStart,Wood,WoodGrad,WoodSol};


/* 
** TEST-PROBLEM
** *****************************************************************
** 
** GaussProblem
** 
*****************************************************************
*/

double GaussStart[30]=
{
-1.03, 1.07, -1.10, 1.13, -1.17, 1.20, -1.23, 1.27, -1.30, 1.33, 
-1.37, 1.40, -1.43, 1.47, -1.50, 1.53, -1.57, 1.60, -1.63, 1.67,
-1.70, 1.73, -1.77, 1.80, -1.83, 1.87, -1.90, 1.93, -1.97, 2.00
};
double GaussSol[30]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

double Gauss(double* p)
{
  double s=0;
  int i;

  for (i=0;i<30;i++) {s+=p[i]*p[i];}
  return 1-exp(-s/60.);
}

void GaussGrad(double* p, double* grad)
{
  
  double s=0;
  int i;

  for (i=0;i<30;i++) {s+=p[i]*p[i];}
  s = 1-exp(-s/60.);
  for (i=0;i<30;i++) { grad[i] = s*p[i]/30.;}
}

const testProblem GaussProblem={30,GaussStart,Gauss,GaussGrad,GaussSol};



/* 
** TEST-PROBLEM
** *****************************************************************
** 
** BigGaussProblem
**
** You have to run fillBigGaussSol() and fillBigGaussStart() before
** testing the problem BigGaussProblem.
** 
*****************************************************************
*/

double BigGaussStart[1000];

void fillBigGaussStart(void){
for (i=0;i<1000;i++) { BigGaussStart[i] = pow(-1,i+1)*(1+i/1000.);}
}

double BigGaussSol[1000];

void fillBigGaussSol(void){
for (i=0;i<1000;i++) { BigGaussStart[i] = 0.;}
}

double BigGauss(double* p)
{
  double s=0;
  int i;

  for (i=0;i<1000;i++) {s+=p[i]*p[i];}
  return 1-exp(-s/2000.);
}

void BigGaussGrad(double* p, double* grad)
{
  
  double s=0;
  int i;

  for (i=0;i<1000;i++) {s+=p[i]*p[i];}
  s = 1-exp(-s/2000.);
  for (i=0;i<1000;i++) { grad[i] = s*p[i]/1000.;}
}

const testProblem BigGaussProblem={1000,BigGaussStart,BigGauss,BigGaussGrad,BigGaussSol};
