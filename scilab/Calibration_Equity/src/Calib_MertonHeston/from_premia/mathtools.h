#ifndef  _MATHTOOLS_H
#define _MATHTOOLS_H
#include "optype.h"

typedef unsigned char boolean;

#define false 0
#define true 1

#define PI 3.14159265358979

#define BIG_DOUBLE 1.0e6
#define PRECISION 1.0e-7 /*Precision for the localization of FD methods*/
#define INC 1.0e-5 /*Relative Increment for Delta-Hedging*/

#define MAXLOOPS 5000
#define ABS(x) (((x) < 0) ? -(x) : (x))
#define POW(x,y) pow( (double) (x), (double) (y))
#define MAX(A,B) ( (A) > (B) ? (A):(B) )
#define MIN(A,B) ( (A) < (B) ? (A):(B) )
#define SQR(X) ((X)*(X))
#define CUB(X) ((X)*(X)*(X))
double dabs(double x);
double nd(double x);
double N(double x);
double erf(double x);
double NN( double a, double b, double r);
void Sort(unsigned long n, double *arr);
int nn(int l);
#endif
