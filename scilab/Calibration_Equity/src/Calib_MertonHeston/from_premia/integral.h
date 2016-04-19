#ifndef  _INTEGRAL_H
#define _INTEGRAL_H
#include "optype.h"
#define FACTOR  1.6
#define NTRY     50
#define JMAX     40
#define EPS 1.0-6
#define FJMAX 14
#define JMAXP (FJMAX+1)
#define KK 5

double *vector(long nl,long nh);
void free_vector(double *v,long nl,long nh);
void polint(double xa[],double ya[],int n, double x,double *y,double *dy);
double midpnt(double (*func)(double), double a, double b, int n);
double midpntbis(double (*func)(double), double a, double b, int n);
double midsql(double (*func2)(double), double aa, double bb, int n);
double midsqlbis(double (*func2)(double), double aa, double bb, int n);
double midinf(double (*func1)(double), double aa, double bb, int n);
double midinfbis(double (*func1)(double), double aa, double bb, int n);
double qromo(double(*func)(double),double a,double b,double(*choose)(double(*)(double),double,double,int));
void polint(double xa[],double ya[],int n, double x,double *y,double *dy);
 #endif
