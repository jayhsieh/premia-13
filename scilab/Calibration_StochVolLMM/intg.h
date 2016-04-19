#ifndef _INTEGRAL_H
#define _INTEGRAL_H
#include "optype.h"
#define FACTOR  1.6
#define NTRY     50
#define JMAX     40

#define FJMAX 14
#define JMAXP (FJMAX+1)
#define KK 5

/* double *vector(long nl,long nh);
   void free_vector(double *v,long nl,long nh); */
void polint(double xa[],double ya[],int n, double x,double *y,double *dystatic );
double midpnt(double (*func)(double), double a, double b, int n);
double midpntbis(double (*func)(double), double a, double b, int n);
double midsql(double (*func2)(double), double aa, double bb, int n);
double midsqlbis(double (*func2)(double), double aa, double bb, int n);
double midinf(double (*func1)(double), double aa, double bb, int n);
double midinfbis(double (*func1)(double), double aa, double bb, int n);
double qromo(double(*func)(double),double a,double b,double(*choose)(double(*)(double),double,double,int));
double integrale_gauss(double (*funcg)(double), double a, double b);
void integrale_gauss_vect(void (*funcg_vect)(double,int,double *), double a, double b, int dimx, double *sum);
void init_gauss(int n);
void free_gauss(void);

void intg_multi(double a, double b, double (*f)(double []), double ea, double er, double *val, double *abserr, double x[], int n);
double f_new(double x);
void intg_f3(double a, double b, double (*f)(double []), double ea, double er, double *val, double *abserr, double x1, double x2);
double f3_new(double x);
void intg(double a, double b, double (*f)(double), double ea, double er, double *val, double *abserr);

#endif
