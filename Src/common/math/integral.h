#ifndef _INTEGRAL_H
#define _INTEGRAL_H
#include "optype.h"
#define FACTOR  1.6
#define NTRY     50
#define JMAX     40

#define FJMAX 14
#define JMAXP (FJMAX+1)
#define KK 5

double integrale_gauss(double (*funcg)(double), double a, double b);
void integrale_gauss_vect(void (*funcg_vect)(double,int,double *), double a, double b, int dimx, double *sum);
void init_gauss(int n);
void free_gauss(void);

#endif
