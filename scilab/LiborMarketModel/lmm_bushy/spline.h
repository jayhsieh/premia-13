#ifndef SPLINE_H
#define SPLINE_H
#include <memory.h>
#include <string.h>
#include <stdio.h>
#include<stdlib.h>
#include<math.h>

double cubspline(int Metode, double xi, double *xx,double *yy, int N);
void spline(double *x, double *y, int N,double  yp1, double ypn ,double *y2);
void splint(double *xa, double *ya,double  *y2a,int N ,double x, double *y);
double SplineX3(double x , double *xx, double *yy, int N);
double dxx(double x1, double x0);
int Sort(double *x,double*y, int size);

#endif
