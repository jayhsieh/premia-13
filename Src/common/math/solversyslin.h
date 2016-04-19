#ifndef SOLVERSYSLIN_H
#define SOLVERSYSLIN_H
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

void multiplytridiag(double **M, double* u, double *r, int n);
void tridiagsolve(double **M, double *u,  double *r, int n);

#endif    

