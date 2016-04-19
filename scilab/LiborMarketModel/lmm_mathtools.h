#ifndef LMM_MATHTOOLS_H
#define LMM_MATHTOOLS_H

#define FACTOR  1.6
#define NTRY     50
#define JMAX     40
#define EPS 1.0-6
#define FJMAX 14
#define JMAXP (FJMAX+1)
#define KK 5

double N(double x);
double integrale_gauss(double (*funcg)(double), double a, double b);
void init_gauss(int n);
void free_gauss(void);
int Cholesky(double *M, int dim);
void Resolution(double *AuxR, double *Res, double *M, int dim);
#endif
