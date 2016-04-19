#ifndef NORMAL_DF
#define NORMAL_DF
#include <math.h>
#include <iostream>
#include <cstdlib>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

double gammln(double xx);

void gser(double *gamser, double a, double x, double *gln);/* Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser. Also returns ln\Gamma(a) as gln*/
void gcf(double *gammcf, double a, double x, double *gln); /*Returns the incomplete gamma funcion Q(a, x) evaluated by its continued fraction representation as gammcf. Also returns ln\Gamma(a) as gln.*/

double gammp(double a, double x); /* returns the imcomplete gamma function P(a, x)*/

double gammq(double a, double x); /*Returns the imcomplete function Q(a, x)=1-P(a, x).*/

double erff(double x);/* Return the error function erf(x)*/

double erffc(double x);

double normal_df(double x);

#endif


