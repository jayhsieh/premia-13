#include <math.h> 
/*
Fast integration routine
*/  
void intg(double a, double b, double (*f)(double), double ea, double er, double *val, double *abserr);
/*
Parameters:
a,b -- interval of integration
f -- integrand 
ea, er -- desired absolute and relative accuracy

Return:
val -- value of the integral
abserr -- estimate of the real absolute error 
*/
