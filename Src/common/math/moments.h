#ifndef  _MOMENTS_H
#define _MOMENTS_H

#define NGAMMA 300
#define EPS_MOMENT 3.e-14

void gauleg(double x1, double x2, double x[], double w[], int n);
double gammadensity(double x, double a, double b);
double factrl(int n);
double bico(int n, int k);
double factln(int n);
double Moments(int n,double r,double sigma,double t);
double logdens(double x,double m,double sg);
double Der1Logdens(double x, double m, double sg);
double Der2Logdens(double x, double m, double sg);
double Der3Logdens(double x, double m, double sg);
double Der4Logdens(double x, double m, double sg);
double momlog(int n, double mean, double var);
double Normdens(double x, double m, double sg);
double Der1Normdens(double x, double m, double sg);
double Der2Normdens(double x, double m, double sg);
double Der3Normdens(double x, double m, double sg);
double Der4Normdens(double x, double m, double sg);

#endif
