#ifndef _BSVANILLAS_H_
#define _BSVANILLAS_H_

#include "utils1.h"

int dCall_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice);
int dPut_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice);
double SigmaImplicite(int TypeOpt,double St,double T,double K, double r,double divid,double cbs,int *nbiters);
//
#endif
