#ifndef  _TRANSOPT_H
#define _TRANSOPT_H

#include "pnl/pnl_finance.h"

int AlgCrayer(int N,double *Z,int  ssl,double *A,double *B,double *C,double*Q,double  *S);
int PutMinAn(double s1,double s2,double k,double t,
             double r,double divid1,double divid2,
             double sigma1,double sigma2,double rho,
             double *ptprice,double *ptdelta1,double *ptdelta2);
int iac_kobol_europut(int ifCall, double lm, double lp, 
			 double alm, double alp, double cm, double cp, 
			 double r, double T, double Strike,
			 double Spot, double eps, double *Price);
#endif
