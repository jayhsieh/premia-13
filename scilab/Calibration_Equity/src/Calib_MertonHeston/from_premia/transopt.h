#ifndef  _TRANSOPT_H
#define _TRANSOPT_H

int Call_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);
int Put_BlackScholes_73 (double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);
int AlgCrayer(int N,double *Z,int  ssl,double *A,double *B,double *C,double*Q,double  *S);
int PutMinAn(double s1,double s2,double k,double t,
             double r,double divid1,double divid2,
             double sigma1,double sigma2,double rho,
             double *ptprice,double *ptdelta1,double *ptdelta2);
#endif
