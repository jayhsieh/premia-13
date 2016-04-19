#include <math.h>
//#include "bs1d_std.h"
#include "parameter.h"



// copied from bs1d_std.h now removed form this file
#include "bs1d.h"
#include "std.h"

#include "mathtools.h"
#include "random.h"
#include "numfunc.h"
#include "transopt.h"
#include "linsys.h"



//
int MyCall_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);
int MyPut_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *ptprice,double *ptdelta);
int dCall_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice);
int dPut_BlackScholes_73(double s,double k,double t,double r,double divid,double sigma,double *dprice);
double SigmaImplicite(double eps,double a,double b,int TypeOpt,double St,double T,double K,
					  double r,double divid,double cbs,int *nbiters);
