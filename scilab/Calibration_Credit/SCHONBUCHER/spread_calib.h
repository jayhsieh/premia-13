#include "initial_calibration.h"

/** this function return the path of the spread of the index CDS between (0,t)**/

double *calib_spread_CDS(double t,double T,int M,int nb,double vol,double r,double R,double *spread);

/** this function return the path of the spread of the CDO (a,b) between (0,t)**/

double *calib_spread_CDO(double t,double T,int M,int nb,double vol,double r,double a,double b,double R,double *spread); 
