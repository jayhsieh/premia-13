#include "tryschon.h"

/** t time of evaluation of the spread, T maturity, M number of firms, nb number of subdivision, r interest rate, R recovery rate, trans matrix of transition rate**/

/** this function give the spread of the index CDS at time t knowing the matrix of transition rates**/

double spread_CDS(double t,double T,int M,int nb,double r,double R,double **trans);

/** this function give the spread of the CDO's tranche (a,b)  at time t knowing the matrix of transition rates**/

double spread_CDO(double t,double T,int m,int nb,double r,double a,double b,double R,double **trans); 

