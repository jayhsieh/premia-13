#include "maths.h"

/** trans is a matrix of transition rates, nb number of discretisation points, M number of firms , T maturity ,t the time of evalaution of the spread**/

/** this function simulate the number of default between 0 and t**/

int simul_perte(double t,double T,int M,int nb,double **trans);

/** Simulate the probability of transition see Schonbucher paper**/

double  ***proba(double t,double T,int M,int nb,double **trans); 
