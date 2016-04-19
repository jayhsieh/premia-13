#include "spread_dynamique.h"


double* initial_rates_CDS(double T,int nb,double r,double *spread);

/** caibration of transition rates a_n[0,s] ,0<s<10 , 0<n<125 **/

double *initial_rates_CDO(double T,int n,double r,double *spread);

/**simulation dynamics path on (0,t) of the transition knowing a_n[0,T],the drift term should be satisfy an equation given in the Shonbucher paper and the volatility should be given by the user**/

double*** calib_rates_CDS(double t,double T,int nb,double vol,double r,double* spread);
				
double*** calib_rates_CDO(double t,double T,int nb,double vol,double r,double* spread);

