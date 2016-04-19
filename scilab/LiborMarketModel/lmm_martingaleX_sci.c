/*************************************************************************************
 *
 *
 *
 *   Example: interface to scilab function 
 *          
 *
 *  
 *
 *************************************************************************************/
#include<stdio.h>

#include"lmm_martingaleX.h"

int lmm_cap_martX_sci(double *price,  double* nb_fac, double *period, double* strike , double* maturity )
{
  
  double  tenor= (*period);
  int numMat;     
  double *caplets;
  double* maturities;
  double swaption_price;
  int nb_MC   =10000;          // number of monte carlo paths
  int nb_time_step=10;         // number of time steps in the euler scheme
  int i;
  
  numMat=(int)( (*maturity)/tenor);
  caplets=(double*)malloc(numMat*sizeof(double));
  maturities=(double*)malloc(numMat*sizeof(double));

  lmm_caplet_terminalX_pricer(caplets, maturities , numMat, nb_MC, *nb_fac, nb_time_step, *strike, *period);
  price[0]=0.0;
  for(i=1;i<numMat;i++)
    {
      price[i]=caplets[i];
    }
  
  free(caplets);
  free(maturities);

  return(1);
}




int lmm_swpt_martX_sci(double *price, double* swpt_mat, double* swap_mat , int* nb_factors, double* strike, double* period )
{
  double  tenor=(*period);  
  double swaption_price;
  int nb_MC   =10000;          // number of monte carlo paths
  int nb_time_step=10;         // number of time steps in the euler scheme
  
  lmm_swaption_payer_terminalX_pricer(price , *swpt_mat , *swap_mat , nb_MC, *nb_factors, nb_time_step, *strike, tenor);
  
return(1);
}
