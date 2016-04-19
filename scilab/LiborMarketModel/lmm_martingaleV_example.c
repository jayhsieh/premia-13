/**********************************************************************************************
 *
 *
 *
 *   Example: using Martingale V approach to price caplets and swaptions
 *
 *
 *
 *
 *********************************************************************************************/
#include<stdio.h>

#include"lmm_martingaleV.h"


int main()
{


  double  tenor=0.5;   // period (in years) of the rate usually 3 or 6 months
  int numMat;     
  double *caplets;
  int i;
  double* maturities;
  double swaption_price;
  double swaption_maturity=3;  // swaption maturity in years
  double swap_maturity=7.;     // swap maturity in years
  int nb_MC   =10000;          // number of monte carlo paths
  int nb_time_step=10;         // number of time steps in the euler scheme
  double  strike=0.05;         // strike 
  int nb_factors=2;            // number of factors

  
  
  numMat=(int)(swap_maturity/tenor);
  maturities=(double*)malloc((numMat+1)*sizeof(double));
  caplets=(double*)malloc((numMat+1)*sizeof(double));
  lmm_caplet_spotV_pricer( caplets ,  maturities , numMat , nb_MC , nb_factors , 
			   nb_time_step , strike ,  tenor);
  lmm_swaption_payer_spotV_pricer(&swaption_price ,  swaption_maturity , swap_maturity   , 
				  nb_MC ,  nb_factors ,  nb_time_step ,  strike , tenor);

  for(i=1;i<numMat;i++)
    {
      printf("Caplet(Ti=%lf,%lf,K=%lf)=%lf \n", tenor*i ,tenor*(i+1) , strike , caplets[i]);
    }
  printf("Spt(T=%lf,%lf,K=%lf)=%lf \n",swaption_maturity, swap_maturity, strike, swaption_price);
  free(caplets);
  free(maturities);

  return(1);
}
