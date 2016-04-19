/************************************************
 *   Bermudean Swaptions Pricer             
 *
 *   -Spot Probability Measure (Numeraire=Roll-Over Bond)
 *   -Possibility to include Martingale Discretization
 *   -Bermudean with fixed ending date.
 *   Nicola Moreni, August 2004 
 *
 ****************************************************/ 
 

#include<stdio.h>

#include"lmm_header.h"
#include"lmm_volatility.h"
#include"lmm_libor.h"
#include"lmm_random_generator.h"
#include"lmm_products.h"
#include"lmm_numerical.h"
#include"lmm_zero_bond.h"
#include"lmm_bermudaprice.h"

/*REMINDER:
-lmm_header.h contains structures definition: volatility, libor, swaption
-lmm_vol.c lmm_products.c lmm_libor.c  contain allocation/initialization routines
-lmm_numerical.c contains evolution routines, european swaption pricers
-lmm_mathtools.c contains random number generator, cholesky sqrt....
In the present file file all parametres are given an initial value and pricing routine 
is called
*/
  
int main()
{

  float  tenor=0.5;               // tenor is the lenght of the rate usually 3 months or 6 months 
  int numberTimeStep=10;
  int numFac=1;
  double  swaptionMat=1.0;        //(years)
  double  swapMat=4.0;            //(years)
  double payoff_as_Regressor=3.0; //(years) Maturity after which payoff is included in regression
  double  priceVal=0.20;
  double K=0.05;                  //strike
  long numberMCPaths=10000;
  int Regr_Basis_Dimension=4 ;    //finite-dimensional approx. of L² 
  char Explanatory='N';           //Explanatory variable for regression B=Brownian, 
                                  //                 S=Nominal Swap Paying Value, N=Numeraire;
  char* basis_name="HerD";        //Hermite basis
  char* measure_name="Spot" ;     // spot "    "    "    "
  double p;

  
  
  p=lmm_swaption_payer_bermudan_LS_pricer(tenor ,numberTimeStep, numFac, swaptionMat , swapMat ,  
                   payoff_as_Regressor , numberMCPaths , Regr_Basis_Dimension , basis_name ,
                   measure_name , Explanatory , K); 

  printf("Bermudean swaption with terminal time %f and exercise\n",swapMat);
  printf("each %f year starting from %f ,is, under %s measure\n %f bps.\n",tenor,swaptionMat 
	                                                                    , measure_name , p*10000);
  
  
  
  return(EXIT_SUCCESS);
}







