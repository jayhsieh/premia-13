/************************************************
 *   Bermudean Swaptions Pricer             
 *
 *   Interface to scilab function
 *   
 *
 *
 ****************************************************/ 
 



#include"lmm_bermudaprice.h"

int lmm_bermuda_LS_sci( double*price , double*tenor , int* numFac , double*  swaptionMat , double*  swapMat , double* K , double* payoff_as_Regressor , int* Regr_Basis_Dimension)
{
  
  int numberTimeStep=10;
  long numberMCPaths=10000;
  char Explanatory='N';           //Explanatory variable for regression B=Brownian, 
                                  //                 S=Nominal Swap Paying Value, N=Numeraire
  char* basis_name="HerD";        //Hermite basis
  char* measure_name="Spot" ;     // spot "    "    "    "

  (*price)=lmm_swaption_payer_bermudan_LS_pricer(*tenor ,numberTimeStep, *numFac, *swaptionMat , *swapMat ,*payoff_as_Regressor , numberMCPaths , *Regr_Basis_Dimension , basis_name ,measure_name , Explanatory , *K); 

  /*
  printf("Bermudean swaption with terminal time %f and exercise\n",swapMat);
  printf("each %f year starting from %f ,is, under %s measure\n %f bps.\n",tenor,swaptionMat 
	                                                                    , measure_name , p*10000);
  
  */
  
  
  return(EXIT_SUCCESS);
}







