/************************************************
 *   Bermudean Swaptions Pricer  (Andersen-MC in one-factor-LMM)           
 *
 *    Sonke Blunck,  2005
 *
 ****************************************************/ 
 

#include<iostream>

using namespace std;

#include "lmm_bermudaprice_andersen.h"



int main()
{

  double p;
  float tenor=0.5;               /*tenor is the lenght of the rate usually 3 months or 6 months */
  float swaptionMat=1.0;        //(years)
  float swapMat=4.0;            //(years)
  float K=0.06;                 //strike
  float flatInitialValue=0.06;  //for the SDE of the LIBORs
  float vol=0.2;                //for the SDE of the LIBORs
  long numberMCPaths1=10000;    //for the comput. of the parameters
  long numberMCPaths2=50000;    //for the comput. of the price
  int AndersenStrategy=1;       //must be 1 or 2 

  srand48(1);

  p=lmm_swaption_payer_bermudan_andersen_pricer(tenor, swaptionMat, swapMat, K, flatInitialValue, vol, numberMCPaths1, numberMCPaths2, AndersenStrategy); 

  printf("Bermudean swaption with Andersen algorithm is %f\n",p*10000);
 
  return(EXIT_SUCCESS);
}







