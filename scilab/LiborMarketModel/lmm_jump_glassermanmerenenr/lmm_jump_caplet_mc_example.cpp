/*--------------------------------------------------------------*/
/*   CF approx. for caplet prices in one-factor LMM with jumps  */
/*   Algorithm of Glasserman/Merener                            */
/*                                                              */
/*--------------------------------------------------------------*/
/*  Sonke Blunck, Premia 2005                                   */
/*--------------------------------------------------------------*/


#include<stdio.h>
#include<stdlib.h>

#include "lmm_jump_capletprice_glassmer.h"



int main()
{

  double p;
  float tenor=0.5;               /*usually 3 months or 6 months */
  float capletMat=2.0;          //(years)
  float K=0.07;                 //strike
  float flatInitialValue=0.06;  //for the SDE of the LIBORs
  float vol=0.1;                //for the SDE of the LIBORs   
  int numberMCPaths=1000;       

  srand48(1);

  p=lmm_jump_caplet_MC_pricer(tenor, capletMat, K, flatInitialValue, vol, numberMCPaths ); 

  printf("Caplet with Monte Carlo is %f\n",p*10000);



  return(EXIT_SUCCESS);
}

