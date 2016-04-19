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
#include <string.h>

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
  Volatility *ptVol;
  Libor *ptLib;
  Swaption *ptSwpt;

  float  tenor=0.5;
  int numMat=16;       //WARNING: including T_0=0, maturities number is 'numMat+1'
  int numberTimeStep=10;
  int numFac=1;
  double  swaptionMat=1.0;//(years)
  double  swapMat=8.0;//(years)
  double payoff_as_Regressor=5.0;//(years) Maturity after which payoff is included in regression
  double  priceVal=0.20;
  double K=0.05;//strike
  long numberMCPaths=10000;
  int Regr_Basis_Dimension=4; //finite-dimensional approx. of L² 
  char Explanatory='N';//Explanatory variable for regression B=Brownian, S=Nominal Swap Paying Value, N=Numeraire;
  char Basis_Choice[10];  //See in the followings
  char Measure_Choice[10];// "    "    "    "

  
  char auxstr[5];
  int I;
  double DT=(tenor/(double)numberTimeStep);
  
 
  strcpy(Basis_Choice,"HerD"); //"CanD" for canonical basis or "HerD" for hermite 
  strcpy(Measure_Choice,"Spot");//"Spot" or "Fwd" cfr numerical.c for evolution
 
  //Memory Allocation and Libor/Swap Structures Initialization
  mallocLibor(&ptLib,numMat,tenor);
  mallocVolatility(&ptVol,numFac);
  mallocSwaption(&ptSwpt,swaptionMat,swapMat,priceVal,K,tenor);
  
  //Swaption Computation
  I=computeBermudeanSwaption(numberMCPaths,numberTimeStep,ptLib,ptSwpt,ptVol,Regr_Basis_Dimension, payoff_as_Regressor,Basis_Choice,Measure_Choice,Explanatory);
  

  //Printing Results to STDOutput
  //printf("s=%d, e=%d, exdates=%d, payoffasregr=%d\n",s,ptLib->numberOfMaturities,numberOfExerciseDates,PayOff_As_Regressor);
  
  printf("Bermudean swaption with terminal time %f and excercise\n",ptSwpt->swapMaturity);
  printf("each %f year starting from %f ,is, under %s measure\n %f bps.\n",ptLib->tenor,ptSwpt->swaptionMaturity,Measure_Choice,ptSwpt->price*10000);
  //initLibor(ptLib);
  //I=computeSwaptionPriceSpot(numberMCPaths,numberTimeStep,DT,ptLib,ptSwpt,ptVol);
  
  //freeing Memory
  freeSwaption(&ptSwpt);
  freeLibor(&ptLib);
  freeVolatility(&ptVol);
  
  return(EXIT_SUCCESS);
}







