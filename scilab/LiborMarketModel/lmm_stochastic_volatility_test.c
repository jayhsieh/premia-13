/*********************************************************************************************************************************************************
 *
 *
 *
 *
 *
 *   Implementation of the article 
 *
 *          Title: "Libor Market Model: From deterministic to stochastic volatility" 
 *          From:  Lixin WU, Fan ZHANG 
 *          Working paper: version October 4, 2002  
 *
 *   Author: Jose Da Fonseca
 *   Date  : 2004
 *
 *
 *
 *********************************************************************************************************************************************************/


#include"stdio.h"
#include"math.h"

#include"lmm_stochastic_volatility.h"

int main()
{

  float  tenor=0.5;
  int numFac=2;
  int numMat=27;
  double  swaptionMat=5.;
  double  swapMat=10.;
  double  priceVal=0.25;
  double atm_Strike;
  int s, i,M;
  double percent=-20.;


  printf(" Payer swaption with maturity: %lf \n",swaptionMat );
  printf(" on a swap rate with maturity: %lf  (tenor equal to %lf)  \n", swapMat , swaptionMat - swapMat );
  printf(" the strike is equal to (1+%lf) of the ATM strike \n", percent/100.); 
  printf(" the period of the underlying libor rate is %lf \n" , tenor );
  printf(" the price in bps is : %lf \n", lmm_swaption_payer_stoVol_pricer(tenor ,numFac ,swaptionMat , swapMat , percent)*10000 ) ;

  
  return(1);

} 
