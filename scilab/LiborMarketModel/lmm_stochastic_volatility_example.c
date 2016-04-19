/*********************************************************************************************************************************************************
 *
 *  example: - using the pricer function in a program
 *           - stochastic volatility model
 *
 *
 *********************************************************************************************************************************************************/


#include"stdio.h"
#include"math.h"

#include"lmm_stochastic_volatility.h"

int main()
{

  float  tenor=0.5;         // period of the rate usually 3 months or 6 months
  int numFac=2;             // number of factors: dim of the brownian motion of the rates   
  double  swaptionMat=5.;   // swaption maturity
  double  swapMat=10.;      // swap maturity -> the tenor of swaption is  swapMat - swaptionMat
  double percent=-20.;      // the strike will be equal to (1+percent/100)*atm_strike
  double price;

  printf(" Payer swaption with maturity: %lf \n",swaptionMat );
  printf(" on a swap rate with maturity: %lf  (tenor equal to %lf)  \n", swapMat , 
	                                                   swapMat - swaptionMat);
  printf(" the strike is equal to (1+%lf) of the ATM strike \n", percent/100.); 
  printf(" the period of the underlying libor rate is %lf \n" , tenor );
  price=lmm_swaption_payer_stoVol_pricer(tenor ,numFac ,swaptionMat , swapMat , percent);
  printf(" the price in bps is : %lf \n", price*10000 ) ;

  
  return(1);

} 






 

