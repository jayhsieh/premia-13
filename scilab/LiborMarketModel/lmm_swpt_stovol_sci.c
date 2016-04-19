/*********************************************************************************************************************************************************
 *
 *  Interface with scilab
 *           
 *
 *
 *********************************************************************************************************************************************************/


#include"stdio.h"
#include"math.h"

#include"lmm_stochastic_volatility.h"

void lmm_swpt_stovol_sci(double *price , double* tenor , int* numFac , double* swaptionMat , double* swapMat , double* percent)
{

  *price=lmm_swaption_payer_stoVol_pricer(*tenor ,*numFac ,*swaptionMat , *swapMat , *percent);
  
} 






 

