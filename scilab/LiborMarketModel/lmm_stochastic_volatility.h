

#ifndef LMM_STOCHASTIC_VOLATILITY_H 
#define LMM_STOCHASTIC_VOLATILITY_H 

#include"stdio.h"
#include"math.h"

#include"lmm_header.h"
#include"lmm_random_generator.h"
#include"lmm_volatility.h"
#include"lmm_products.h"
#include"lmm_libor.h"
#include"lmm_mathtools.h"
#include"complex.h"

double lmm_swaption_payer_stoVol_pricer(double period , int number_of_factors ,double swaption_maturity , double swap_maturity , double percent_of_ATM_strike );

#endif
