
#ifndef _LIBOR_AFFINE_PRICING_H
#define _LIBOR_AFFINE_PRICING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "math/read_market_zc/InitialYieldCurve.h"

double cf_swaption_fourier_libaff(StructLiborAffine *LiborAffine, double first_reset_date, double contract_maturity, double period, double Nominal, double swaption_strike, int swaption_payer_receiver);

// Value of the swaption when strike=0
double SwapValue(double T_start, double T_end, double period, double strike, double nominal, ZCMarketData* ZCMarket);

#endif
