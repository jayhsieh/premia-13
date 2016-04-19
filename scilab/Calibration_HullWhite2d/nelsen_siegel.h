#ifndef NELSENSIEGEL_H_INCLUDED
#define NELSENSIEGEL_H_INCLUDED

#include "pnl/pnl_vector.h"
#include "read_market_data.h"

double BondPrice_NelsonSiegel(double T, double beta0, double beta1, double beta2, double tau);

void NelsonSiegelConstraints(const PnlVect *x, PnlVect *res, void *UpperLower);

double fitting_error_zc(const PnlVect *x, void *ZCMarket);

int NelsonSiegel_Fitting(ZCMarketData *ZCMarket);

#endif // FUNCTION_TO_MINIMISE_H_INCLUDED
