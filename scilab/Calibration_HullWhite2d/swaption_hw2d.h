#ifndef SWAPTION_H_INCLUDED
#define SWAPTION_H_INCLUDED

#include "read_market_data.h"

double cf_swaption_hw2d(ZCMarketData* ZCMarket, int payer_receiver, double Nominal, double periodicity, double option_maturity, double contract_maturity, double swaption_strike, double a, double b, double sigma, double eta, double rho);

void HW2dparams_to_G2dparams(double a, double b, double *sigma, double *eta, double *rho);

void G2dparams_to_HW2dparams(double a, double b, double *sigma, double *eta, double *rho);

void HW2dparams_to_G2dparams_vect(PnlVect *HW2dparams, PnlVect *G2dparams);


#endif // SWAPTION_H_INCLUDED
