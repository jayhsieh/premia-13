#ifndef SWAPTION_BLACK_H_INCLUDED
#define SWAPTION_BLACK_H_INCLUDED

#include "InitialYieldCurve.h"

double black_swaption_price(ZCMarketData *ZCMarket, int payer_receiver, double option_mat, double periodicity, double swap_mat, double Nominal, double swaption_strike, double vol);

double bk_swaption_vol_implied_newton(ZCMarketData* ZCMarket, double swaption_price, int payer_receiver, double option_mat, double periodicity, double swap_mat, double Nominal, double swaption_strike);

#endif // SWAPTION_BLACK_H_INCLUDED
