#ifndef CAPFLOOR_BLACK_H_INCLUDED
#define CAPFLOOR_BLACK_H_INCLUDED

#include "read_market_data.h"

double atm_strike_capfloor(ZCMarketData *ZCMarket, double option_mat, double periodicity, double swap_mat);

double black_capfloor_price(ZCMarketData *ZCMarket, int cap_floor, double first_reset_date, double periodicity, double contract_maturity, double Nominal, double capfloor_strike, double vol);

double bk_capfloor_vol_implied_newton(ZCMarketData* ZCMarket, double capfloor_price, int cap_floor, double first_reset_date, double periodicity, double contract_maturity, double Nominal, double capfloor_strike);

double bk_capfloor_vol_implied_bisection(ZCMarketData* ZCMarket, double capfloor_price, int cap_floor, double first_reset_date, double periodicity, double contract_maturity, double Nominal, double capfloor_strike);

#endif // CAPFLOOR_HW2D_H_INCLUDED
