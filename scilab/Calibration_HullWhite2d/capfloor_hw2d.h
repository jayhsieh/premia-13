#ifndef CAPFLOOR_HW2D_H_INCLUDED
#define CAPFLOOR_HW2D_H_INCLUDED

#include "read_market_data.h"

double cf_ZBOvolatility_hw2d(double a, double sigma1, double b, double sigma2, double rho, double T, double S);

double cf_zcbond_option_hw2d(ZCMarketData* ZCMarket, int call_put, double a, double sigma1, double b, double sigma2, double rho, double S, double T, double X);

double capfloor_price_hw2d(ZCMarketData* ZCMarket, int cap_floor,double a,double sigma1,double b,double sigma2,double rho,double Nominal,double K,double periodicity,double first_payement,double contract_maturity);

#endif // CAPFLOOR_HW2D_H_INCLUDED
