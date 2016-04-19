
#ifndef HULLWHITE1DINCLUDES_H
#define HULLWHITE1DINCLUDES_H

#include<stdlib.h>
#include "hullwhite1d_stdi.h"
#include "math/read_market_zc/InitialYieldCurve.h"

/// Computation of the two coefficient A(t,T) and B(t,T) used in the price of a Zero Coupon.
/// These two coefficients are independants of r and u.
void ZCPrice_CoefficientHW1D(ZCMarketData* ZCMarket, double a, double sigma, double t, double T, double* A_tT, double* B_tT);

/// Price of a ZC using the two coefficient A(t,T) and B(t,T) . H&W is a affine model.
double ZCPrice_Using_CoefficientHW1D(double r_t, double A_tT, double B_tT);

/// Price at date t of a ZC maturing at T, knowing that r(t)=r_t
double cf_hw1d_zcb(ZCMarketData* ZCMarket, double a, double sigma, double t, double r_t, double T);

/// Price at date t of a european put on a ZC maturing at T, knowing that r(t)=r_t.
double cf_hw1d_zbput(ZCMarketData* ZCMarket, double a,double sigma, double S, double T, double K);

/// Price at date t of a european call on a ZC maturing at T, knowing that r(t)=r_t.
double cf_hw1d_zbcall(ZCMarketData* ZCMarket, double a,double sigma, double S, double T, double K);

#endif
