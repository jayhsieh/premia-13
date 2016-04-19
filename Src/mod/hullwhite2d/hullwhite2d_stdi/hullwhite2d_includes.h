
#ifndef HULLWHITE2DINCLUDES_H
#define HULLWHITE2DINCLUDES_H

#include<stdlib.h>
#include "hullwhite2d_stdi.h"
#include "math/read_market_zc/InitialYieldCurve.h"

/// Computation of the three coefficient A(t,T), B(t,T) and C(tT) used in the price of a Zero Coupon.
/// These three coefficients are independants of r and u.
void ZCPrice_Coefficient(ZCMarketData* ZCMarket, double a, double sigma1, double b, double sigma2, double rho, double t,
                                    double T, double* A_tT, double* B_tT, double* C_tT);

/// Price of a ZC using the three coefficient A(t,T), B(t,T) and C(tT). H&W is a affine model.
double ZCPrice_Using_Coefficient(double r_t, double u_t, double A_tT, double B_tT, double C_tT);

/// Price at date t of a ZC maturing at T, knowing that r(t)=r_t and u(t)=u_t.
double cf_hw2d_zcb(ZCMarketData* ZCMarket, double a, double sigma1, double b, double sigma2, double rho, double t, double r_t, double u_t, double T);

#endif
