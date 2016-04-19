
#ifndef QUADRATICMODEL_H
#define QUADRATICMODEL_H

#include "math/read_market_zc/InitialYieldCurve.h"

/**
The instantaneous sport interest rate r_t is described by r_t = 0.5*x_t^2, with x_t following SDE :

dx_t = (alpha_t - beta * x_t) * dt + sigma * dW_t
x_0 = sqrt(2*r_0)
**/

// Strucure which contains information on the the T-maturity zc bond at time t=0 under, Quadratic model.
// P(0, T) = exp(- ( B*r + b * sqrt(2*r) + c))
typedef struct
{
    double T;
    double P; // Price of the T-maturity bond at time t=0
    double f0_T; // T-maturity forward rate at time t=0
    double B; // Coefficients of the T-maturity bond at time t=0 : P(0,T)=exp(-(.5*B*x^2+b*x+c))
    double b;
    double c;
    double dB; // Derivatives of B and b with respect to T
    double db;
    double V; // Variance of x under the T-forward probability
} Data;


// Coefficents of the omega distribution : .5*B*x^2+bx+c where x is normally distributed with mean mu and variance V
typedef struct
{
  double B;
  double b;
  double c;
  double mu;
  double V;
} Omega;



// coefficients of the chi^2 distribution : alpha+beta X where X is non centrally chi^2 distributed with nu degree of freedom and non-centrality parameter lambda
typedef struct
{
  double nu;
  double lambda;
  double beta;
  double alpha;
} Chn;

// Computes the strcuture data at time T
void bond_coeffs(ZCMarketData* ZCMarket, Data *data, double T,  double beta, double sigma, double x0);

// Gives the omega distribution of the zero-coupon bond P(T, S) data1 contains the coefficients of bond P(0,T), data2 contains the coefficients of bond P(0,S).
void transport(Omega *om, Data data1, Data data2, double beta, double sigma, double x0);

// Transforms Omega distribution to a chi^2 distribution
void om2chn(Omega om, Chn *chn);

// Compute the initial rate r_0 and corresponding value x_0
void initial_short_rate(ZCMarketData* ZCMarket, double *r0, double *x0);


/* Price of a european option on zero coupon bond*/
double zb_call_quad1d(ZCMarketData *ZCMarket, double beta, double sigma, double T, double S, double strike);
double zb_put_quad1d(ZCMarketData *ZCMarket, double beta, double sigma, double T, double S, double strike);


#endif
