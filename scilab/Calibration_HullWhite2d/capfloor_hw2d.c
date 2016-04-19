#include <stdio.h>
#include <stdlib.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_cdf.h"
#include "capfloor_hw2d.h"

static double cf_ZBOvolatility2d(double a,double sigma1,double b,double sigma2,double rho, double t, double T, double S)
{
    double sigma_p;
    double U, V, B_TS;
    double exp_atT, exp_btT, exp_aTS, exp_bTS;
    double sigma3, eta, rhoG2;

    sigma3 = sqrt(sigma1*sigma1 + sigma2*sigma2/((b-a)*(b-a)) + 2*rho*sigma1*sigma2/(b-a));
    eta = sigma2 / (a-b);
    rhoG2 = (sigma1*rho - eta)/sigma3 ;

    exp_atT = exp(-a*(T-t));
    exp_btT = exp(-b*(T-t));

    exp_aTS = exp(-a*(S-T));
    exp_bTS = exp(-b*(S-T));

    B_TS = (1 - exp_aTS) / a;
    U = (exp_aTS - 1) * exp_atT/(a*(a-b)); //(1/exp_aS - 1/exp_aT)/(a*(a-b));
    V = (exp_bTS - 1) * exp_btT/(b*(a-b)); // (1/exp_bS - 1/exp_bT)/(b*(a-b));

    sigma_p  = sigma3*sigma3*(1-exp_aTS)*(1-exp_aTS)*(1-exp_atT*exp_atT)/(2*a*a*a) ;

    sigma_p += eta*eta*(1-exp_bTS)*(1-exp_bTS)*(1-exp_btT*exp_btT)/(2*b*b*b);

    sigma_p += 2*rhoG2*sigma3*eta*(1-exp_aTS)*(1-exp_bTS)*(1-exp_atT*exp_btT)/(a*b*(a+b)) ;

    sigma_p = sqrt(sigma_p);

    return sigma_p;
}


static double cf_zbput2d(ZCMarketData* ZCMarket,double t,double r,double u,double a,double sigma1,double b,double sigma2,double rho,double S,double T, double X)
{
    double PtS,PtT;
    double h, sigma_p;
    double price;

    sigma_p = cf_ZBOvolatility2d( a, sigma1, b, sigma2, rho, t, T, S);

    PtT=BondPrice(T, ZCMarket);

    PtS=BondPrice(S, ZCMarket);

    h= log(PtS/(PtT*X)) / sigma_p + 0.5 * sigma_p ;

    price = PtS * (cdf_nor(h)-1) - X * PtT * (cdf_nor(h-sigma_p)-1);

    return price;
}

double capfloor_price_hw2d(ZCMarketData* ZCMarket, int cap_floor,double a,double sigma1,double b,double sigma2,double rho,double Nominal,double K,double periodicity,double first_payement,double contract_maturity)
{

    double sum, tim, tip, strike_put, price;
    int i, nb_payement;

    strike_put =  1./(1 + periodicity*K);
    nb_payement=(int)((contract_maturity-first_payement)/periodicity);

    /*Cap=Portfolio of zero-bond Put options*/
    sum=0.;
    for(i=0; i<nb_payement; i++)
    {
        tim   = first_payement + (double)i*periodicity;
        tip   = tim + periodicity;
        sum  += cf_zbput2d(ZCMarket, 0, 0., 0., a, sigma1, b, sigma2, rho, tip, tim, strike_put);
    }

    sum = Nominal*(1.+K*periodicity)*sum;

    /*Price*/
    price = sum;

    return price;
}
