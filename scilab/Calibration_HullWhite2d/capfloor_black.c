#include <stdio.h>
#include <stdlib.h>

#include "pnl/pnl_cdf.h"
#include "pnl/pnl_mathtools.h"

#include "capfloor_black.h"

double atm_strike_capfloor(ZCMarketData *ZCMarket, double first_reset_date, double periodicity, double contract_maturity)
{
    double sum_zc_tau;
    double P0_opmat, P0_swapmat, Ti, SwapRate;
    int i, nb_payement;
    double P0_Ti;

    P0_opmat = BondPrice(first_reset_date, ZCMarket);
    P0_swapmat = BondPrice(contract_maturity, ZCMarket);

    nb_payement=(int)((contract_maturity-first_reset_date)/periodicity);
    Ti = first_reset_date;
    sum_zc_tau = 0.;
    for (i=0;i<nb_payement;i++)
    {
        Ti += periodicity;
        P0_Ti = BondPrice(Ti, ZCMarket);
        sum_zc_tau += periodicity*P0_Ti;
    }

    /*Compute Swap Rate*/
    SwapRate = (P0_opmat-P0_swapmat)/sum_zc_tau;

    return SwapRate;
}

/*Cap cap_floor=1, Floor_receiver=-1*/
double black_capfloor_price(ZCMarketData *ZCMarket, int cap_floor, double first_reset_date, double periodicity, double contract_maturity, double Nominal, double capfloor_strike, double vol)
{
    double sum, d1, d2, Sigma, black_forumla;
    double Tim, Tip, P0_Tim, P0_Tip, F0_Tim_Tip;
    int i, nb_payement, omega;

    omega = cap_floor;

    /*Cap Black's formula pag.17 Brigo-Mercurio*/
    nb_payement=(int)((contract_maturity-first_reset_date)/periodicity);
    Tim = first_reset_date-periodicity;
    sum = 0.;
    for (i=0; i<nb_payement; i++)
    {
        Tim += periodicity;
        Tip = Tim + periodicity;

        P0_Tim = BondPrice(Tim, ZCMarket);
        P0_Tip = BondPrice(Tip, ZCMarket);
        F0_Tim_Tip = (1./periodicity)*(P0_Tim/P0_Tip-1.);

        Sigma = vol * sqrt(Tim);

        d1 = (log(F0_Tim_Tip/capfloor_strike)/Sigma+0.5*Sigma);
        d2 = (log(F0_Tim_Tip/capfloor_strike)/Sigma-0.5*Sigma);
        black_forumla = F0_Tim_Tip*omega*cdf_nor(omega*d1)- capfloor_strike*omega*cdf_nor(omega*d2);
        sum += P0_Tip*black_forumla;
    }

    return Nominal * sum * periodicity;
}


/* Retrieve the black volatility from a cap/floor price*/
double bk_capfloor_vol_implied_bisection(ZCMarketData* ZCMarket, double capfloor_price, int cap_floor, double first_reset_date, double periodicity, double contract_maturity, double Nominal, double capfloor_strike)
{
    int i, MAX_ITERATIONS;
    double ACCURACY;
    double vol_black, vol_black_low, vol_black_high, price, diff;

    MAX_ITERATIONS = 100;
    ACCURACY    = 1.0e-8;
    vol_black_low = 0.00001;

    vol_black_high = 0.3;

    price = black_capfloor_price(ZCMarket, cap_floor, first_reset_date, periodicity, contract_maturity, Nominal, capfloor_strike, vol_black_high);

    while (price<capfloor_price)
    {
        vol_black_high *=2.;

        price = black_capfloor_price(ZCMarket, cap_floor, first_reset_date, periodicity, contract_maturity, Nominal, capfloor_strike, vol_black_high);
    }

    for (i=0; i<MAX_ITERATIONS; i++)
    {
        vol_black = 0.5*(vol_black_high+vol_black_low);

        // Prix du cap/floor
        price = black_capfloor_price(ZCMarket, cap_floor, first_reset_date, periodicity, contract_maturity, Nominal, capfloor_strike, vol_black);

        diff = capfloor_price - price;

        if (fabs(diff)<ACCURACY)
            return vol_black;

        if (diff<0) vol_black_high = vol_black;
        else vol_black_low = vol_black;
    }

    return vol_black;  // something screwy happened, should throw exception
}


double bk_capfloor_vol_implied_newton(ZCMarketData* ZCMarket, double capfloor_price, int cap_floor, double first_reset_date, double periodicity, double contract_maturity, double Nominal, double capfloor_strike)
{
    int i, MAX_ITERATIONS, nb_payement, omega;
    double ACCURACY;
    double sum_price, sum_vega, d1, d2, Sigma, black_forumla, vega, vol_black, price, diff;
    double Tim, Tip, P0_Tim, P0_Tip, F0_Tim_Tip;

    MAX_ITERATIONS = 100;
    ACCURACY    = 1.0e-8;
    vol_black = 0.15;

    omega = cap_floor;
    nb_payement=(int)((contract_maturity-first_reset_date)/periodicity);
    Tim = first_reset_date-periodicity;

    for (i=0; i<MAX_ITERATIONS; i++)
    {
        // Prix du cap/floor
        Tim = first_reset_date-periodicity;
        sum_price = 0.;
        sum_vega = 0.;
        for (i=0; i<nb_payement; i++)
        {
            Tim += periodicity;
            Tip = Tim + periodicity;

            P0_Tim = BondPrice(Tim, ZCMarket);
            P0_Tip = BondPrice(Tip, ZCMarket);
            F0_Tim_Tip = (1./periodicity)*(P0_Tim/P0_Tip-1.);

            Sigma = vol_black * sqrt(Tim);

            d1 = (log(F0_Tim_Tip/capfloor_strike)/Sigma+0.5*Sigma);
            d2 = (log(F0_Tim_Tip/capfloor_strike)/Sigma-0.5*Sigma);
            black_forumla = F0_Tim_Tip*omega*cdf_nor(omega*d1)- capfloor_strike*omega*cdf_nor(omega*d2);
            sum_price += P0_Tip * black_forumla;
            sum_vega += P0_Tip * F0_Tim_Tip * sqrt(Tim) * cdf_nor(d1);
        }

        price = Nominal * periodicity * sum_price;

        diff = capfloor_price - price;
        if (fabs(diff)<ACCURACY)
            return vol_black;

        vega = Nominal * periodicity * sum_vega;
        vol_black = vol_black + diff/vega;
    }

    return vol_black;  // something screwy happened, should throw exception
}
