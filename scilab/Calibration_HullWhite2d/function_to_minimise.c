#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_random.h"

#include "read_market_data.h"
#include "capfloor_black.h"
#include "capfloor_hw2d.h"
#include "swaption_black.h"
#include "swaption_hw2d.h"
#include "function_to_minimise.h"

///************************ Cap Fitting ********************************///
// Quadratic error function for market cap prices
double fitting_error_cap(const PnlVect *x, void *MktATMCapVol)
{
    int i, NbrOfData;
    double Nominal, periodicity, option_maturity, contract_maturity, contract_strike;
    double a, sigma1, b, sigma2, rho;
    double mkt_option_price, mkt_cap_black_vol, mod_option_price, error_fitting, factor;

    a = GET(x, 0);
    sigma1 = GET(x, 1);
    b = GET(x, 2);
    sigma2 = GET(x, 3);
    rho = GET(x, 4);

    Nominal = 100.;
    NbrOfData = ((MktCapZCData*)MktATMCapVol)->cap_vol_market->NbrData;

    error_fitting = 0;
    factor = 1e6;
    for (i=0; i<NbrOfData; i++)
    {
        periodicity = GET(((MktCapZCData*)MktATMCapVol)->cap_vol_market->CapPeriodicity, i);
        option_maturity = GET(((MktCapZCData*)MktATMCapVol)->cap_vol_market->CapFirstResetDate, i);

        contract_maturity = GET(((MktCapZCData*)MktATMCapVol)->cap_vol_market->CapMaturity, i);
        contract_strike = atm_strike_capfloor(((MktCapZCData*)MktATMCapVol)->zc_market, option_maturity, periodicity, contract_maturity);

        mkt_cap_black_vol = GET(((MktCapZCData*)MktATMCapVol)->cap_vol_market->CapVolatility, i);
        mkt_option_price = black_capfloor_price(((MktCapZCData*)MktATMCapVol)->zc_market, 1, option_maturity, periodicity, contract_maturity, Nominal, contract_strike, mkt_cap_black_vol);

        mod_option_price = capfloor_price_hw2d(((MktCapZCData*)MktATMCapVol)->zc_market, 1, a, sigma1, b, sigma2, rho, Nominal, contract_strike, periodicity, option_maturity, contract_maturity);

        error_fitting += factor*SQR(mod_option_price-mkt_option_price)/SQR(mkt_option_price);
    }

    return error_fitting;
}

// print the result of calibration form the parameter x
void  print_model_cap_vol(const PnlVect *x, MktCapZCData *MktATMCapVol)
{
    int i, NbrOfData;
    double Nominal, periodicity, option_maturity, contract_maturity, contract_strike;
    double a, sigma1, b, sigma2, rho;
    double mod_option_price, mkt_cap_black_vol, mod_cap_black_vol, error_fitting, percentage_diff;

    a = GET(x, 0);
    sigma1 = GET(x, 1);
    b = GET(x, 2);
    sigma2 = GET(x, 3);
    rho = GET(x, 4);

    Nominal = 1.;
    NbrOfData = (MktATMCapVol->cap_vol_market)->NbrData;

    error_fitting = 0;
    for (i=0; i<NbrOfData; i++)
    {
        periodicity = GET((MktATMCapVol->cap_vol_market)->CapPeriodicity, i);
        option_maturity = GET((MktATMCapVol->cap_vol_market)->CapFirstResetDate, i);

        contract_maturity = GET((MktATMCapVol->cap_vol_market)->CapMaturity, i);
        mkt_cap_black_vol = GET((MktATMCapVol->cap_vol_market)->CapVolatility, i);

        contract_strike = atm_strike_capfloor(MktATMCapVol->zc_market, option_maturity, periodicity, contract_maturity);

        mod_option_price = capfloor_price_hw2d(MktATMCapVol->zc_market, 1, a, sigma1, b, sigma2, rho, Nominal, contract_strike, periodicity, option_maturity, contract_maturity);

        mod_cap_black_vol = bk_capfloor_vol_implied_newton(MktATMCapVol->zc_market, mod_option_price, 1, option_maturity, periodicity, contract_maturity, Nominal, contract_strike);

        percentage_diff = 100*(mod_cap_black_vol-mkt_cap_black_vol)/mkt_cap_black_vol;

        printf("T0 = %1.2f, Tn = %2.0f, mkt_vol = %.4f mod_vol= %.4f percentage_diff = %.2f %%\n",  option_maturity, contract_maturity, mkt_cap_black_vol, mod_cap_black_vol, percentage_diff);

    }

}

///********************* Swaption Fitting *****************************///
// Quadratic error function for market swaption volatilities
double fitting_error_swaption(const PnlVect *x, void *MktATMSwaptionVol)
{
    int i, NbrOfData;
    double Nominal, periodicity, option_maturity, contract_maturity, contract_strike;
    double a, sigma1, b, sigma2, rho;
    double mkt_option_price, mkt_swaption_black_vol, mod_option_price, error_fitting, factor;

    a = GET(x, 0);
    sigma1 = GET(x, 1);
    b = GET(x, 2);
    sigma2 = GET(x, 3);
    rho = GET(x, 4);

    Nominal = 100.;
    periodicity = (((MktSwpZCData*)MktATMSwaptionVol)->swaption_vol_market)->Periodicity;
    NbrOfData   = (((MktSwpZCData*)MktATMSwaptionVol)->swaption_vol_market)->NbrData;

    error_fitting = 0;
    factor = 1e2;
    for (i=0; i<NbrOfData; i++)
    {
        option_maturity   = GET((((MktSwpZCData*)MktATMSwaptionVol)->swaption_vol_market)->SwaptionMaturity, i);
        contract_maturity = GET((((MktSwpZCData*)MktATMSwaptionVol)->swaption_vol_market)->SwaptionTenor, i);

        contract_strike = atm_strike_swaption(((MktSwpZCData*)MktATMSwaptionVol)->zc_market, option_maturity, periodicity, contract_maturity);

        mkt_swaption_black_vol = GET((((MktSwpZCData*)MktATMSwaptionVol)->swaption_vol_market)->SwaptionVolatility, i);
        mkt_option_price = black_swaption_price(((MktSwpZCData*)MktATMSwaptionVol)->zc_market, 1, option_maturity, periodicity, contract_maturity, Nominal, contract_strike, mkt_swaption_black_vol);

        mod_option_price = cf_swaption_hw2d(((MktSwpZCData*)MktATMSwaptionVol)->zc_market, 1, Nominal, periodicity, option_maturity, contract_maturity, contract_strike, a, b, sigma1, sigma2, rho);

        error_fitting += factor * SQR(mod_option_price-mkt_option_price)/SQR(mkt_option_price);
    }

    return error_fitting;
}

// print the result of calibration form the parameter x
void  print_model_swaption_vol(const PnlVect *x, MktSwpZCData *MktATMSwaptionVol)
{
    int i, NbrOfData;
    double Nominal, periodicity, option_maturity, contract_maturity, contract_strike;
    double a, sigma1, b, sigma2, rho;
    double mkt_option_price, mkt_swaption_black_vol, mod_option_price, mod_swaption_black_vol, percentage_diff;

    a = GET(x, 0);
    sigma1 = GET(x, 1);
    b = GET(x, 2);
    sigma2 = GET(x, 3);
    rho = GET(x, 4);

    Nominal = 1.;
    periodicity = (MktATMSwaptionVol->swaption_vol_market)->Periodicity;
    NbrOfData = (MktATMSwaptionVol->swaption_vol_market)->NbrData;

    printf("\nSwaption calibration results:\n");

    for (i=0; i<NbrOfData; i++)
    {
        option_maturity = GET((MktATMSwaptionVol->swaption_vol_market)->SwaptionMaturity, i);
        contract_maturity = GET((MktATMSwaptionVol->swaption_vol_market)->SwaptionTenor, i);

        contract_strike = atm_strike_swaption(MktATMSwaptionVol->zc_market, option_maturity, periodicity, contract_maturity);

        mkt_swaption_black_vol = GET((MktATMSwaptionVol->swaption_vol_market)->SwaptionVolatility, i);
        mkt_option_price = black_swaption_price(MktATMSwaptionVol->zc_market, 1, option_maturity, periodicity, contract_maturity, Nominal, contract_strike, mkt_swaption_black_vol);

        mod_option_price = cf_swaption_hw2d(MktATMSwaptionVol->zc_market, 1, Nominal, periodicity, option_maturity, contract_maturity, contract_strike, a, b, sigma1, sigma2, rho);

        mod_swaption_black_vol = bk_swaption_vol_implied_newton(MktATMSwaptionVol->zc_market, mod_option_price, 1, option_maturity, periodicity, contract_maturity, Nominal, contract_strike);

        percentage_diff = 100*(mod_swaption_black_vol-mkt_swaption_black_vol)/mkt_swaption_black_vol;

        printf("T0=%1.2f, Tn=%2.0f, mkt_vol = %.4f mod_vol= %.4f percentage_diff = %.2f %%\n",  option_maturity, contract_maturity-option_maturity, mkt_swaption_black_vol, mod_swaption_black_vol, percentage_diff);

    }

}

///************** Choose initial parameters to start the optimization routine with ****************///
void ChooseInputParameters(PnlVect* ParamsInit_HW2d, PnlVect*lower_bounds, PnlVect* upper_bounds, int type_generator)
{
    int i, bad_choice, random_point_or_choose;
    double x;
    char parameters_names[][8] = {"a", "sigma1", "b", "sigma2", "rho", };

    random_point_or_choose=0;
    x=0.;

    pnl_vect_resize(ParamsInit_HW2d, 5);

    printf("\nUse random initial point (=0) Or choose (=1) ? ");
    scanf("%i", &random_point_or_choose);

    if (random_point_or_choose==1)
    {
        for (i=0; i<5; i++)
        {
            do
            {
                printf("Enter the value of the parameter '%s':", parameters_names[i]);
                scanf("%lf", &x);
                LET(ParamsInit_HW2d, i) = x;
                bad_choice = (GET(lower_bounds, i)>=x ||GET(upper_bounds, i)<=x);
                if (bad_choice)
                    printf("Initial point should be between lower_bound and upper_bound! Please Choose \n");
            }
            while (bad_choice);
        }
    }

    else
    {
        ChooseRandomInputParameters(ParamsInit_HW2d, lower_bounds, upper_bounds, type_generator);
    }
}

void ChooseRandomInputParameters(PnlVect* x_input, PnlVect*lower_bounds, PnlVect* upper_bounds, int type_generator)
{
    int i, N;
    double theta, x_min, x_max;

    N = x_input->size;

    for (i=0; i<N; i++)
    {
        theta = 0.25*pnl_rand_uni(type_generator)+0.5;
        x_min = GET(lower_bounds, i);
        x_max = GET(upper_bounds, i);

        LET(x_input, i) = theta*x_min + (1-theta)*x_max;
    }
}
