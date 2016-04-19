#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_optim.h"
#include "read_market_data.h"
#include "function_to_minimise.h"
#include "nelsen_siegel.h"

double BondPrice_NelsonSiegel(double T, double beta0, double beta1, double beta2, double tau)
{
    double exp_T, R;

    if(T<1e-8)
    {
        return 1.0;
    }

    exp_T = exp(-T/tau);
    R = beta0 + beta1*(1-exp_T)*tau/T + beta2*((1-exp_T)*tau/T - exp_T);

    return exp(-0.01*R*T);
}

void NelsonSiegelConstraints(const PnlVect *x, PnlVect *res, void *pt)
{
    pnl_vect_resize(res, 1);
    LET(res, 0) = GET(x, 0)+GET(x, 1);
}

double fitting_error_zc(const PnlVect *x, void *ZCMarket)
{
    int i, NbrOfData;
    double beta0, beta1, beta2, tau;

    double zc_maturity, mkt_price, mod_price, error_fitting, factor;

    beta0 = GET(x, 0);
    beta1 = GET(x, 1);
    beta2 = GET(x, 2);
    tau = GET(x, 3);

    NbrOfData = ((ZCMarketData*)ZCMarket)->Nvalue;

    error_fitting = 0;
    factor = 1e4;
    for (i=0; i<NbrOfData; i++)
    {

        zc_maturity = GET(((ZCMarketData*)ZCMarket)->tm, i);

        mkt_price = GET(((ZCMarketData*)ZCMarket)->Pm, i);

        mod_price = BondPrice_NelsonSiegel(zc_maturity, beta0, beta1, beta2, tau);

        error_fitting += factor*SQR(mod_price-mkt_price)/SQR(mkt_price);
    }

    return error_fitting;
}

int NelsonSiegel_Fitting(ZCMarketData *ZCMarket)
{
    double tolerance;
    int iter_max, CONVERGENCE, print_algo_steps;

    PnlVect *ParamsInit;
    PnlVect *lower_bounds, *upper_bounds;

    PnlRnFuncR FuncToMinimize;
    PnlRnFuncRm NL_Constraints;

    print_algo_steps=0;

    ParamsInit = pnl_vect_create_from_double(4, 1.);
    pnl_vect_resize(ZCMarket->NeslenSiegelParams, 4);

    lower_bounds = pnl_vect_create_from_list(4, 1e-8, -30.0, -30.0, 1e-8);
    upper_bounds = pnl_vect_create_from_list(4, 30.0, 30.0, 30.0, 30.0);

    FuncToMinimize.function = &fitting_error_zc;
    FuncToMinimize.params = ZCMarket;

    NL_Constraints.function = &NelsonSiegelConstraints;
    NL_Constraints.params = NULL;

    tolerance = 1e-10;
    iter_max = 1000;

    CONVERGENCE = pnl_optim_intpoints_bfgs_solve(&FuncToMinimize, NULL, &NL_Constraints, lower_bounds, upper_bounds, ParamsInit, tolerance, iter_max, print_algo_steps, ZCMarket->NeslenSiegelParams);

    if (CONVERGENCE==0) abort();

    pnl_vect_free(&ParamsInit);
    pnl_vect_free(&lower_bounds);
    pnl_vect_free(&upper_bounds);

    return 0;
}
