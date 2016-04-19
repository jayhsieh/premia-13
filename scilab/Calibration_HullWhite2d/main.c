#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_random.h"

#include "read_market_data.h"
#include "swaption_black.h"
#include "swaption_hw2d.h"
#include "capfloor_black.h"
#include "capfloor_hw2d.h"
#include "function_to_minimise.h"

#include "pnl/pnl_optim.h"

// Constraints on variables
void Constraints_HW2d(const PnlVect *x, PnlVect *res, void *p)
{
    pnl_vect_resize(res, 1);
    LET(res, 0) = SQR(GET(x, 0)-GET(x, 2))-1e-8; // To prevent a==b
}

int calibration_cap(int type_generator, int print_algo_steps)
{
    double tolerance, a, sigma1, b, sigma2, rho, temps_cpu;
    int iter_max, CONVERGENCE;

    clock_t temps_initial, temps_final;
    PnlVect *ParamsOutput_HW2d, *ParamsOutput_G2pp, *ParamsInit_HW2d;
    PnlVect *lower_bounds, *upper_bounds;

    ZCMarketData ZCMarket;
    MktATMCapVolData MktATMCapVol;
    MktCapZCData MktCapZC;
    PnlRnFuncRm Constraints; // Structure that contains the constraints

    PnlRnFuncR FuncToMinimize;

    ZCMarket.FlatOrMarket = 1;
    ReadMarketData(&ZCMarket);
    ReadCapMarketData(&MktATMCapVol);
    MktCapZC.cap_vol_market = &MktATMCapVol;
    MktCapZC.zc_market=&ZCMarket;

    ParamsInit_HW2d = pnl_vect_create(5);
    ParamsOutput_HW2d = pnl_vect_create(5);
    ParamsOutput_G2pp = pnl_vect_create(5);

    lower_bounds = pnl_vect_create_from_list(5, -3.0, 1e-8, -3.0, 1e-8, -0.9999);
    upper_bounds = pnl_vect_create_from_list(5, 3.0, 1.0, 3.0, 1.0, 0.9999);

    printf("x = [a, sigma1, b, sigma2, rho]\n");
    printf("bound_min : ");
    pnl_vect_print(lower_bounds);
    printf("bound_max : ");
    pnl_vect_print(upper_bounds);

    ChooseInputParameters(ParamsInit_HW2d, lower_bounds, upper_bounds, type_generator);

    FuncToMinimize.function = &fitting_error_cap;
    FuncToMinimize.params = &MktCapZC;

    Constraints.function = Constraints_HW2d;
    Constraints.params = NULL;

    tolerance = 1e-10;
    iter_max = 2000;

    temps_initial = clock();
    CONVERGENCE = pnl_optim_intpoints_bfgs_solve(&FuncToMinimize, NULL, &Constraints, lower_bounds, upper_bounds, ParamsInit_HW2d, tolerance, iter_max, print_algo_steps, ParamsOutput_HW2d);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;

    if (CONVERGENCE==0) abort();

    a = GET(ParamsOutput_HW2d, 0);
    sigma1 = GET(ParamsOutput_HW2d, 1);
    b = GET(ParamsOutput_HW2d, 2);
    sigma2 = GET(ParamsOutput_HW2d, 3);
    rho = GET(ParamsOutput_HW2d, 4);
    HW2dparams_to_G2dparams(a, b, &sigma1, &sigma2, &rho);
    LET(ParamsOutput_G2pp, 0) = a;
    LET(ParamsOutput_G2pp, 1) = sigma1;
    LET(ParamsOutput_G2pp, 2) = b;
    LET(ParamsOutput_G2pp, 3) = sigma2;
    LET(ParamsOutput_G2pp, 4) = rho;

    printf("--------------------------- Calibration results --------------------------- \n");

    printf("****** Time elapsed = %i mn %i s \n\n", ((int)temps_cpu) / 60,  ((int)temps_cpu) % 60);

    printf("****** Objective_function(Input  point)= %.10f \n", PNL_EVAL_RNFUNCR(&FuncToMinimize, ParamsInit_HW2d));
    printf("****** Objective_function(Output point)= %.10f \n", PNL_EVAL_RNFUNCR(&FuncToMinimize, ParamsOutput_HW2d));

    printf("****** Output-Parameters of the model Hull&White 2-factor: \n");
    pnl_vect_print(ParamsOutput_HW2d);

    printf("****** Equivalent parameters of the model G2++: \n");
    pnl_vect_print(ParamsOutput_G2pp);

    printf("\n****** Model_cap_vol : \n");
    print_model_cap_vol(ParamsOutput_HW2d, &MktCapZC);
    printf("\n*****************************\n");

    pnl_vect_free(&ParamsOutput_HW2d);
    pnl_vect_free(&ParamsOutput_G2pp);
    pnl_vect_free(&ParamsInit_HW2d);

    pnl_vect_free(&lower_bounds);
    pnl_vect_free(&upper_bounds);

    DeleteZCMarketData(&ZCMarket);
    DeleteMktATMCapVolData(&MktATMCapVol);

    return 0;
}


int calibration_swaption(int type_generator, int print_algo_steps)
{
    double tolerance, a, sigma1, b, sigma2, rho, temps_cpu;
    int iter_max, CONVERGENCE;

    clock_t temps_initial, temps_final;
    PnlVect *ParamsOutput_HW2d, *ParamsOutput_G2pp, *ParamsInit_HW2d;
    PnlVect *lower_bounds, *upper_bounds;

    ZCMarketData ZCMarket;
    MktATMSwaptionVolData MktATMSwaptionVol;
    MktSwpZCData MktSwpZC;
    PnlRnFuncR FuncToMinimize;
    PnlRnFuncRm Constraints; // Structure that contains the constraints

    ZCMarket.FlatOrMarket = 1;
    ReadMarketData(&ZCMarket);
    ReadSwaptionMarketData(&MktATMSwaptionVol);
    MktSwpZC.swaption_vol_market = &MktATMSwaptionVol;
    MktSwpZC.zc_market=&ZCMarket;

    ParamsInit_HW2d = pnl_vect_create(5);
    ParamsOutput_HW2d = pnl_vect_create(5);
    ParamsOutput_G2pp = pnl_vect_create(5);

    lower_bounds = pnl_vect_create_from_list(5, -3.0, 1e-8, -3.0, 1e-8, -0.9999);
    upper_bounds = pnl_vect_create_from_list(5, 3.0, 1.0, 3.0, 1.0, 0.9999);

    printf("x = [a, sigma1, b, sigma2, rho]\n");
    printf("bound_min : ");
    pnl_vect_print(lower_bounds);
    printf("bound_max : ");
    pnl_vect_print(upper_bounds);

    ChooseInputParameters(ParamsInit_HW2d, lower_bounds, upper_bounds, type_generator);

    FuncToMinimize.function = &fitting_error_swaption;
    FuncToMinimize.params = &MktSwpZC;

    Constraints.function = Constraints_HW2d;
    Constraints.params = NULL;

    tolerance = 1e-10;
    iter_max = 150;

    temps_initial = clock();
    CONVERGENCE = pnl_optim_intpoints_bfgs_solve(&FuncToMinimize, NULL, &Constraints, lower_bounds, upper_bounds, ParamsInit_HW2d, tolerance, iter_max, print_algo_steps, ParamsOutput_HW2d);
    temps_final = clock();
    temps_cpu = ((double) (temps_final - temps_initial)) / CLOCKS_PER_SEC;

    if (CONVERGENCE==0) abort();

    a = GET(ParamsOutput_HW2d, 0);
    sigma1 = GET(ParamsOutput_HW2d, 1);
    b = GET(ParamsOutput_HW2d, 2);
    sigma2 = GET(ParamsOutput_HW2d, 3);
    rho = GET(ParamsOutput_HW2d, 4);
    HW2dparams_to_G2dparams(a, b, &sigma1, &sigma2, &rho);
    LET(ParamsOutput_G2pp, 0) = a;
    LET(ParamsOutput_G2pp, 1) = sigma1;
    LET(ParamsOutput_G2pp, 2) = b;
    LET(ParamsOutput_G2pp, 3) = sigma2;
    LET(ParamsOutput_G2pp, 4) = rho;

    printf("--------------------------- Calibration results --------------------------- \n");

    printf("****** Time elapsed = %i mn %i s \n\n", ((int)temps_cpu) / 60,  ((int)temps_cpu) % 60);

    printf("****** Objective_function(Input  point)= %.10f \n", PNL_EVAL_RNFUNCR(&FuncToMinimize, ParamsInit_HW2d));
    printf("****** Objective_function(Output point)= %.10f \n", PNL_EVAL_RNFUNCR(&FuncToMinimize, ParamsOutput_HW2d));

    printf("****** Output-Parameters of the model Hull&White 2-factor: \n");
    pnl_vect_print(ParamsOutput_HW2d);

    printf("****** Equivalent parameters of the model G2++: \n");
    pnl_vect_print(ParamsOutput_G2pp);

    printf("****** Model_swaption_vol : \n");
    print_model_swaption_vol(ParamsOutput_HW2d, &MktSwpZC);

    pnl_vect_free(&ParamsOutput_HW2d);
    pnl_vect_free(&ParamsOutput_G2pp);
    pnl_vect_free(&ParamsInit_HW2d);

    pnl_vect_free(&lower_bounds);
    pnl_vect_free(&upper_bounds);

    DeleteZCMarketData(&ZCMarket);
    DeleteMktATMSwaptionVolData(&MktATMSwaptionVol);

    return 0;
}

int main()
{
    int print_algo_steps, type_generator, calib_cap_or_swaption, exit_or_continue;

    calib_cap_or_swaption=0;
    exit_or_continue=0;
    print_algo_steps=0;

    type_generator = PNL_RNG_MERSENNE_RANDOM_SEED;
    pnl_rand_init(type_generator, 1, 5);

Execution:

    printf("\nCalibration over cap (=0) or swaption (=1) ? ");
    scanf("%d", &calib_cap_or_swaption);

    if (calib_cap_or_swaption==1)
    {
        printf("\n********** Calibration of Hull&White 2-factors model over swaption volatilities\n");
        printf("\nDon't print optimization algorithm steps (=0) or print (=1) ? ");
        scanf("%d", &print_algo_steps);
        print_algo_steps *=2;

        calibration_swaption(type_generator, print_algo_steps);
    }

    else
    {
        printf("\n********** Calibration of Hull&White 2-factors model over cap volatilities\n");
        printf("\nDon't print optimization algorithm steps (=0) or print (=1) ? ");
        scanf("%d", &print_algo_steps);
        print_algo_steps *=2;

        calibration_cap(type_generator, print_algo_steps);
    }

    printf("\nBack to calibration (=0) or Exit (=1) ? ");
    scanf("%d", &exit_or_continue);

    if (exit_or_continue==0) goto Execution;

    else return 0;
}
