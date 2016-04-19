#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_integration.h"
#include "pnl/pnl_optim.h"

#include "lmm_volsto_piterbarg.h"
#include "lmm_volsto_piterbarg_calib.h"

#define MAX_PATH_LEN 1000

static int print_or_no_print;
static int n_params_vol;
static int swaption_index;
static int _n_swap;
static int _m_swap;
static int NbrSteps=100;
static double skew_mkt_mean;
static PnlMat *weight_mat;


static double func_to_intg_1(double t, void *LmmPiterbarg)
{
    return pow(SwapRate_vol((StructLmmPiterbarg*)LmmPiterbarg, t, _n_swap, _m_swap), 2);
}

static double func_to_intg_2(double s, void *LmmPiterbarg)
{
    double theta = ((StructLmmPiterbarg*)LmmPiterbarg)->Var_SpeedMeanReversion;
    return pow(SwapRate_vol((StructLmmPiterbarg*)LmmPiterbarg, s, _n_swap, _m_swap), 2)*(exp(theta*s)-exp(-theta*s))/(2*theta);
}

static double func_to_discretize(double t, StructLmmPiterbarg *LmmPiterbarg)
{
    PnlFunc func;
    int NbrPts = 20;
    double v1, v2, v_nm2;
    double theta=LmmPiterbarg->Var_SpeedMeanReversion, eta=LmmPiterbarg->Var_Volatility;

    func.params = LmmPiterbarg;
    func.function = func_to_intg_1;
    v1 = pnl_integration(&func, 0.0, t, NbrPts, "simpson");

    func.function = func_to_intg_2;
    v2 = pnl_integration(&func, 0.0, t, NbrPts, "simpson");

    v_nm2 = v1 + SQR(eta)*exp(-theta*t)*v2;

    return v_nm2*pow(SwapRate_vol(LmmPiterbarg, t, _n_swap, _m_swap), 2);
}

// uses swaption_index
static double func_to_intg_3(double t, void *LmmPiterbarg)
{

    double Tn = MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, swaption_index, 0);
    double dt = Tn/(NbrSteps-1);
    int j_t = floor(t/dt);

    return (MGET(weight_mat, swaption_index, MIN(j_t+1, weight_mat->n-1))*(t-j_t*dt) + MGET(weight_mat, swaption_index, j_t)*((j_t+1)*dt-t)) / dt;
}

// uses swaption_index
static double func_to_intg_4(double t, void *LmmPiterbarg)
{
    return func_to_intg_3(t, LmmPiterbarg)*SwapRate_skew(LmmPiterbarg, t, _n_swap, _m_swap);
}

static void SwapRate_weightmat(StructLmmPiterbarg *LmmPiterbarg)
{
    int i, j, NbrData;
    double t, Tn, Tm, dt, sumw;
    PnlFunc func;
    int NbrPts = 100;

    NbrData = LmmPiterbarg->NbrSwptData;
    weight_mat = pnl_mat_create(NbrData, NbrSteps);

    func.params = LmmPiterbarg;
    func.function = func_to_intg_3;

    for (i=0; i<NbrData; i++)
    {
        Tn = MGET(LmmPiterbarg->SwaptionGrid, i, 0);
        Tm = Tn+MGET(LmmPiterbarg->SwaptionGrid, i, 1);
        swaption_index = i;

        _n_swap = indiceTimeGrid(LmmPiterbarg->TimeDates, Tn);
        _m_swap = indiceTimeGrid(LmmPiterbarg->TimeDates, Tm);

        dt = Tn/(NbrSteps-1);
        t = 0.;
        for (j=0; j<NbrSteps; j++)
        {
            MLET(weight_mat, i, j) = func_to_discretize(t, LmmPiterbarg);
            t += dt;
        }

        sumw = pnl_integration(&func, 0.0, Tn, NbrPts, "simpson");

        for (j=0; j<NbrSteps; j++) MLET(weight_mat, i, j) /= sumw;
    }

}

static double SwapRate_skew_avg_tmp_mat(StructLmmPiterbarg *LmmPiterbarg, int i)
{
    PnlFunc func;
    int NbrPts = 100;
    double result, Tn, Tm;

    Tn = MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 0);
    Tm = Tn+MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 1);

    _n_swap = indiceTimeGrid(LmmPiterbarg->TimeDates, Tn);
    _m_swap = indiceTimeGrid(LmmPiterbarg->TimeDates, Tm);

    swaption_index = i;

    func.params = LmmPiterbarg;
    func.function = func_to_intg_4;
    result = pnl_integration(&func, 0.0, Tn, NbrPts, "simpson");

    return result;
}

void PnlVect_To_InterpObj(StructInterpObj *InterpObj, const PnlVect *x)
{
    int i, j, l, n, m;

    n = InterpObj->x->size;
    m = InterpObj->y->size;

    l=0;
    for (i=0; i<n; i++)
    {
        for (j=0; j<m; j++)
        {
            if (pnl_mat_int_get(InterpObj->zero_index,i,j)==1)
            {
                MLET(InterpObj->f_x_y, i, j) = GET(x, l);
                l++;
            }
        }
    }
}

double func_to_min_vols(const PnlVect *x, void *LmmPiterbarg)
{
    int i, j, k, n_params_1, n_params_2, n_swap, m_swap, NbrData;
    double Tn, Tm, lambda_mod, lambda_mkt, error;
    PnlVect v;
    NbrData = ((StructLmmPiterbarg*)LmmPiterbarg)->NbrSwptData;

    n_params_1 = 0;
    n_params_2 = 0;
    for (i=0; i<((StructLmmPiterbarg*)LmmPiterbarg)->NbrVolFactors; i++)
    {
        n_params_2 += pnl_mat_int_sum((((StructLmmPiterbarg*)LmmPiterbarg)->Vols_factor[i])->zero_index);
        v = pnl_vect_wrap_subvect_with_last(x, n_params_1, n_params_2-1);
        PnlVect_To_InterpObj(((StructLmmPiterbarg*)LmmPiterbarg)->Vols_factor[i], &v);
        n_params_1 = n_params_2;
    }

    error = 0.;
    for (i=0; i<NbrData; i++)
    {
        Tn = MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 0);
        Tm = Tn+MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 1);
        lambda_mkt = MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 3);

        n_swap = indiceTimeGrid(((StructLmmPiterbarg*)LmmPiterbarg)->TimeDates, Tn);
        m_swap = indiceTimeGrid(((StructLmmPiterbarg*)LmmPiterbarg)->TimeDates, Tm);

        lambda_mod = SwapRate_vol_avg(LmmPiterbarg, n_swap, m_swap, skew_mkt_mean);

        error += pow(lambda_mod-lambda_mkt, 2);
        if (print_or_no_print==1)
        printf("Tn=%f, Tm=%f, lambda_mkt = %f, lambda_mod = %f, error=%f\n", Tn, Tm, lambda_mkt, lambda_mod, fabs(lambda_mod-lambda_mkt));
    }
    if (print_or_no_print==1)
    {
        for (k=0; k<((StructLmmPiterbarg*)LmmPiterbarg)->NbrVolFactors; k++)
        {
            printf("++++++++++++++++++++++++++ Instantaneous Vols matrix[%i]++++++++++++++++++++++++++ \n", k);
            for (i=0; i<((StructLmmPiterbarg*)LmmPiterbarg)->Vols_factor[k]->f_x_y->m; i++)
            {
                for (j=0; j<((StructLmmPiterbarg*)LmmPiterbarg)->Vols_factor[k]->f_x_y->n; j++)
                {
                    if (pnl_mat_int_get(((StructLmmPiterbarg*)LmmPiterbarg)->Vols_factor[k]->zero_index, i, j)==1)
                        printf("%f ", MGET(((StructLmmPiterbarg*)LmmPiterbarg)->Vols_factor[k]->f_x_y, i, j));
                }
                printf("\n");
            }
        }
    }
    return error;
}

double func_to_min_skews(const PnlVect *x, void *LmmPiterbarg)
{
    int i, j, n_swap, m_swap, NbrData;
    double Tn, Tm, skew_mod, skew_mkt, error;
    NbrData = ((StructLmmPiterbarg*)LmmPiterbarg)->NbrSwptData;

    PnlVect_To_InterpObj(((StructLmmPiterbarg*)LmmPiterbarg)->Skews, x);

    error = 0.;
    for (i=0; i<NbrData; i++)
    {
        Tn = MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 0);
        Tm = Tn+MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 1);
        skew_mkt = MGET(((StructLmmPiterbarg*)LmmPiterbarg)->SwaptionGrid, i, 2);

        n_swap = indiceTimeGrid(((StructLmmPiterbarg*)LmmPiterbarg)->TimeDates, Tn);
        m_swap = indiceTimeGrid(((StructLmmPiterbarg*)LmmPiterbarg)->TimeDates, Tm);

        //skew_mod = SwapRate_skew_avg(LmmPiterbarg, n_swap, m_swap);
        skew_mod = SwapRate_skew_avg_tmp_mat(LmmPiterbarg, i);
        error += pow(skew_mod-skew_mkt, 2);


        if (print_or_no_print==1)
        printf("Tn=%f, Tm=%f, skew_mkt = %f, skew_mod = %f, error=%f\n", Tn, Tm, skew_mkt, skew_mod, fabs(skew_mod-skew_mkt));


    }

    if (print_or_no_print==1)
    {
        printf("++++++++++++++++++++++++++ Instantaneous Skews matrix++++++++++++++++++++++++++ \n");
        for (i=0; i<((StructLmmPiterbarg*)LmmPiterbarg)->Skews->f_x_y->m; i++)
        {
            for (j=0; j<((StructLmmPiterbarg*)LmmPiterbarg)->Skews->f_x_y->n; j++)
            {
                if (pnl_mat_int_get(((StructLmmPiterbarg*)LmmPiterbarg)->Skews->zero_index, i, j)==1)
                    printf("%f ", MGET(((StructLmmPiterbarg*)LmmPiterbarg)->Skews->f_x_y, i, j));
            }
            printf("\n");
        }
    }
    return error;
}

void LmmPiterbarg_Calib_vol(StructLmmPiterbarg *LmmPiterbarg)
{
    double tolerance, vol_mkt_mean;
    int i, CONVERGENCE, iter_max, print_algo_steps;
    int NbrVolFactors = ((StructLmmPiterbarg*)LmmPiterbarg)->NbrVolFactors;
    int NbrData = ((StructLmmPiterbarg*)LmmPiterbarg)->NbrSwptData;
    PnlVect *x_input, *x_output;
    PnlVect *lower_bounds, *upper_bounds;
    PnlRnFuncR FuncToMinimize;

    n_params_vol = 0;
    for (i=0; i<NbrVolFactors; i++)
        n_params_vol += pnl_mat_int_sum((LmmPiterbarg->Vols_factor[i])->zero_index);

    skew_mkt_mean=0.;
    for (i=0; i<NbrData; i++) skew_mkt_mean += MGET(LmmPiterbarg->SwaptionGrid, i, 2)/NbrData;
    vol_mkt_mean=0.;
    for (i=0; i<NbrData; i++) vol_mkt_mean += MGET(LmmPiterbarg->SwaptionGrid, i, 3)/NbrData;

    for (i=0; i<NbrVolFactors; i++)
        pnl_mat_set_double(LmmPiterbarg->Vols_factor[i]->f_x_y, 1./sqrt(NbrVolFactors)*vol_mkt_mean);

    print_or_no_print = 0;
    tolerance = 1e-5;
    iter_max = 150;
    print_algo_steps = 3;

    lower_bounds = pnl_vect_create_from_double(n_params_vol, 0.);
    upper_bounds = pnl_vect_create_from_double(n_params_vol, 1.);

    FuncToMinimize.function = &func_to_min_vols;
    FuncToMinimize.params = LmmPiterbarg;

    x_input = pnl_vect_create_from_double(n_params_vol, 1./sqrt(NbrVolFactors)*vol_mkt_mean);
    x_output = pnl_vect_create(n_params_vol);

    CONVERGENCE = pnl_optim_intpoints_bfgs_solve(&FuncToMinimize, NULL, NULL, lower_bounds, upper_bounds, x_input, tolerance, iter_max, print_algo_steps, x_output);

    print_or_no_print = 1;
    printf(" ++++++++++++++++++++++++++ Volatilities calibration results ++++++++++++++++++++++++++ \n");
    func_to_min_vols(x_output, LmmPiterbarg);

    pnl_mat_free(&weight_mat);
    pnl_vect_free(&x_input);
    pnl_vect_free(&x_output);
    pnl_vect_free(&lower_bounds);
    pnl_vect_free(&upper_bounds);
}

void LmmPiterbarg_Calib_skew(StructLmmPiterbarg *LmmPiterbarg)
{
    int i, NbrData = ((StructLmmPiterbarg*)LmmPiterbarg)->NbrSwptData;
    double tolerance;
    int CONVERGENCE, iter_max, print_algo_steps;
    int n_params_skew = pnl_mat_int_sum((LmmPiterbarg->Skews)->zero_index);
    PnlVect *x_input, *x_output;
    PnlVect *lower_bounds, *upper_bounds;
    PnlRnFuncR FuncToMinimize;

    skew_mkt_mean=0.;
    for (i=0; i<NbrData; i++) skew_mkt_mean += MGET(LmmPiterbarg->SwaptionGrid, i, 2)/NbrData;
    pnl_mat_set_double(LmmPiterbarg->Skews->f_x_y, skew_mkt_mean);

    SwapRate_weightmat(LmmPiterbarg);

    print_or_no_print = 0;
    tolerance = 1e-6;
    iter_max = 200;
    print_algo_steps = 3;

    lower_bounds = pnl_vect_create_from_double(n_params_skew, -1.);
    upper_bounds = pnl_vect_create_from_double(n_params_skew, 1.);

    FuncToMinimize.function = &func_to_min_skews;
    FuncToMinimize.params = LmmPiterbarg;

    x_input = pnl_vect_create_from_double(n_params_skew, skew_mkt_mean);
    x_output = pnl_vect_create_from_double(n_params_skew, skew_mkt_mean);

    CONVERGENCE = pnl_optim_intpoints_bfgs_solve(&FuncToMinimize, NULL, NULL, lower_bounds, upper_bounds, x_input, tolerance, iter_max, print_algo_steps, x_output);

    print_or_no_print = 1;
    printf(" ++++++++++++++++++++++++++ Skew calibration results ++++++++++++++++++++++++++ \n");
    func_to_min_skews(x_output, LmmPiterbarg);

    pnl_mat_free(&weight_mat);
    pnl_vect_free(&x_input);
    pnl_vect_free(&x_output);
    pnl_vect_free(&lower_bounds);
    pnl_vect_free(&upper_bounds);
}
