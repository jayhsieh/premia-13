#ifndef _LMM_PITERBARG
#define _LMM_PITERBARG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "math/read_market_zc/InitialYieldCurve.h"


// 0<T(0)<T(1)<....<T(N)
// L(0), L(1), ..., L(N-1)
// beta(0), beta(1), ..., beta(N-1)
// sigma(0), sigma(1), ..., sigma(N-1)
typedef struct StructLmmPiterbarg
{
    ZCMarketData *ZCMarket; // Structure that constraints initial yield curve

    PnlVect *TimeDates;

    double Var_SpeedMeanReversion;
    double Var_Volatility;
    int NbrVolFactors;

    PnlMat* SkewsParams;
    PnlMat* VolsParams;

}StructLmmPiterbarg;


StructLmmPiterbarg* SetLmmPiterbarg(int InitYieldCurve_flag, double R_flat, char *curve, double period, double T_last, double Var_SpeedMeanReversion, double Var_Volatility, int NbrVolFactors, PnlMat* SkewsParams, PnlMat* VolsParams);

double LiborRate_skew(StructLmmPiterbarg *LmmPiterbarg, double t, double Tn);

void FreeLmmPiterbarg(StructLmmPiterbarg **LmmPiterbarg);

int indiceTimeGrid(PnlVect *TimeGrid, double s);
double SwapRate_skew(StructLmmPiterbarg *LmmPiterbarg, double t, int n_swap, int m_swap);
double SwapRate_vol(StructLmmPiterbarg *LmmPiterbarg, double t, int n_swap, int m_swap);

double SwapRate_skew_avg(StructLmmPiterbarg *LmmPiterbarg, int n_swap, int m_swap);

double SwapRate_vol_avg(StructLmmPiterbarg *LmmPiterbarg, int n_swap, int m_swap, double skew_n_m);

#endif
