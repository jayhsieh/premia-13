#ifndef _LMM_VOLSTO_PITERBARG
#define _LMM_VOLSTO_PITERBARG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "InitialYieldCurve.h"


typedef struct StructInterpObj
{
    PnlVect *x;
    PnlVect *y;
    PnlMat *f_x_y; // f[i,j] = f(x[i], y[j])
    PnlMatInt *zero_index;

}StructInterpObj;

// 0<T(0)<T(1)<....<T(N)
// L(0), L(1), ..., L(N-1)
// beta(0), beta(1), ..., beta(N-1)
// sigma(0), sigma(1), ..., sigma(N-1)
typedef struct StructLmmPiterbarg
{
    int NbrSwptData;
    PnlMat *SwaptionGrid;
    ZCMarketData *ZCMarket; // Structure that constraints initial yield curve

    PnlVect *TimeDates;

    double Var_SpeedMeanReversion;
    double Var_Volatility;
    int NbrVolFactors;

    StructInterpObj *Skews;
    StructInterpObj **Vols_factor;

}StructLmmPiterbarg;


StructLmmPiterbarg* SetLmmPiterbarg(char *FileSwptGridData, int InitYieldCurve_flag, double R_flat, double period, double T_last, double Var_SpeedMeanReversion, double Var_Volatility, int NbrVolFactors);

double LiborRate_skew(StructLmmPiterbarg *LmmPiterbarg, double t, double Tn);

void FreeLmmPiterbarg(StructLmmPiterbarg **LmmPiterbarg);

int indiceTimeGrid(PnlVect *TimeGrid, double s);
double SwapRate_skew(StructLmmPiterbarg *LmmPiterbarg, double t, int n_swap, int m_swap);
double SwapRate_vol(StructLmmPiterbarg *LmmPiterbarg, double t, int n_swap, int m_swap);

double SwapRate_skew_avg(StructLmmPiterbarg *LmmPiterbarg, int n_swap, int m_swap);

double SwapRate_vol_avg(StructLmmPiterbarg *LmmPiterbarg, int n_swap, int m_swap, double skew_n_m);

#endif
