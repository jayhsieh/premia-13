#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "InitialYieldCurve.h"
#include "lmm_volsto_piterbarg.h"
#include "lmm_volsto_piterbarg_calib.h"
#include "lmm_volsto_piterbarg_pricing.h"

int main()
{
    double period, T_last, Var_SpeedMeanReversion, Var_Volatility, R_flat;
    int NbrVolFactors, InitYieldCurve_flag;
    char FileSwptGridData[]="market_skews.dat";

    StructLmmPiterbarg *LmmPiterbarg;

    InitYieldCurve_flag = 0;
    R_flat = 0.05;

    NbrVolFactors = 1;
    Var_SpeedMeanReversion = 0.15;
    Var_Volatility = 1.3;

    period = 0.5;
    T_last = 40.;

    LmmPiterbarg = SetLmmPiterbarg(FileSwptGridData, InitYieldCurve_flag, R_flat, period, T_last, Var_SpeedMeanReversion, Var_Volatility, NbrVolFactors);


    printf("--------------------------------- Model Calibration --------------------------------- \n");
    printf("------------- First step: Calibrate time-dependant volatilities functions ------------- \n");
    printf("------- We try to minimize the error function using an optimization algorithm ------- \n");
    LmmPiterbarg_Calib_vol(LmmPiterbarg);
    printf("--------------------------------- Model Calibration --------------------------------- \n");
    printf("------------- Second step: Calibrate time-dependant skews functions ------------- \n");
    printf("------- We try to minimize the error function using an optimization algorithm ------- \n");
    LmmPiterbarg_Calib_skew(LmmPiterbarg);

    FreeLmmPiterbarg(&LmmPiterbarg);

    return 0;
}


