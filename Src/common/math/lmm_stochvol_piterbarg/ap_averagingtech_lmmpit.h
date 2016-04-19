#ifndef _Ap_Averaging_LMM_PITERBARG
#define _Ap_Averaging_LMM_PITERBARG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lmm_stochvol_piterbarg.h"
#include "math/read_market_zc/InitialYieldCurve.h"


double cf_lmm_stochvol_piterbarg_swpt(StructLmmPiterbarg *LmmPiterbarg, double Tn, double Tm, double period, double swaption_strike, double Nominal, int Payer_Receiver, int FlagClosedFormula_in);

double cf_lmm_stochvol_piterbarg_capfloor(StructLmmPiterbarg *LmmPiterbarg, double Tn, double Tm, double period, double cap_strike, double Nominal, int flag_capfloor, int FlagClosedFormula_in);


#endif
