#ifndef _LMM_VOLSTO_PITERBARG_PRICING
#define _LMM_VOLSTO_PITERBARG_PRICING

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "InitialYieldCurve.h"
#include "lmm_volsto_piterbarg.h"

double cf_lmm_volsto_piterbarg_swpt(StructLmmPiterbarg *LmmPiterbarg, double Tn, double Tm, double period, double swaption_strike, double Nominal, int Payer_Receiver, int FlagClosedFormula_in);

#endif
