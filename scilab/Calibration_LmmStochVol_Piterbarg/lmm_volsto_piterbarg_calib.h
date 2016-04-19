#ifndef _LMM_VOLSTO_PITERBARG_CALIB
#define _LMM_VOLSTO_PITERBARG_CALIB

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "pnl/pnl_mathtools.h"
#include "InitialYieldCurve.h"

double func_to_min_vols(const PnlVect *x, void *LmmPiterbarg);

double func_to_min_skews(const PnlVect *x, void *LmmPiterbarg);

void LmmPiterbarg_Calib_vol(StructLmmPiterbarg *LmmPiterbarg);

void LmmPiterbarg_Calib_skew(StructLmmPiterbarg *LmmPiterbarg);
#endif
