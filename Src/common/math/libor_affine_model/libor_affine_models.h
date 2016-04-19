
#ifndef _LIBOR_AFFINE_MODELS_H
#define _LIBOR_AFFINE_MODELS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_vector.h"
#include "pnl/pnl_mathtools.h"

///******************* CIR 1d Model*******************///
void phi_psi_cir1d(PnlVect *ModelParams, double t, dcomplex u, dcomplex *phi_i, dcomplex *psi_i);
double MaxMgfArg_cir1d(PnlVect *ModelParams, double T);

///******************* Gamma-OU 1d Model*******************///
void phi_psi_gou1d(PnlVect *ModelParams, double t, dcomplex u, dcomplex *phi_i, dcomplex *psi_i);
double MaxMgfArg_gou1d(PnlVect *ModelParams, double T);


#endif
