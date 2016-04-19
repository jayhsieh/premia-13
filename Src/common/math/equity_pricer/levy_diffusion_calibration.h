#ifndef LEVY_DIFFUSION_CALIBRATION_H_
#define LEVY_DIFFUSION_CALIBRATION_H_
#include <stdio.h>
#include <stdlib.h>

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include "finance_tool_box.h"


extern void HestonConstraints(PnlVect *res, const PnlVect *x, void *params);
extern double QuadraticError_ForHeston(const PnlVect *GenerationParams,void *data);
extern void BatesConstraints(PnlVect *res, const PnlVect *x, void *params);
extern double QuadraticError_ForBates(const PnlVect *GenerationParams,void *data);
extern void BNSConstraints(PnlVect *res, const PnlVect *x, void *params);
extern double QuadraticError_ForBNS(const PnlVect *GenerationParams,void *data);


extern void NIGGammaOUConstraints(PnlVect *res, const PnlVect *x, void *params);
extern double QuadraticError_ForNIGGammaOU(const PnlVect *GenerationParams,void *data);

#endif // LEVY_DIFFUSION_CALIBRATION_H

