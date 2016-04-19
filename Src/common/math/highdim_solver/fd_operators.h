#ifndef FD_OPERATORS
#define FD_OPERATORS

#include "fd_operators_common.h"

int FDOperatorJamInit(FDOperatorJam *, unsigned);
void FDOperatorJamFree(FDOperatorJam *);
void FDOperatorJamCoMatricesFillerSet(
                             FDSolverCoMatricesFiller *cmf,
                             FDOperatorJamCoMatricesFillerData *data,
                             FDOperatorJamCoMatricesFillerEqDef_t def,
                             FDOperatorJamCoMatricesFillerEqApply_t apply,
                             void *eq_data
                                     );

#endif

