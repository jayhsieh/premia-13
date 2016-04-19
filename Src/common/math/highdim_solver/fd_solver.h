#ifndef FD_SOLVER_H
#define FD_SOLVER_H

#include "fd_solver_common.h"

extern int FDSolverInit(FDSolver *solver,
                        FDSolverVectorFiller *icf,
                        FDSolverCoMatricesFiller *AcBcf,
                        FDSolverCoMatricesFiller *AnBnf);

extern void FDSolverFree(FDSolver *);

extern int FDSolverGet(FDSolver *, double *);

int FDSolverResetMatrices(FDSolver *solver, FDSolverCoMatricesFiller *AcBcf,
                          FDSolverCoMatricesFiller *AnBnf);

int FDSolverStep(FDSolver *solver);

#endif

