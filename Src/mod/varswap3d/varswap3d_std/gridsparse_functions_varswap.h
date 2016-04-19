#ifndef GRIDSPARSE_FUNCTIONS_VARSWAP_H
#define GRIDSPARSE_FUNCTIONS_VARSWAP_H
#include "varswap3d_std.h"
#include "math/equity_pricer/pde_tools.h"
#include "math/equity_pricer/gridsparse_constructor.h"

typedef struct SVSSparseOp{
  VARSWAP3D_MOD * Model;
  double theta_time_scheme;
  PremiaPDETimeGrid * TG;
  GridSparse *G;
  PnlVect  *PC,* V_tmp0,* V_tmp1,* V_tmp2,* V_tmp3,
    * V_tmp4,* V_tmp5,* V_tmp6,* V_tmp7,*V_volatility;
}SVSSparseOp;


//////////////////////////////////////////////////////////
// Stochastic Variance Swap Model Operator on Sparse Grid
//////////////////////////////////////////////////////////

extern SVSSparseOp * svs_sparse_operator_create(int lev,int N_T,VARSWAP3D_MOD * M);
extern void svs_sparse_operator_free(SVSSparseOp ** Op);
extern void GridSparse_Solve_svs(SVSSparseOp * Op,PnlVect * Vres);
#endif
