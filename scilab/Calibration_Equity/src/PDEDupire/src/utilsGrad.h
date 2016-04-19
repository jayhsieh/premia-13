#ifndef _UTILS_GRAD
#define _UTILS_GRAD

#include "solveSystem.h"
#include "data.h"
#include "sparse.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

void buildTranspDerivPhi(struct sparseMat **pt_transp_deriv_phi_j, int j, int n_d, int N, int M, int nbMaturites, double *y_fineGrid, double *T_fineGrid, struct marketData **marketOptionPrices);

void buildAdjointTridiagSyst(struct tridiagSystem *tridiagSyst, int j, struct tridiag *A, double *u_prev, double theta, double timeStep, double timeStep2, struct sparseMat *transp_deriv_phi, double *phi_U, double *U_tilde, int n_d, int M);

void findClosestIndices(int *ok, int *k, int *l, double K, double T, double *y_grid, double *t_grid, int N, int M);

void matMultVect(double *A, double **B, double *C, int nbLinesB, int nbColB);
/* A = B*C */


void VectPlus(double *A, double *B, int nbLines);
/* A = A+B */

double g(double y, double **U, struct marketData **data, int nbMaturites, int N, int M, double T_max, double *y_grid,  double *t_grid);

void compute_R(struct sparseMat **R, int N, int M, double **invA_times_B, double **invC_times_D, int **B_indices, int **D_indices, struct tridiag *D, double *T_coarseGrid);

void compute_Li(struct sparseMat **pt_Li, int i, int nbParamSigma, int n, int m, double **invA_times_B, int **B_indices);

void compute_Pi(struct sparseMat **pt_Pi, int i, int nbParamSigma, int n, int m, int *H0_indices, int *Hm_indices, double **invA_times_B);
void compute_deriv_Au_j(struct sparseMat **deriv_Au_j, int j, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, int N, int M, int n, int m, struct sparseMat *R, double **sigmaFineGrid, double **sigmaCoarseGrid, double **U);

void compute_transp_Vkl(struct sparseMat ***transp_Vkl, int k, int l, int N, int M, double delta_y, double delta_T);

#endif
