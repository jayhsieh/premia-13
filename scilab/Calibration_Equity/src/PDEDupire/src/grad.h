#ifndef _GRAD_
#define _GRAD_ 

#include "data.h"
#include "solveSystem.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

void computeGrad_F(double *grad_F, double *sigma_param, double lambda, int N, int M, int n, int m, double *y_fineGrid, double *T_fineGrid, double *y_coarseGrid, double *T_coarseGrid, double r, double q, double S_0, int optionType, struct marketData **marketOptionPrices, int nbMaturites);



void computeGrad_F1(double **matGradF1, double **invA_times_B, double **invC_times_D, int **D_indices, int **B_indices, struct tridiag *D, int n, int m, double *y_coarseGrid, double *T_coarseGrid);


void computeGrad_G(double *pt_grad_G, double *sigma_param, int N, int M, int n, int m, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, double r, double q, int optionType, double S_0, double **invA_times_B, double **invC_times_D, int **D_indices, int **B_indices, struct tridiag *D, struct marketData **marketOptionPrices, int nbMaturites);

#endif 
