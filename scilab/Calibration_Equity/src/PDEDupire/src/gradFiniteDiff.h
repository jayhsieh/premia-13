#ifndef _GRAD_FINITE_DIFF_
#define _GRAD_FINITE_DIFF_ 

#include "data.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

void computeGrad_F_finite_diff(double *grad_F_finite_diff, double *sigma_param, double lambda, int n, int m, int N, int M, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, double r, double q, double S_0, double theta, int optionType, struct marketData *marketOptionPrices);
/* computes the gradient of F with finite diff method
   OUTPUT:
   - grad_F_finite_diff (gradient evaluated in sigma_param)
   INPUTS:
   - sigma_param : point where we evaluate F
   - lambda : coeff of F1 (F=G+lambda*F1)
   - y_coarseGrid : discretized values of y  (for the coarse grid)                     
   - T_coarseGrid : discretized values of T  (for the coarse grid)        
   - y_fineGrid : discretized values of y  (for the fine grid)                     
   - T_fineGrid : discretized values of T  (for the fine grid)
   - (n,m) = size of the coarse grid
   - (N,M) = size of the fine grid
   - S_0 = price of the underlying asset at t_0
   - r : RF rate                                             
   - q : dividends
   - theta : parameter of the finite difference scheme       
   - optionType : type of the option (1 for call, 0 for put) 
   - optionPrices : market prices of options
*/


double gradScalarDir(double *sigma_param, double lambda, double *direction, double h, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, int n, int m, int N, int M, double S_0, double r, double q, double theta, int optionType, struct marketData *optionPrices);
/* computes <gradF(sigma_param),direction> = (1/h)(F(sigma_param + direction)-F(sigma_param))
   OUTPUT:
   - returns <gradF(sigma_param),direction>
   INPUTS:
   - sigma_param : point where we evaluate F
   - lambda : coeff of F1 (F=G+lambda*F1)
   - direction : direction (vector) in which we want to compute the grad
   - h : step 
   - y_coarseGrid : discretized values of y  (for the coarse grid)                     
   - T_coarseGrid : discretized values of T  (for the coarse grid)        
   - y_fineGrid : discretized values of y  (for the fine grid)                     
   - T_fineGrid : discretized values of T  (for the fine grid)
   - (n,m) = size of the coarse grid
   - (N,M) = size of the fine grid
   - S_0 = price of the underlying asset at t_0
   - r : RF rate                                             
   - q : dividends
   - theta : parameter of the finite difference scheme       
   - optionType : type of the option (1 for call, 0 for put) 
   - optionPrices : market prices of options
*/


#endif
