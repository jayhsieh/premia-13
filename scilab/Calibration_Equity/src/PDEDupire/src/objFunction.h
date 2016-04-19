#include "data.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

double F(double *sigma_param, double lambda, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, int n, int m, int N, int M, double S_0, double r, double q, double theta, int optionType, struct marketData *optionPrices);
/* evaluates the cost function in sigma_param 
   OUTPUT:
   - returns F(sigma_param)
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
   


double F1(double *sigma_param, double *y_coarseGrid, double *T_coarseGrid, int n, int m, int optionType);
/* evaluates F1 in sigma_param 
   OUTPUT:
   - returns F1(sigma_param)
   INPUTS:
   - sigma_param : point where we evaluate F1
   - y_coarseGrid : discretized values of y  (for the coarse grid)                     
   - T_coarseGrid : discretized values of T  (for the coarse grid)        
   - (n,m) = size of the coarse grid
   - optionType : type of the option (1 for call, 0 for put) 
*/
   


double G(double *sigma_param, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, int n, int m, int N, int M, double S_0, double r, double q, double theta, int optionType, struct marketData *optionPrices);
/* evaluates G in sigma_param 
   OUTPUT:
   - returns G(sigma_param) :  G(sigma_param) = 1/2 || U(sigma_param) - optionPrices ||_2^2
   INPUTS:
   - sigma_param : point where we evaluate G
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
   
