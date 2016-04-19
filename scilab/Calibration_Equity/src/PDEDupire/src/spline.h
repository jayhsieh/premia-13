#ifndef _SPLINE
#define _SPLINE 

#include "data.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.
*/

/**************types definitions**********************/

struct derivData{
  int n; /*number of space steps*/
  int m; /*number of time steps*/ 
  double *deriv_y_0; /*array of derivatives of sigma wrt y in i=0 and j=0...m*/
  double *deriv_y_n; /*array of derivatives of sigma wrt y in i=n and j=0...m*/
  double *deriv_T_0; /*array of derivatives of sigma wrt t in j=0 and i=0...n*/
  double *deriv_T_m; /*array of derivatives of sigma wrt t in j=m and i=0...n*/
  double deriv_yT_00; /*cross derivative in i=0,j=0 */
  double deriv_yT_0m; /*cross derivative in i=0,j=m */
  double deriv_yT_n0; /*cross derivative in i=n,j=0 */
  double deriv_yT_nm; /*cross derivative in i=n,j=m */
};



/**************functions definitions**********************/

double sigma_func(double S, double t);
/* OUTPUT
   - returns sigma(S,t)
   INPUTS
   - S
   - t   */

void buildCoarseGrids(double *y_coarseGrid, double *T_coarseGrid,int n,int m, double t_0, double T_max, double y_min, double y_max);
/* builds the arrays of discretized values of y and T for the coarse grid
   OUTPUTS:
   - y_coarseGrid = array of discretized values of y for the coarse grid
   - T_coarseGrid = array of discretized values of y for the coarse grid
   - (n,m) = size of the coarse grid
   - t_0 = time origin
   - T_max = time horizon
   - y_min = log(S_min)
   - y_max = log(S_max)
*/


void buildFineGrid(double *y_fineGrid, double *T_fineGrid, int N, int M, double t_0, double T_max, double y_min, double y_max, double S_0, int gridType);
/* builds the arrays of discretized values of y and T for the fine grid
   OUTPUTS:
   - y_fineGrid = array of discretized values of y for the fine grid
   - T_fineGrid = array of discretized values of y for the fine grid
   - (N,M) = size of the fine grid
   - t_0 = time origin
   - T_max = time horizon
   - y_min = log(S_min)
   - y_max = log(S_max)
   - gridType = type of the y_fineGrid (0 for regular, 1 for tanh)
*/


double invtanh(double x, double a, double b);
/*fonction invtanh inverse de a*tanh(x-b) */
/* returns y so that a*tanh(y-b) = x     */


void derivativesGrids (double **deriv_y, double **deriv_T, double **deriv_yT, double **sigmaCoarseGrid, double *y_coarseGrid, double *T_coarseGrid, struct derivData *data);
/* AIM : computes the derivatives wrt y, T, and the cross derivatives wrt y and T on all the points of the coarse grid 
   OUTPUTS:
   deriv_y = matrix of the derivatives wrt y                     
   deriv_T = matrix of the derivatives wrt T                     
   deriv_yT =  matrix of the cross derivatives wrt y and T    
   INPUT PARAMETERS :                                                                                              
   sigmaCoarseGrid = coarse grid of the values of sigma                                            
   y_coarseGrid = array of the discretized values of y                       
   t_coarseGrid = array of the discretized values of t                     
   data = data (derivatives on the borders, cross derivatives on the corners) of the interpolation problem */ 



void interpCoeffs(double ****coeffs, double *y_coarseGrid, double *T_coarseGrid, double **deriv_y, double **deriv_T, double **deriv_yT, double **sigmaCoarseGrid, int n, int m);
/* OUTPUTS
   coeffs = 4-dim array containing the coefficients of the interp. function    
   INPUT PARAMETERS :                                                                              
   y_coarseGrid = array of the discretized values of y                             
   T_coarseGrid = array of the discretized values of T                              
   deriv_y = matrix of the derivatives wrt y                                                 
   deriv_T = matrix of the derivatives wrt T                                                 
   deriv_yT =  matrix of the cross derivatives wrt y and T                                         
   sigmaCoarseGrid = coarse grid of the values of sigma                                                       
   n = number of space (price) steps                                                         
   m = number of time steps                                                                  */






void computeK(double **K, int i, int j, double **deriv_t, double **deriv_y, double **deriv_yT, double **sigmaCoarseGrid);
/*OUTPUT
  - K (=Kij) =  matirx of size (4,4)
  INPUT
  - i,j = indices of the matrix Kij to be computed
  - deriv_y = matrix of the derivatives wrt y                                                 
  - deriv_T = matrix of the derivatives wrt T                                                 
  - deriv_yT =  matrix of the cross derivatives wrt y and T 
  - sigmaCoarseGrid = coarse grid of the values of sigma */


void computeH(double **H, double delta_T);
/*OUTPUT
  - H =  matirx of size (4,4)
  INPUT
  - delta_T = time step */


void computeP(double **P, double delta_y);
/*OUTPUT
  - P =  matirx of size (4,4)
  INPUT
  - delta_y = space step */


void matMult(double **AB, double **A, double **B, int n, int m, int p);
/* OUPUT:
   - AB = A*B of size (m,p)
   INPUTS
   - A = matrix of size (m,n)
   - B = matrix of size (n,p) */



int find_index(double x, double *y_grid, int size, int start);
/* OUTPUT
   - returns the index i such that grid[i]<=x<grid[i+1]
   INPUT
   - x (x = y or T)
   - grid = grid of discretized values of x, grid = x_0,...,x_n (x = y or T)
   - size = size of grid - 1 
   - start = index at which we start the search of i */


double interpEval2_dy(double ****coeffs, double y, double T, int i, int j, double y_i, double T_j);
/* OUTPUT
   - returns the evaluation of the derivative wrt to y at the point (y,T) 
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_i = y_coarseGrid[i]
   - T_j = T_coarseGrid[j]          */

  double interpEval2_dT(double ****coeffs, double y, double T, int i, int j, double y_i, double T_j);
/* OUTPUT
   - returns the evaluation of the derivative wrt to T at the point (y,T) 
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_i = y_coarseGrid[i]
   - T_j = T_coarseGrid[j]          */

double interpEval(double ****coeffs, double y, double T, double *y_coarseGrid, double *T_coarseGrid, int n, int m);
/* OUTPUT
   - returns the evaluation of the interpolation function at the point (y,T)
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_coarseGrid = grid of discretized values of y, grid_y = y_0,...,y_n
   - T_coarseGrid = grid of discretized values of T, grid_T = T_0,...,T_m
   - n = size of y_coarseGrid - 1 
   - m : size of T_coarseGrid - 1 */ 


double interpEval2(double ****coeffs, double y, double T, int i, int j, double y_i, double T_j);
/* OUTPUT
   - returns the evaluation of the interpolation function at the point (y,T)
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_i = 
   - T_j =           */



void discretizeSigma(double (*volatility)(double,double), double **sigmaCoarseGrid, int n, int m, double *y_coarseGrid, double *T_coarseGrid);
/* OUTPUT
   - sigmaCoarseGrid = coarse grid of discretized values of sigma
   INPUTS
   - volatility = function defining the volatility
   - n,m = size of the coarse grid
   - y_coarseGrid = array of dicretized values y (for the coarse grid)
   - T_coarseGrid = array of dicretized values T (for the coarse grid) */
   


void interpole(double **sigmaFineGrid, double **sigmaCoarseGrid, int n, int m, int N, int M, struct derivData *data, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid);
/* OUTPUT
   - sigmaFineGrid = fine grid of the interpolated values of sigma
   INPUTS
   - sigmaCoarseGrid = coarse grid of discretized values of sigma
   - (N,M) = size of the fine grid
   - (n,m) = size of the coarse grid
   - data = data (derivatives on the borders, cross derivatives on the corners) of the interpolation problem
   - y_coarseGrid = array of dicretized values y (for the coarse grid)
   - T_coarseGrid = array of dicretized values T (for the coarse grid)
   - y_fineGrid = array of dicretized values y (for the fine grid)
   - T_fineGrid = array of dicretized values T (for the fine grid)  */ 


void ddlGrid1ToGrid2(double **pt_ddl2, double *ddl1, int n1, int m1, int n2, int m2, double *y_coarseGrid1, double *T_coarseGrid1, double *y_coarseGrid2, double *T_coarseGrid2);

void printDerivData(struct derivData *d);
/* prints the data in d */


#endif
