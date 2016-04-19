#include "objFunction.h"
#include "spline.h"
#include "DupirePDE.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

double F(double *sigma_param, double lambda, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, int n, int m, int N, int M, double S_0, double r, double q, double theta, int optionType, struct marketData *optionPrices){
/* evaluates the cost function in sigma_param : F(sigma) = G(sigma)+lambda*F1(sigma)
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

  return G(sigma_param,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid,n,m,N,M,S_0,r,q,theta,optionType,optionPrices) + lambda*F1(sigma_param,y_coarseGrid,T_coarseGrid,n,m,optionType);

}



double F1(double *sigma_param, double *y_coarseGrid, double *T_coarseGrid, int n, int m, int optionType){
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
   



  double **deriv_y;
  double **deriv_T;
  double **deriv_yT;
  double **sigmaCoarseGrid;
  struct derivData *interpolData;
  double sum;
  int i,j,k,l;

  sigmaCoarseGrid = (double **) malloc((n+1)*sizeof(double *));

  deriv_y = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    deriv_y[i] = (double *) malloc((m+1)*sizeof(double));

  deriv_T = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    deriv_T[i] = (double *) malloc((m+1)*sizeof(double));

  deriv_yT = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    deriv_yT[i] = (double *) malloc((m+1)*sizeof(double));


  interpolData = (struct derivData *) malloc(sizeof(struct derivData));
  interpolData->deriv_y_0 = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_y_n = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_T_0 = (double *) malloc((n+1)*sizeof(double));
  interpolData->deriv_T_m = (double *) malloc((n+1)*sizeof(double));

  for (k=0;k<n+1;k++)
    sigmaCoarseGrid[k] = (double *) malloc((m+1)*sizeof(double));
  for (k=0;k<m+1;k++)
    for (l=0;l<n+1;l++)
      sigmaCoarseGrid[l][k] = sigma_param[k*(n+1)+l];
  
  for (k=(n+1)*(m+1);k<(n+1)*(m+1)+m+1;k++)
    interpolData->deriv_y_0[k-(n+1)*(m+1)] = sigma_param[k];

  for (k=(n+2)*(m+1);k<(n+2)*(m+1)+m+1;k++)
    interpolData->deriv_y_n[k-(n+2)*(m+1)] = sigma_param[k];
  
  for (k=(n+3)*(m+1);k<(n+3)*(m+1)+n+1;k++)
    interpolData->deriv_T_0[k-(n+3)*(m+1)] = sigma_param[k];

  for (k=(n+3)*(m+1)+n+1;k<(n+3)*(m+1)+2*(n+1);k++)
    interpolData->deriv_T_m[k-((n+3)*(m+1)+n+1)] = sigma_param[k];

  interpolData->deriv_yT_00 = sigma_param[(n+3)*(m+1)+2*(n+1)];
  interpolData->deriv_yT_n0 = sigma_param[(n+3)*(m+1)+2*(n+1)+1];
  interpolData->deriv_yT_0m = sigma_param[(n+3)*(m+1)+2*(n+1)+2];
  interpolData->deriv_yT_nm = sigma_param[(n+3)*(m+1)+2*(n+1)+3];

  interpolData->n = n;
  interpolData->m = m;
  
  /* sigma_param is now stored in sigmaCoarseGrid (values of sigma) and interpolData (derivatives) */

  /* computation of the derivatives on the coarse grid : deriv_y,deriv_T,deriv_yT */
  derivativesGrids(deriv_y,deriv_T,deriv_yT,sigmaCoarseGrid,y_coarseGrid,T_coarseGrid,interpolData);

  /* computation of the sum of squares of derivatives */
  sum = 0;
  for (i=0;i<n+1;i++)
    for (j=0;j<m+1;j++)
      sum = sum + pow(deriv_y[i][j],2) + pow(deriv_T[i][j],2);

  /* free memory */
  for (i=0;i<n+1;i++)
    free(deriv_y[i]);
  free(deriv_y);

  for (i=0;i<n+1;i++)
    free(deriv_T[i]);
  free(deriv_T);

 for (i=0;i<n+1;i++)
    free(deriv_yT[i]);
  free(deriv_yT);
 
  for (i=0;i<n+1;i++)
    free(sigmaCoarseGrid[i]);
  free(sigmaCoarseGrid);

  free(interpolData->deriv_y_0);
  free(interpolData->deriv_y_n);
  free(interpolData->deriv_T_0);
  free(interpolData->deriv_T_m);
  free(interpolData);

  /* return F1(sigma_param) = sum */
  return sum;  

}



double G(double *sigma_param, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, int n, int m, int N, int M, double S_0, double r, double q, double theta, int optionType, struct marketData *optionPrices){
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
   
  
  double **prices;
  double **sigmaFineGrid,**sigmaCoarseGrid;
  struct derivData *interpolData;
  int i,j,k,l;
  int ok;
  double sum,price,price1,price2;
  struct marketData *pt_data;
  
  /* memory allocation and initialization */
  /* reorganization of sigma_param (vector of parameters) in 
     - sigmaCoarseGrid for the values of sigma on the fine grid
     - interpolData for the derivatives on the borders and on the corners */

  sigmaCoarseGrid = (double **) malloc((n+1)*sizeof(double *));
  for (k=0;k<n+1;k++){
    sigmaCoarseGrid[k] = (double *) malloc((m+1)*sizeof(double));
    for (l=0;l<m+1;l++)
      sigmaCoarseGrid[k][l] = sigma_param[l*(n+1)+k];
  }
  sigmaFineGrid = (double **) malloc((N+1)*sizeof(double *));
  for (k=0;k<N+1;k++)
    sigmaFineGrid[k] = (double *) malloc((M+1)*sizeof(double));

  prices = (double **) malloc((N+1)*sizeof(double *));
  for (k=0;k<N+1;k++)
    prices[k] = (double *) malloc((M+1)*sizeof(double));


  interpolData = (struct derivData *) malloc(sizeof(struct derivData));
  interpolData->deriv_y_0 = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_y_n = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_T_0 = (double *) malloc((n+1)*sizeof(double));
  interpolData->deriv_T_m = (double *) malloc((n+1)*sizeof(double));

  for (k=(n+1)*(m+1);k<(n+1)*(m+1)+m+1;k++)
    interpolData->deriv_y_0[k-(n+1)*(m+1)] = sigma_param[k];

  for (k=(n+2)*(m+1);k<(n+2)*(m+1)+m+1;k++)
    interpolData->deriv_y_n[k-(n+2)*(m+1)] = sigma_param[k];
  
  for (k=(n+3)*(m+1);k<(n+3)*(m+1)+n+1;k++)
    interpolData->deriv_T_0[k-(n+3)*(m+1)] = sigma_param[k];

  for (k=(n+3)*(m+1)+n+1;k<(n+3)*(m+1)+2*(n+1);k++)
    interpolData->deriv_T_m[k-((n+3)*(m+1)+n+1)] = sigma_param[k];

  interpolData->deriv_yT_00 = sigma_param[(n+3)*(m+1)+2*(n+1)];
  interpolData->deriv_yT_n0 = sigma_param[(n+3)*(m+1)+2*(n+1)+1];
  interpolData->deriv_yT_0m = sigma_param[(n+3)*(m+1)+2*(n+1)+2];
  interpolData->deriv_yT_nm = sigma_param[(n+3)*(m+1)+2*(n+1)+3];

  interpolData->n = n;
  interpolData->m = m;

  /* interpolation of sigmaCoarseGrid to sigmaFineGrid */
  interpole(sigmaFineGrid,sigmaCoarseGrid,n,m,N,M,interpolData,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid);

  /* with sigmaFineGrid, we solve Dupire PDE */
  solve(optionType,prices,S_0,N,M,r,q,theta,f,sigmaFineGrid,y_fineGrid,T_fineGrid);

  sum = 0;
  pt_data = optionPrices;
  while (pt_data != NULL){
    /* computation of the closest indices */
    k = find_index(log(pt_data->strike),y_fineGrid,N,0);
    l = find_index(pt_data->maturity,T_fineGrid,M,0);
    /* bilinear interpolation */
    price1 = prices[k+1][l] - (y_fineGrid[k+1]-log(pt_data->strike))*(prices[k+1][l]-prices[k][l])/(y_fineGrid[k+1]-y_fineGrid[k]);
    price2 = prices[k+1][l+1] - (y_fineGrid[k+1]-log(pt_data->strike))*(prices[k+1][l+1]-prices[k][l+1])/(y_fineGrid[k+1]-y_fineGrid[k]);
    price = price2 - (T_fineGrid[l+1]-pt_data->maturity)*(price2-price1)/(T_fineGrid[l+1]-T_fineGrid[l]);
    sum = sum + pow(pt_data->price-price,2);
    /* update of pt_data */
    pt_data = pt_data->next;
  }

  /* free memory */
  for (k=0;k<N+1;k++)
    free(sigmaFineGrid[k]);
  free(sigmaFineGrid);
  
  for (k=0;k<N+1;k++)
    free(prices[k]);
  free(prices);
  
  for (k=0;k<n+1;k++)
    free(sigmaCoarseGrid[k]);
  free(sigmaCoarseGrid);
  
  free(interpolData->deriv_y_0);
  free(interpolData->deriv_y_n);
  free(interpolData->deriv_T_0);
  free(interpolData->deriv_T_m);
  free(interpolData);


  /* return G(sigma_param) */
  return sum/2;

}
