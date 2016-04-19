#ifndef DupirePDE
#define DupirePDE 

#include "solveSystem.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.
*/

/************************functions definitions****************************/


void solve(int optionType, double **res, double S_0, int N, int M, double r, double q, double theta, double (*f)(int,double,double), double **sigmaFineGrid, double *y_fineGrid, double *T_fineGrid);
/* solves the EDP using the finite difference method 
   OUTPUT :
   - res : fine grid of prices to be computed (size N*M)
   INPUTS :
   - optionType : type of the option (1 for call, 0 for put) 
   - S_0 : price of the asset at t_0                         
   - N : number of space steps of the fine grid                              
   - M : number of time steps of the fine grid                               
   - r : RF rate                                             
   - q : dividends                                           
   - theta : parameter of the finite difference scheme       
   - sigma : fine grid of values of sigma                          
   - y_fineGrid : discretized values of y  (for the fine grid)                     
   - T_fineGrid : discretized values of T  (for the fine grid)        */



double f(int optionType, double y, double S_0);
/* condition for T=t_0                    */
/* OUTPUT :
   - returns f(y,S_0) (condition for T=t_0
   INPUTS :
   - optionType = 1 for a call, 0 for a put
   - y : log of the price of the asset
   - S_0 : price of the asset at t_0  */




void buildOperator(struct tridiag *A, double r, double q, double S_0, int j, int N, double **sigmaFineGrid, double *y_fineGrid);
/*  builds the discretized operator A            
    OUPUT:
    - A : discretized operator (tridiagonal matrix represented by a struct of 3 vectors) of size N-1
    INPUTS:
    - r : RF rate
    - q : dividends 
    - S_0 : price of the asset at t_0 
    - j : index of time step
    - N : number of space (price) steps for the fine grid
    - sigmaFineGrid : fine grid of the values of sigma
    - y_fineGrid :discretized values of y  (for the fine grid)  */




void buildTridiagSystem(int optionType, struct tridiagSystem *S, struct tridiag *A, struct tridiag *A_prev, double *u_prev, double k, double theta, double condLim); 
/* builds the data of the tridiagonal system S  (Mat*X = b, with Mat tridiagonal matrix, b right hand side vector)     
   OUTPUT:
   - S : tridiagonal system 
   INPUTS:
   - optionType: type of the option (1 for call, 0 for put)
   - A : current discretized operator                              
   - A_prev : previous discretized operator                         
   - u_prev : vector of prices computed in the previous iteration   
   - k : size of time step                                          
   - theta : parameter of the finite difference scheme     
   - condLim : limit condition (depending on the time of the option) */



void optionPricer(int *pt_i, int *pt_j, double T, double K, double **res, double *y_fineGrid, double *T_fineGrid, int N, int M);
/*returns the indices of Rij of the fine grid corresponding to maturity=T and strike=K */
/* OUTPUT:
   - pt_i,pt_j:indices of Rij corresponding to maturity=T and strike=K 
   INPUTS:
   - T,K : maturity,strike of the option we want to price
   - res : gird of prices (solution of Dupire PDE)
   - y_fineGrid : array of discretized values of y for the fine grid
   - T_fineGrid : array of discretized values of T for the fine grid 
   - N,M : size of the fine grid */

double computeImpsig(double S_0, double r, double q, int optionType, double t_0, double K, double T, double V, double tol);


double BSprice(int optionType, double S, double t, double sigma, double r, double q, double T, double K);
/*       computes the price of a european  option      */
/*         (with sigma constant) with BS model         */
/* PARAMETERS :                                        */
/* optionType = 1 for a call, 0 for a put              */
/* S = spot price at time t                            */
/* t = present time                                    */
/* sigma = constant volatility of the underlying asset */
/* r = RF rate                                         */
/* q = dividends                                       */
/* T = maturity of the option                          */
/* K strike price of the option                        */




double BSpriceTimeDepVol(double S, double t, double r, double q, double T, double K, double (*intSigma)(double,double));
/* AIM : computes the price of an european call option        */
/*         (with sigma time-dependent) with BS model          */
/* PARAMETERS :                                               */
/* S = spot price at time t                                   */
/* t = present time                                           */
/* sigma = constant volatility of the underlying asset        */
/* r = RF rate                                                */
/* q = dividends                                              */
/* T = maturity of the option                                 */
/* K = strike price of the option                             */
/* intSigma = function that computes the integral of sigma^2  */


double distrib(double d);
/* AIM : approximates the distribution function of a random variable ~ N(0,1) */


double intSigma(double t, double T);
/* AIM : computes the integral of sigma^2 over [t,T], divided */
/* by (T-t), when sigma is time-dependent                     */


double rand_nb();
double priceSimul(double t, double t_0, double S_0);
double priceEstim(double t, double t_0, double S_0, double strike, int nbSimul);

#endif
