#include "DupirePDE.h"
#include "solveSystem.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.
*/

void solve(int optionType, double **res, double S_0, int N, int M, double r, double q, double theta, double (*f)(int,double,double), double **sigmaFineGrid, double *y_fineGrid, double *T_fineGrid){
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
   - T_fineGrid : discretized values of t  (for the fine grid)        */


  /* declarations */
  double *u_prev,*u_next;
  double h_i, h_i_1;
  int i,j;
  struct tridiag *op, *op_prev;
  struct tridiagSystem *S1;
  struct bidiagSystem *S2; 
  double x,h_0,h_1,alpha,alpha2,gamma,gamma2,condLim,k;
  double t_0,T_max,y_min,y_max;

  
  y_min = y_fineGrid[0];
  y_max = y_fineGrid[N];
  t_0 = T_fineGrid[0];
  T_max = T_fineGrid[M];

  /* memory allocation for the operator (of type struct tridiag *) */
  op = (struct tridiag *) malloc(sizeof(struct tridiag));
  op->subdiag = (double *) malloc((N-2)*sizeof(double));
  op->diag = (double *) malloc((N-1)*sizeof(double));
  op->updiag = (double *) malloc((N-2)*sizeof(double));
  op->size = N-1;
 
  /* memory allocation for the operator (of type struct tridiag *) */
  op_prev = (struct tridiag *) malloc(sizeof(struct tridiag));
  op_prev->subdiag = (double *) malloc((N-2)*sizeof(double));
  op_prev->diag = (double *) malloc((N-1)*sizeof(double));
  op_prev->updiag = (double *) malloc((N-2)*sizeof(double));
  op_prev->size = N-1;
  
  /* memory allocation for the tridiagonal system S1*/
  S1 = (struct tridiagSystem *) malloc(sizeof(struct tridiagSystem));
  S1->T = (struct tridiag *) malloc(sizeof(struct tridiag));
  S1->T->subdiag = (double *) malloc((N-2)*sizeof(double));
  S1->T->diag = (double *) malloc((N-1)*sizeof(double));
  S1->T->updiag = (double *) malloc((N-2)*sizeof(double));
  S1->T->size = N-1;
  S1->b = (double *) malloc((N-1)*sizeof(double));
  S1->size = N-1;

  /* memory allocation for the bidiagonal system S2*/
  S2 = (struct bidiagSystem *) malloc(sizeof(struct bidiagSystem));
  S2->T = (struct bidiag *) malloc(sizeof(struct bidiag));
  S2->T->subdiag = (double *) malloc((N-2)*sizeof(double));
  S2->T->diag = (double *) malloc((N-1)*sizeof(double));
  S2->T->size = N-1;
  S2->b = (double *) malloc((N-1)*sizeof(double));
  S2->size = N-1;

 
  /* memory allocation for u_prev and u_next */
  u_prev = (double *) malloc((N-1)*sizeof(double));
  u_next = (double *) malloc((N-1)*sizeof(double));


  /* initialization of u */
 
  if (optionType == 1)
    res[0][0] = S_0;
  else
    res[0][0] = 0; 
  
  for (i=1;i<N;i++){
    u_prev[i-1] = f(optionType,y_fineGrid[i],S_0);
    res[i][0] = u_prev[i-1];
  }
  res[N][0] = f(optionType,y_max,S_0);



  /*builds the initial discretized operator*/
  buildOperator(op_prev,r,q,S_0,0,N,sigmaFineGrid,y_fineGrid);
  for (j=1;j<=M;j++){
    buildOperator(op,r,q,S_0,j,N,sigmaFineGrid,y_fineGrid); /*builds the discretized operator*/ 
    /* initial condition for y=y_min*/
    h_0 = y_fineGrid[1]-y_fineGrid[0];
    h_1 = y_fineGrid[2]-y_fineGrid[1];
    k = T_fineGrid[j]-T_fineGrid[j-1];
    x = .5*(r-q+pow(sigmaFineGrid[1][j-1],2)/2);
    alpha = -x/h_0 - pow(sigmaFineGrid[1][j-1],2)/((h_1+h_0)*h_0);
    gamma = x/h_0 - pow(sigmaFineGrid[1][j-1],2)/((h_1+h_0)*h_0); 
    x = .5*(r-q+pow(sigmaFineGrid[1][j],2)/2);
    alpha2 = -x/h_0 - pow(sigmaFineGrid[1][j],2)/((h_1+h_0)*h_0);
    gamma2 =  x/h_0 - pow(sigmaFineGrid[1][j],2)/((h_1+h_0)*h_0);
    if (optionType==1) //call
      condLim = - S_0*k*(theta*alpha*exp(-q*(T_fineGrid[j-1]-t_0)) + (1-theta)*alpha2*exp(-q*(T_fineGrid[j]-t_0)));
    else //put
      condLim = - k*theta*gamma*(-S_0*exp(-q*(T_fineGrid[j-1]-t_0))+exp(y_max)*exp(-r*(T_fineGrid[j-1]-t_0))) + k*(1-theta)*gamma2*(-S_0*exp(-q*(T_fineGrid[j]-t_0))+exp(y_max)*exp(-r*(T_fineGrid[j]-t_0)));
      
    buildTridiagSystem(optionType,S1,op,op_prev,u_prev,k,theta,condLim); /*builds the tridiagonal system*/

    tridiagToBidiagSyst(S2,S1); /* changes the tridiag system to a bidiagonal one */
    solveSyst(u_next,S2); /* solves the bidiagonal system*/

  

    /*we update u_prev and stock the results in the matrix res */ 
    if (optionType == 1){ //call
      res[0][j] = S_0*exp(-q*(T_fineGrid[j]-t_0)); /* = limit condition for y = y_min for a call option*/
      res[N][j] = 0; /* = limit condition for y = y_max for a call option*/
    }else{
      res[0][j] = 0;/* = limit condition for y = y_min for a put option*/
      res[N][j] = -S_0*exp(-q*(T_fineGrid[j]-t_0))+exp(y_max)*exp(-r*(T_fineGrid[j]-t_0)); /* = limit condition for y = y_max for a put option*/
    }

    for (i=1;i<N;i++){
      h_i = y_fineGrid[i+1]-y_fineGrid[i];
      h_i_1 = y_fineGrid[i]-y_fineGrid[i-1];
      res[i][j] = u_next[i-1];
      u_prev[i-1] = u_next[i-1];
    }
    
    /* update of op_prev */
    affectOperator(op,op_prev);

  }

 
  /*free memory space*/
  free(S1->T->subdiag);
  S1->T->subdiag = NULL;
  free(S1->T->diag);
  S1->T->diag = NULL;
  free(S1->T->updiag);
  S1->T->updiag = NULL;
  free(S1->T);
  S1->T = NULL;
  free(S1->b);
  S1->b = NULL;
  free(S1);
  S1 = NULL;
  
  free(S2->T->subdiag);
  S2->T->subdiag = NULL;
  free(S2->T->diag);
  S2->T->diag = NULL;
  free(S2->T);
  S2->T = NULL;
  free(S2->b);
  S2->b = NULL;
  free(S2);
  S2 = NULL;

  free(op->subdiag);
  op->subdiag = NULL;
  free(op->diag);
  op->diag = NULL;
  free(op->updiag);
  op->updiag = NULL;
  free(op);
  op = NULL;

  free(op_prev->subdiag);
  op_prev->subdiag = NULL;
  free(op_prev->diag);
  op_prev->diag = NULL;
  free(op_prev->updiag);
  op_prev->updiag = NULL;
  free(op_prev);
  op_prev = NULL;


  free(u_prev);
  u_prev = NULL;
  free(u_next);
  u_next = NULL;

}


double f(int optionType, double y, double S_0){
/* condition for T=t_0                   
   OUTPUT :
   - returns f(y,S_0) (condition for T=t_0
   INPUTS :
   - optionType = 1 for a call, 0 for a put
   - y : log of the price of the asset
   - S_0 : price of the asset at t_0  */

  if (optionType==1)
    return (S_0-exp(y)>=0) ? S_0-exp(y) : 0;
  else
    return (exp(y)-S_0>=0) ? exp(y)-S_0 : 0; 
}




void buildOperator(struct tridiag *A, double r, double q, double S_0, int j, int N, double **sigmaFineGrid, double *y_fineGrid){
/*  builds the discretized operator A            
    OUPUT:
    - A : discretized operator (tridiagonal matrix) of size N-1
    INPUTS:
    - r : RF rate
    - q : dividends 
    - S_0 : price of the asset at t_0 
    - j : index of time step
    - N : number of space (price) steps for the fine grid
    - sigmaFineGrid : fine grid of the values of sigma
    - y_fineGrid :discretized values of y  (for the fine grid)  */                   

  
  double x,h_i,h_i_1;
  int i;

  i=1;
  h_i = y_fineGrid[i+1]- y_fineGrid[i];
  h_i_1 = y_fineGrid[i] - y_fineGrid[i-1];
  x = .5*(r-q+pow(sigmaFineGrid[i][j],2)/2);
  A->diag[i-1] = pow(sigmaFineGrid[i][j],2)/(h_i+h_i_1)*(1/h_i + 1/h_i_1) + x*(1/h_i_1 - 1/h_i) + q;
  A->updiag[i-1] = x/h_i - pow(sigmaFineGrid[i][j],2)/((h_i+h_i_1)*h_i);
  for (i=2;i<=N-2;i++){
    h_i = y_fineGrid[i+1]- y_fineGrid[i];
    h_i_1 = y_fineGrid[i] - y_fineGrid[i-1];
    x = .5*(r-q+pow(sigmaFineGrid[i][j],2)/2);
    A->subdiag[i-2] = -x/h_i_1 - pow(sigmaFineGrid[i][j],2)/((h_i+h_i_1)*h_i_1);
    A->diag[i-1] = pow(sigmaFineGrid[i][j],2)/(h_i+h_i_1)*(1/h_i + 1/h_i_1) + x*(1/h_i_1 - 1/h_i) + q;
    A->updiag[i-1] = x/h_i - pow(sigmaFineGrid[i][j],2)/((h_i+h_i_1)*h_i);
  }

  /* i=N-1 */
  h_i = y_fineGrid[N]- y_fineGrid[N-1];
  h_i_1 = y_fineGrid[N-1] - y_fineGrid[N-2];
  x = .5*(r-q+pow(sigmaFineGrid[N-1][j],2)/2);
  A->subdiag[N-3] = -x/h_i_1 - pow(sigmaFineGrid[N-1][j],2)/((h_i+h_i_1)*h_i_1);
  A->diag[N-2] = pow(sigmaFineGrid[N-1][j],2)/(h_i+h_i_1)*(1/h_i + 1/h_i_1) + x*(1/h_i_1 - 1/h_i) + q;


}





void buildTridiagSystem(int optionType, struct tridiagSystem *S, struct tridiag *A, struct tridiag *A_prev, double *u_prev, double k, double theta, double condLim){
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


  int size,i;


  size = A->size;

  S->b[0] = (1-theta*k*A_prev->diag[0])*u_prev[0] - theta*k*A_prev->updiag[0]*u_prev[1];

  if (optionType == 1) // call option
    S->b[0] = S->b[0] + condLim;
  
  S->T->subdiag[0] = (1-theta)*k*A->subdiag[0];
  S->T->diag[0] = 1+(1-theta)*k*A->diag[0];
  S->T->updiag[0] = (1-theta)*k*A->updiag[0];
  
  for (i=1;i<size-1;i++){
    S->b[i] = (-theta*k*A_prev->subdiag[i-1])*u_prev[i-1] + (1-theta*k*A_prev->diag[i])*u_prev[i] - theta*k*(A_prev->updiag[i])*u_prev[i+1];
    S->T->subdiag[i] = (1-theta)*k*A->subdiag[i];
    S->T->diag[i] = 1+(1-theta)*k*A->diag[i];
    S->T->updiag[i] = (1-theta)*k*A->updiag[i];
  }   

  S->T->diag[size-1] = 1+(1-theta)*k*A->diag[size-1];
  S->b[size-1] = -theta*k*A_prev->subdiag[size-2]*u_prev[size-2] + (1-theta*k*A_prev->diag[size-1])*u_prev[size-1];
  if (optionType == 0) // put option
    S->b[size-1] = S->b[size-1] + condLim;
} 


void optionPricer(int *pt_i, int *pt_j, double T, double K, double **res, double *y_fineGrid, double *T_fineGrid, int N, int M){
/*returns the indices of Rij of the fine grid corresponding to maturity=T and strike=K */
/* OUTPUT:
   - pt_i,pt_j:pointers on the indices of Rij corresponding to maturity=T and strike=K 
   INPUTS:
   - T,K : maturity,strike of the option we want to price
   - res : gird of prices (solution of Dupire PDE)
   - y_fineGrid : array of discretized values of y for the fine grid
   - T_fineGrid : array of discretized values of T for the fine grid 
   - N,M : size of the fine grid */


  int i,j;
  i=1;
  j=1;
  while (i<=N && exp(y_fineGrid[i])<K)
    i++;
  if (i==N+1 || K<0)
	printf("You chose a strike price that is not within the bounds\n");
  else{
    while (j<=M && T_fineGrid[j]<=T)
      j++;
    if (j==M+1 && T_fineGrid[M]==T){
      *pt_i = i-1;
      *pt_j = M; 
    }else if (j==M+1 || T<0)
      printf("You chose a maturity that is not within the bounds\n");
    else{
      *pt_i = i-1;
      *pt_j = j-1;
    }
  }
}

double computeImpsig(double S_0, double r, double q, int optionType, double t_0, double K, double T, double V, double tol){

  double sigma_min,sigma_max;
  double sigma;
  double price;

  sigma_min = 0;
  sigma_max = 10;
  sigma = .5*(sigma_max+sigma_min);
  price = BSprice(optionType,S_0,t_0,sigma,r,q,T,K);
  while (fabs(price-V)>tol){
    if (price-V > 0)
      sigma_max = sigma;
    else
      sigma_min = sigma;
    sigma = .5*(sigma_max+sigma_min);
    price = BSprice(optionType,S_0,t_0,sigma,r,q,T,K);
  }

  return sigma;

}

/*computes the price of a european option (with sigma constant) with BS model*/
/* at time t, when the price of the asset is S, maturity=T, strike price=K,..etc...*/
double BSprice(int optionType, double S, double t, double sigma, double r, double q, double T, double K){

  double d,res;

  if (T == t){
    if (optionType==1){
      return (S-K>=0) ? S-K : 0;
    } else{
      return (K-S>=0) ? K-S : 0; 
    } 
  } else{
    d = (log(S/K) + (r-q+pow(sigma,2)/2)*(T-t))/(sigma*sqrt(T-t));
    res = S*exp(-q*(T-t))*distrib(d) - K*exp(-r*(T-t))*distrib(d-sigma*sqrt(T-t));
    if (optionType==1){
      return res;
    } else{
      res = res - S*exp(-q*(T-t))+K*exp(-r*(T-t));
      return res;
    }
  }
}



/*computes the price of a european call option (with sigma time-dependent) with BS model*/
/* at time t, when the price of the asset is S, maturity=T, strike price=K,..etc...     */
double BSpriceTimeDepVol(double S, double t, double r, double q, double T, double K, double (*intSigma)(double,double)){

  double d,sol,sum;

  if (T == t){
    return (S-K>=0) ? S-K : 0;
  } else{
    sum = intSigma(t,T);
    d = (log(S/K) + (r-q+sum/2)*(T-t))/(sqrt(sum*(T-t)));
    sol = S*exp(-q*(T-t))*distrib(d) - K*exp(-r*(T-t))*distrib(d-sqrt(sum*(T-t)));
    return sol;
  }
}



/* approximates the distribution function of a random variable ~ N(0,1) */
double distrib(double d){

  double p,b1,b2,b3,b4,b5,t,pi,res;
 
  pi = 3.14159;
  p = 0.2316419;
  b1 = 0.319381530;
  b2 = -0.356563782;
  b3 = 1.781477937;
  b4 = -1.821255978;
  b5 = 1.330274429;
  t = 1/(1+p*d);


  if (d>=0)
    res = 1 - (1/sqrt(2*pi))*exp(-pow(d,2)/2)*(b1*t + b2*pow(t,2) + b3*pow(t,3) + b4*pow(t,4) + b5*pow(t,5));
  else
    res = 1 - distrib(-d);
  
  return res;
}



/* computes the integral of sigma^2 over [t,T], divided by (T-t), when sigma is time-dependent */
/* in this case (for a test), we suppose sigma*(t)=t/2 */
double intSigma(double t, double T){
  //sigma = 1.0 - (t / 2.);
  //return 1 - 0.5 * (T+t) + 1/12. * (pow(T,2)+T*t+pow(t,2));
  //sigma = t / 2.;
  return (pow(T,3)-pow(t,3))/(12*(T-t));
  //sigma = t;
  //return (pow(T,3)-pow(t,3))/(3*(T-t));
}


/* computes a pseudo-random number between -1 and 1 */
double rand_nb(){

  int a,b;
  double s;
  time_t t;

  t=time(NULL);
  a=(int) t;
  srand(a);
  b=rand();
  s=2*((double)b)/((double)RAND_MAX)-1;
  return s;

}


/* simulates the price of the asset at time t, when the price at t_0 is S_0 */
double priceSimul(double t, double t_0, double S_0){

  double g1,g2,g3,res;

  g1=rand_nb();
  g2=rand_nb();
  g3=rand_nb();
  res = 1/sqrt((t-t_0)*(pow(g1+S_0/sqrt(t),2)+pow(g2,2)+pow(g3,2)));
  return res;

}


/* simulates nbSimul times the price of the asset and computes the mean value of the payoffs of the option*/
double priceEstim(double t, double t_0, double S_0, double strike, int nbSimul){

  double sum,assetPrice,optionPayoff;
  int i;

  sum=0;
  for(i=0;i<nbSimul;i++){
    assetPrice = priceSimul(t,t_0,S_0);
    optionPayoff = (assetPrice - strike >= 0) ? assetPrice-strike : 0;
    sum = sum + optionPayoff;
  }
  
  return sum/nbSimul;
 
}

