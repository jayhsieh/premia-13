#include <math.h>
#include "spline.h"
#include "solveSystem.h"
#include <stdlib.h>
#include <stdio.h>

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.
*/

double sigma_func(double S, double t){
/* OUTPUT
   - returns sigma(S,t)
   INPUTS
   - S
   - t   */


  
  //if (0<=t && t<=0.25)
  //  return 15/S;
  //else if (0.25<t && t<=0.5)
  //  return 40/S;
  //else if (0.5<t && t<=1)
  //  return 20/S;
  //else
  //  return -1;

  return 0.01+0.1*exp(-S/100)+0.01*t;
  //return .2;
  //return 1.0 - (t / 2.);
  //return t / 2.;
  //return t;
  //return 15 / S;
  //return 0.05+0.1*exp(-S/100)+0.5*t;
}


void buildCoarseGrids(double *y_coarseGrid, double *T_coarseGrid,int n, int m, double t_0, double T_max, double y_min, double y_max){
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

  double k,h;
  int i,j;

  /* steps for the coarse grid : */
  k = (T_max-t_0)/m; /* size of time step*/
  h = (y_max-y_min)/n; /* size of space step */

 /******************************initialization of y_coarseGrid and t_coarseGrid***********************/
  for (i=0;i<n+1;i++)
    y_coarseGrid[i] = y_min+h*i;
 
  for (i=0;i<m+1;i++)
    T_coarseGrid[i] = t_0+k*i;
  
}


void buildFineGrid(double *y_fineGrid, double *T_fineGrid, int N, int M, double t_0, double T_max, double y_min, double y_max, double S_0, int gridType){
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

  double K,H,valInvtanh;
  int i,l;
  struct marketData *pt_data;


  K = (T_max-t_0)/M; /* size of time step*/
  for (i=0;i<M+1;i++)
      T_fineGrid[i] = t_0+K*i;
   
  /* steps for the fine grid : */
  H = (y_max-y_min)/N;   
  
  for (i=0;i<N+1;i++)
    y_fineGrid[i] = y_min+H*i;
   
  
  if (gridType){ /* tanh grid */
    i=1;
    valInvtanh = invtanh(y_min+H,y_max,log(S_0));
    while (i<N && valInvtanh <= y_min){
      i = i+1;
      valInvtanh = invtanh(y_min+i*H,y_max,log(S_0)); 
    }
    
    /* the i  first values (from 0 to i-1) of invtanh(y_min+i*H,y_max,log(S_0))    */
    /* are below y_min so for the i first values of y_tanhGrid[.], we do a uniform */
    /* discretization so that y_coarseGrid[0] = y_min                              */
    for (l=i-1;l>=0;l--)
      y_fineGrid[l] =  valInvtanh - (i-l)*(valInvtanh - y_min)/i;
    
    while (i<N && valInvtanh < y_max){
      y_fineGrid[i] = valInvtanh;
      i = i+1;
      valInvtanh = invtanh(y_min+i*H,y_max,log(S_0));
    }
    /* once the image of y_min+i*H is greater than y_max, we do a uniform */ 
    /*   discretization until N, so that y_coarseGrid[N] = y_max          */
    for (l=i;l<=N;l++)
      y_fineGrid[l] =  y_fineGrid[i-1] + (l-i+1)*(y_max -  y_fineGrid[i-1])/(N-i+1);
    
    
    
  }  
  
}


double invtanh(double x, double a, double b){
  /*fonction invtanh (inverse de a*tanh(y-b)) */
  /* returns y so that a*tanh(y-b) = x        */
  return b+.5*log((a+x)/(a-x));
}




/****************************************************************************************************************************/
/**                                               DERIVATIVESGRIDS                                                         **/
/****************************************************************************************************************************/



void derivativesGrids(double **deriv_y, double **deriv_T, double **deriv_yT, double **sigmaCoarseGrid, double *y_coarseGrid, double *T_coarseGrid, struct derivData *data){
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


  int i,j,n,m;
  struct tridiag *V,*U; //matrices of the tridiagonal systems
  struct bidiagSystem *B,*C; //bidiagonal systems
  double aux;
  double *sol1,*sol2;/* will temporary contain the solutions of the bidiag systems at each step*/


  
  n = data->n;
  m = data->m;


  /*******************************************************************/
  /*       memory allocation                                         */
  /*******************************************************************/

  /*memory allocation for sol1 and sol2*/
  sol1 = (double *) malloc((n+1)*sizeof(double));
  sol2 = (double *) malloc((m+1)*sizeof(double));


  /* memory allocation for the tridiagonal matrix V */
  V = (struct tridiag *) malloc(sizeof(struct tridiag));
  V->subdiag = (double *) malloc(n*sizeof(double));
  V->diag = (double *) malloc((n+1)*sizeof(double));
  V->updiag = (double *) malloc(n*sizeof(double));
  V->size = n+1;


  /* memory allocation for the bidiagonal system B*/
  B = (struct bidiagSystem *) malloc(sizeof(struct bidiagSystem));
  B->T = (struct bidiag *) malloc(sizeof(struct bidiag));
  B->T->subdiag = (double *) malloc(n*sizeof(double));
  B->T->diag = (double *) malloc((n+1)*sizeof(double));
  B->T->size = n+1;
  B->b = (double *) malloc((n+1)*sizeof(double));
  B->size = n+1;
  
  /* memory allocation for the tridiagonal matrix U */
  U = (struct tridiag *) malloc(sizeof(struct tridiag));
  U->subdiag = (double *) malloc(m*sizeof(double));
  U->diag = (double *) malloc((m+1)*sizeof(double));
  U->updiag = (double *) malloc(m*sizeof(double));
  U->size = m+1;


  /* memory allocation for the bidiagonal system C*/
  C = (struct bidiagSystem *) malloc(sizeof(struct bidiagSystem));
  C->T = (struct bidiag *) malloc(sizeof(struct bidiag));
  C->T->subdiag = (double *) malloc(m*sizeof(double));
  C->T->diag = (double *) malloc((m+1)*sizeof(double));
  C->T->size = m+1;
  C->b = (double *) malloc((m+1)*sizeof(double));
  C->size = m+1;



  /*******************************************************************************/
  /*          derivatives wrt y                                                  */
  /*******************************************************************************/

  /* initialisation of V */
  V->diag[0] = 1;
  V->updiag[0] = 0;
  for (i=1;i<n;i++){
    V->updiag[i] = y_coarseGrid[i]-y_coarseGrid[i-1];
    V->diag[i] = 2*(y_coarseGrid[i+1]-y_coarseGrid[i-1]);
    V->subdiag[i-1] = y_coarseGrid[i+1]-y_coarseGrid[i];
  }
  V->subdiag[n-1] = 0;
  V->diag[n] = 1;


  /*builds the bidiagonal matrix B->T from the tridiagonal matrix V*/
  tridiagToBidiagMat(B->T,V);

  for (j=0;j<m+1;j++){
    /*computation of the right-hand side vector B->b of the bidiagonal system*/
    B->b[n] = data->deriv_y_n[j];
    B->b[0] = data->deriv_y_0[j];
    for (i=n-1;i>0;i--){
      aux = 3*((sigmaCoarseGrid[i+1][j]-sigmaCoarseGrid[i][j])*(y_coarseGrid[i]-y_coarseGrid[i-1])/(y_coarseGrid[i+1]-y_coarseGrid[i]) + (sigmaCoarseGrid[i][j]-sigmaCoarseGrid[i-1][j])*(y_coarseGrid[i+1]-y_coarseGrid[i])/(y_coarseGrid[i]-y_coarseGrid[i-1]));
      /*aux = ith value of the right-hand side vector of the tridiagonal system*/
      B->b[i] = aux -  V->updiag[i]*B->b[i+1]/(B->T->diag[i+1]);
    }

    /*solves the bidiag system B and stocks the result in deriv_y */
    solveSyst(sol1,B);
    
    for (i=0;i<=n;i++)
      deriv_y[i][j] = sol1[i];

  }

  /*******************************************************************************/
  /*          derivatives wrt t                                                  */
  /*******************************************************************************/
  
 
  /* initialization of U */
  U->diag[0] = 1;
  U->updiag[0] = 0;
  for (j=1;j<m;j++){
    U->updiag[j] = T_coarseGrid[j]-T_coarseGrid[j-1];
    U->diag[j] = 2*(T_coarseGrid[j+1]-T_coarseGrid[j-1]);
    U->subdiag[j-1] = T_coarseGrid[j+1]-T_coarseGrid[j];
  }
  U->subdiag[m-1] = 0;
  U->diag[m] = 1;

  /*build the bidiagonal matrix C from U*/
  tridiagToBidiagMat(C->T,U);
 
 

  for (i=0;i<=n;i++){
    
    /*computation of the right-hand side vector C->b  of the bidiagonal system*/
    C->b[m] = data->deriv_T_m[i];
    C->b[0] = data->deriv_T_0[i];
   
  for (j=m-1;j>0;j--){
      aux = 3*((sigmaCoarseGrid[i][j+1]-sigmaCoarseGrid[i][j])*(T_coarseGrid[j]-T_coarseGrid[j-1])/(T_coarseGrid[j+1]-T_coarseGrid[j]) + (sigmaCoarseGrid[i][j]-sigmaCoarseGrid[i][j-1])*(T_coarseGrid[j+1]-T_coarseGrid[j])/(T_coarseGrid[j]-T_coarseGrid[j-1]));
      /*aux = jth value of the right-hand side vector of the tridiagonal system corresponding to C*/
      C->b[j] = aux -  U->updiag[j]*C->b[j+1]/(C->T->diag[j+1]);

  }


  /*solves the bidiag system C and stocks the result in deriv_T */
  solveSyst(sol2,C);



  for (j=0;j<=m;j++)
    deriv_T[i][j] = sol2[j];
  
  } 


  /*******************************************************************************/
  /*          cross derivatives for j=0,m and i=1...n-1                          */
  /*******************************************************************************/ 

  /*the tridiagonal matrix of the system is V, the bidiagonal matrix is then B->T */
  /*we keep the bidiag system B, we will only have to change the right-hand side  */

  /*j=0*/
  /*computation of the right-hand side vector B->b of the bidiagonal system*/
  B->b[n] = data->deriv_yT_n0;
  B->b[0] = data->deriv_yT_00;
  for (i=n-1;i>0;i--){
    aux = 3*((deriv_T[i+1][0]-deriv_T[i][0])*(y_coarseGrid[i]-y_coarseGrid[i-1])/(y_coarseGrid[i+1]-y_coarseGrid[i])+(deriv_T[i][0]-deriv_T[i-1][0])*(y_coarseGrid[i+1]-y_coarseGrid[i])/(y_coarseGrid[i]-y_coarseGrid[i-1]));
    /*aux = ith value of the right-hand side vector of the tridiagonal system*/
    B->b[i] = aux -  V->updiag[i]*B->b[i+1]/(B->T->diag[i+1]);
  }

  /*solves the bidiag system B and stocks the result in deriv_y */
  solveSyst(sol1,B);
  for (i=0;i<=n;i++)
    deriv_yT[i][0] = sol1[i];

  /*j=m*/
  /*computation of the right-hand side vector B->b of the bidiagonal system*/
  B->b[n] = data->deriv_yT_nm;
  B->b[0] = data->deriv_yT_0m;
  for (i=n-1;i>0;i--){
    aux = 3*((deriv_T[i+1][m]-deriv_T[i][m])*(y_coarseGrid[i]-y_coarseGrid[i-1])/(y_coarseGrid[i+1]-y_coarseGrid[i])+(deriv_T[i][m]-deriv_T[i-1][m])*(y_coarseGrid[i+1]-y_coarseGrid[i])/(y_coarseGrid[i]-y_coarseGrid[i-1]));
    /*aux = ith value of the right-hand side vector of the tridiagonal system*/
    B->b[i] = aux -  V->updiag[i]*B->b[i+1]/(B->T->diag[i+1]);
  }
  
  /*solves the bidiag system B and stocks the result in deriv_yT */
  solveSyst(sol1,B);
  for (i=0;i<=n;i++)
    deriv_yT[i][m] = sol1[i];
   



  /*******************************************************************************/
  /*          cross derivatives for j=1,...,M-1 and i=0...n                      */
  /*******************************************************************************/ 
  
  /*the tridiagonal matrix of the system is U, the bidiagonal matrix is then C->T */
  /*we keep the bidiag system C, we will only have to change the right-hand side  */

  
  for (i=0;i<=n;i++){
    
    /*computation of the right-hand side vector C->b of the bidiagonal system*/
    C->b[m] = deriv_yT[i][m];
    C->b[0] = deriv_yT[i][0];


  for (j=m-1;j>0;j--){
      aux = 3*((deriv_y[i][j+1]-deriv_y[i][j])*(T_coarseGrid[j]-T_coarseGrid[j-1])/(T_coarseGrid[j+1]-T_coarseGrid[j])+(deriv_y[i][j]-deriv_y[i][j-1])*(T_coarseGrid[j+1]-T_coarseGrid[j])/(T_coarseGrid[j]-T_coarseGrid[j-1]));
      /*aux = jth value of the right-hand side vector of the tridiagonal system corresponding to C*/
      C->b[j] = aux -  U->updiag[j]*C->b[j+1]/(C->T->diag[j+1]);

  }

  
  /*solves the bidiag system C and stocks the result in deriv_yT */
  solveSyst(sol2,C);
  for (j=0;j<=m;j++)
    deriv_yT[i][j] = sol2[j];
 
  } 


  /***********************************************************************/
  /*                       memory desallocation                          */
  /***********************************************************************/
  
  free(V->subdiag);
  V->subdiag = NULL;
  free(V->diag);
  V->diag = NULL;
  free(V->updiag); 
  V->updiag = NULL;
  free(V);
  V = NULL;

  free(U->subdiag);
  U->subdiag = NULL;
  free(U->diag);
  U->diag = NULL;
  free(U->updiag); 
  U->updiag = NULL;
  free(U);
  U = NULL;
 
  free(B->T->subdiag);
  B->T->subdiag = NULL;
  free(B->T->diag);
  B->T->diag = NULL;
  free(B->T);
  B->T = NULL;
  free(B->b);
  B->b = NULL;
  free(B);
  B = NULL;

  free(C->T->subdiag);
  C->T->subdiag = NULL;
  free(C->T->diag);
  C->T->diag = NULL;
  free(C->T);
  C->T = NULL;
  free(C->b);
  C->b = NULL;
  free(C);
  C = NULL;

  free(sol1);
  sol1 = NULL;
  free(sol2);
  sol2 = NULL;

}
  

/**********************************************************************************************************************************************/
/**                                                         INTERPCOEFFS                                                                     **/
/**********************************************************************************************************************************************/

  
void interpCoeffs(double ****coeffs, double *y_coarseGrid, double *T_coarseGrid, double **deriv_y, double **deriv_T, double **deriv_yT, double **sigmaCoarseGrid, int n, int m){
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



  int i,j,line,col;
  double **L; /* K_ij . A(delta_t)' */
  double **K; /* matrix of given values (cf report) */
  double **H; /* 2 last columns of  A(delta_t)' */
  double **P; /* 2 last lines of  A(delta_y) */
  double **KH; /* KH=K.H */
  double **PL; /* PL=P.L */
  double delta_t, delta_y;


  /*memory allocation for K,L,H,P,KH,PL) */
  L = (double **) malloc(4*sizeof(double *));
  for (line=0;line<4;line++)
    L[line] = (double *) malloc(4*sizeof(double));

  K = (double **) malloc(4*sizeof(double *));
  for (line=0;line<4;line++)
    K[line] = (double *) malloc(4*sizeof(double));
 
  H = (double **) malloc(4*sizeof(double *));
  for (line=0;line<4;line++)
    H[line] = (double *) malloc(2*sizeof(double));

  P = (double **) malloc(2*sizeof(double *));
  for (line=0;line<2;line++)
    P[line] = (double *) malloc(4*sizeof(double));

  KH = (double **) malloc(4*sizeof(double *));
  for (line=0;line<4;line++)
    KH[line] = (double *) malloc(2*sizeof(double));

  PL = (double **) malloc(2*sizeof(double *));
  for (line=0;line<2;line++)
    PL[line] = (double *) malloc(4*sizeof(double));




  for (i=1;i<=n;i++){
    for (j=1;j<=m;j++){

      delta_t = T_coarseGrid[j]-T_coarseGrid[j-1];
      delta_y = y_coarseGrid[i]-y_coarseGrid[i-1];

      /* computation of K */
      computeK(K,i,j,deriv_T,deriv_y,deriv_yT,sigmaCoarseGrid);

      /* computation of H */
      computeH(H,delta_t);

      /*computation of KH */
      matMult(KH,K,H,4,4,2);

      /*computation of L : the two first columns of L are the same than the */
      /* two first columns of K. The 2 last columns of L are equal to KH    */
      for (line=0;line<4;line++)
	for (col=0;col<2;col++)
	  L[line][col] = K[line][col];
      
      for (line=0;line<4;line++)
	for (col=2;col<4;col++)
	  L[line][col] = KH[line][col-2];

      
      /* computation of P */
      computeP(P,delta_y);

      /* computation of PL */
      matMult(PL,P,L,4,2,4);

      /*computation of coeff_ij : the two first lines are the same than the */
      /* two first lines of L. The 2 last columns are equal to PL           */
      for (line=0;line<2;line++)
	for (col=0;col<4;col++){
	  coeffs[i-1][j-1][line][col] = L[line][col];
	  //printf("%lf\n",coeffs[i-1][j-1][line][col]);
	}
      
      for (line=2;line<4;line++)
	for (col=0;col<4;col++){
	   coeffs[i-1][j-1][line][col] = PL[line-2][col];
	   //printf("%lf\n",coeffs[i-1][j-1][line][col]);
	}


    
    }
    
  }



/* memory desallocation */

  for (line=0;line<4;line++){
    free(L[line]);
    L[line] = NULL;
  }
  free(L);
  L = NULL;

  for (line=0;line<4;line++){
    free(K[line]);
    K[line] = NULL;
  }
  free(K);
  K = NULL;

  for (line=0;line<4;line++){
    free(H[line]);
    H[line] = NULL;
  }
  free(H);
  H = NULL;

  for (line=0;line<2;line++){
    free(P[line]);
    P[line] = NULL;
  }
  free(P);
  P = NULL;

  for (line=0;line<4;line++){
    free(KH[line]);
    KH[line] = NULL;
  }
  free(KH);
  KH = NULL;

  for (line=0;line<2;line++){
    free(PL[line]);
    PL[line] = NULL;
  }
  free(PL);
  PL = NULL;

}





void computeK(double **K, int i, int j, double **deriv_T, double **deriv_y, double **deriv_yT, double **sigmaCoarseGrid){
/*OUTPUT
  - K (=Kij) =  matrix of size (4,4)
  INPUT
  - i,j = indices of the matrix Kij to be computed
  - deriv_y = matrix of the derivatives wrt y                                                 
  - deriv_T = matrix of the derivatives wrt T                                                 
  - deriv_yT =  matrix of the cross derivatives wrt y and T 
  - sigmaCoarseGrid = coarse grid of the values of sigma */

  K[0][0] = sigmaCoarseGrid[i-1][j-1];
  K[0][1] = deriv_T[i-1][j-1]; 
  K[0][2] = sigmaCoarseGrid[i-1][j];
  K[0][3] = deriv_T[i-1][j];
  K[1][0] = deriv_y[i-1][j-1];
  K[1][1] = deriv_yT[i-1][j-1]; 
  K[1][2] = deriv_y[i-1][j];
  K[1][3] = deriv_yT[i-1][j];
  K[2][0] = sigmaCoarseGrid[i][j-1];
  K[2][1] = deriv_T[i][j-1]; 
  K[2][2] = sigmaCoarseGrid[i][j];
  K[2][3] = deriv_T[i][j];
  K[3][0] = deriv_y[i][j-1];
  K[3][1] = deriv_yT[i][j-1]; 
  K[3][2] = deriv_y[i][j];
  K[3][3] = deriv_yT[i][j];

}



void computeH(double **H, double delta_T){
/*OUTPUT
  - H =  matrix of size (4,4)
  INPUT
  - delta_T = time step */



  double delta_T_2,delta_T_3;
  
  delta_T_2 = pow(delta_T,2);
  delta_T_3 = pow(delta_T,3);

  H[0][0] = -3/delta_T_2;
  H[1][0] = -2/delta_T;
  H[2][0] = 3/delta_T_2;
  H[3][0] = -1/delta_T;
  H[0][1] = 2/delta_T_3;
  H[1][1] = 1/delta_T_2;
  H[2][1] = -2/delta_T_3;
  H[3][1] = 1/delta_T_2;
  
}

void computeP(double **P, double delta_y){
/*OUTPUT
  - P =  matirx of size (4,4)
  INPUT
  - delta_y = space step */

  double delta_y_2,delta_y_3;
  
  delta_y_2 = pow(delta_y,2);
  delta_y_3 = pow(delta_y,3);

  P[0][0] = -3/delta_y_2;
  P[0][1] = -2/delta_y;
  P[0][2] = 3/delta_y_2;
  P[0][3] = -1/delta_y;
  P[1][0] = 2/delta_y_3;
  P[1][1] = 1/delta_y_2;
  P[1][2] = -2/delta_y_3;
  P[1][3] = 1/delta_y_2;
  
}



void matMult(double **AB,double **A, double **B, int n, int m, int p){
/* OUPUT:
   - AB = A*B of size (m,p)
   INPUTS
   - A = matrix of size (m,n)
   - B = matrix of size (n,p) */


  int i,j,k;
  double sum;

  for (i=0;i<m;i++)
    for (j=0;j<p;j++){
      sum = 0;
      for (k=0;k<n;k++)
	sum = sum + A[i][k]*B[k][j];
      AB[i][j] = sum;
    }

}


int find_index(double x, double *grid, int size, int start){
/* OUTPUT
   - returns the index i such that grid[i]<=x<grid[i+1]
   INPUT
   - x (x = y or T)
   - grid = grid of discretized values of x, grid = x_0,...,x_n (x = y or T)
   - size = size of grid - 1 (M for T_grid, N for y_grid)
   - start = index at which we start the search of i */

  int k;

  if (grid[size] <= x)
    return size-1;
  else{
    k=start;
    while (k<size && grid[k]<=x)
      k++;
    return k-1;
  }
}




double interpEval(double ****coeffs, double y, double T, double *y_coarseGrid, double *T_coarseGrid, int n, int m){
/* OUTPUT
   - returns the evaluation of the interpolation function at the point (y,T)
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_coarseGrid = grid of discretized values of y, grid_y = y_0,...,y_n
   - T_coarseGrid = grid of discretized values of T, grid_T = T_0,...,T_m
   - n = size of y_coarseGrid - 1 
   - m : size of T_coarseGrid - 1 */ 

  int l,k,i,j;
  double sum;

  i = find_index(y,y_coarseGrid,n,0);
  j = find_index(T,T_coarseGrid,m,0);

  sum = 0;
  for (k=0;k<4;k++)
    for (l=0;l<4;l++)
      sum = sum + coeffs[i][j][k][l]*pow(y-y_coarseGrid[i],k)*pow(T-T_coarseGrid[j],l);
 
  return sum;

}
  
double interpEval2(double ****coeffs, double y, double T, int i, int j, double y_i, double T_j){
/* OUTPUT
   - returns the evaluation of the interpolation function at the point (y,T)
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_i = y_coarseGrid[i]
   - T_j = T_coarseGrid[j]          */

  int l,k;
  double sum;

  sum = 0;
  for (k=0;k<4;k++)
    for (l=0;l<4;l++)
      sum = sum + coeffs[i][j][k][l]*pow(y-y_i,k)*pow(T-T_j,l);

  return sum;

}

double interpEval2_dy(double ****coeffs, double y, double T, int i, int j, double y_i, double T_j){
/* OUTPUT
   - returns the evaluation of the derivative wrt to y at the point (y,T) 
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_i = y_coarseGrid[i]
   - T_j = T_coarseGrid[j]          */

  int l,k;
  double sum;

  sum = 0;
  for (k=1;k<4;k++)
    for (l=0;l<4;l++)
      sum = sum + k*coeffs[i][j][k][l]*pow(y-y_i,k-1)*pow(T-T_j,l);

  return sum;

}

double interpEval2_dT(double ****coeffs, double y, double T, int i, int j, double y_i, double T_j){
/* OUTPUT
   - returns the evaluation of the derivative wrt to T at the point (y,T) 
   INPUTS
   - coeffs = 4-dim array containing the coefficients of the interp. function  
   - (y,T) = point at which the interpolation function will be evaluated
   - y_i = y_coarseGrid[i]
   - T_j = T_coarseGrid[j]          */

  int l,k;
  double sum;

  sum = 0;
  for (k=0;k<4;k++)
    for (l=1;l<4;l++)
      sum = sum + l*coeffs[i][j][k][l]*pow(y-y_i,k)*pow(T-T_j,l-1);

  return sum;

}


void discretizeSigma(double (*volatility)(double,double), double **sigmaCoarseGrid, int n, int m, double *y_coarseGrid, double *T_coarseGrid){
/* OUTPUT
   - sigmaCoarseGrid = coarse grid of discretized values of sigma
   INPUTS
   - volatility = function defining the volatility
   - n,m = size of the coarse grid
   - y_coarseGrid = array of dicretized values y (for the coarse grid)
   - T_coarseGrid = array of dicretized values T (for the coarse grid) */
   
  int i,j;
  double y_min,y_max,t_0,T_max;
  
  y_min = y_coarseGrid[0];
  y_max = y_coarseGrid[n];
  t_0 = T_coarseGrid[0];
  T_max = T_coarseGrid[m];
  


  for (i=0;i<=n;i++)
    for(j=0;j<=m;j++)
      sigmaCoarseGrid[i][j] = volatility(exp(y_coarseGrid[i]),T_coarseGrid[j]);

}


void interpole(double **sigmaFineGrid, double **sigmaCoarseGrid, int n, int m, int N, int M, struct derivData *data, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid){
/*sigmaFineGrid : N*M */
/*sigmaCoarseGrid : n*m , n<<N, m<<M*/

  int i,j,l,I,J;
  double ****coeffs;
  double **deriv_y,**deriv_T,**deriv_yT;
  double y_min,y_max,t_0,T_max;
  
  y_min = y_coarseGrid[0];
  y_max = y_coarseGrid[n];
  t_0 = T_coarseGrid[0];
  T_max = T_coarseGrid[m];
  
  /* memory allocation for the 4-dim array coeffs*/ 
  coeffs = (double ****) malloc(n*sizeof(double ***));
  for (i=0;i<n;i++){
    coeffs[i] = (double ***) malloc(m*sizeof(double **));
    for (j=0;j<m;j++){
      coeffs[i][j] = (double **) malloc(4*sizeof(double *));
      for (l=0;l<4;l++)
	coeffs[i][j][l] = (double *) malloc(4*sizeof(double));
    }
  }

  /* memory allocation for the derivatives grids */
  deriv_y = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    deriv_y[i] = (double *) malloc((m+1)*sizeof(double));

  deriv_T = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    deriv_T[i] = (double *) malloc((m+1)*sizeof(double));

  deriv_yT = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    deriv_yT[i] = (double *) malloc((m+1)*sizeof(double));


  /* computation of the 1st order derivatives and cross derivatives at all the points of the coarse grid*/
  derivativesGrids(deriv_y,deriv_T,deriv_yT,sigmaCoarseGrid,y_coarseGrid,T_coarseGrid,data);

  /* computation of the coeffs of the interpolation function */
  interpCoeffs(coeffs,y_coarseGrid,T_coarseGrid,deriv_y,deriv_T,deriv_yT,sigmaCoarseGrid,n,m);

   /* for each point of the fine grid, we compute i and j so that sigmaFineGrid[I][J] is in R_ij */
  /* and we evaluate the interpolating function in this point                                */


  i = 0;
  for (I=0;I<=N;I++){
    i = find_index(y_fineGrid[I],y_coarseGrid,n,0);
    j = 0;
    for (J=0;J<=M;J++){
       j = find_index(T_fineGrid[J],T_coarseGrid,m,0);
       sigmaFineGrid[I][J] = interpEval2(coeffs,y_fineGrid[I],T_fineGrid[J],i,j,y_coarseGrid[i],T_coarseGrid[j]);
    }
  }

  /* " k_coarseGrid = exp(y_coarseGrid); " */
/*   i = 0; */
/*   for (I=0;I<=N;I++){ */
/*     i = find_index(exp(y_fineGrid[I]),k_coarseGrid,n,0); */
/*     j = 0; */
/*     for (J=0;J<=M;J++){ */
/*        j = find_index(T_fineGrid[J],T_coarseGrid,m,0); */
/*        sigmaFineGrid[I][J] = interpEval2(coeffs,exp(y_fineGrid[I]),T_fineGrid[J],i,j,k_coarseGrid[i],T_coarseGrid[j]); */
/*     } */
/*   } */


  for (i=0;i<n;i++){
    for (j=0;j<m;j++){
      for (l=0;l<4;l++){
	free(coeffs[i][j][l]);
	coeffs[i][j][l] = NULL;
      }
      free(coeffs[i][j]);
      coeffs[i][j] = NULL;
    }
    free(coeffs[i]);
    coeffs[i] = NULL;
  }
  free(coeffs);
  coeffs = NULL;

  for (i=0;i<=n;i++){
    free(deriv_y[i]);
    deriv_y[i] = NULL;
  }
  free(deriv_y);
  deriv_y = NULL;

  for (i=0;i<=n;i++){
    free(deriv_T[i]);
    deriv_T[i] = NULL;
  }
  free(deriv_T);
  deriv_T = NULL;

  for (i=0;i<=n;i++){
    free(deriv_yT[i]);
    deriv_yT[i] = NULL;
  }
  free(deriv_yT);
  deriv_yT = NULL;

}

void ddlGrid1ToGrid2(double **pt_ddl2, double *ddl1, int n1, int m1, int n2, int m2, double *y_coarseGrid1, double *T_coarseGrid1, double *y_coarseGrid2, double *T_coarseGrid2){
/*sigmaCoarseGrid2 : n2*m2 */
/*sigmaCoarseGrid1 : n1*m1 , n1<n2, m1<m2*/

  int i,j,k,l,I,J;
  double ****coeffs;
  double **deriv_y1,**deriv_T1,**deriv_yT1;
  double y_min,y_max,t_0,T_max;
  double **sigmaCoarseGrid1;
  struct derivData *data;

  (*pt_ddl2) = (double *) malloc((n2+3)*(m2+3)*sizeof(double));

  /* memory alloc for sigmaCoarseGrid1 */
  sigmaCoarseGrid1 = (double **) malloc((n1+1)*sizeof(double *));
  for (j=0;j<n1+1;j++)
    sigmaCoarseGrid1[j] = (double *) malloc((m1+1)*sizeof(double));

  data = (struct derivData *) malloc(sizeof(struct derivData));
  data->deriv_y_0 = (double *) malloc((m1+1)*sizeof(double));
  data->deriv_y_n = (double *) malloc((m1+1)*sizeof(double));
  data->deriv_T_0 = (double *) malloc((n1+1)*sizeof(double));
  data->deriv_T_m = (double *) malloc((n1+1)*sizeof(double));


  /* memory allocation for the 4-dim array coeffs*/ 
  coeffs = (double ****) malloc(n1*sizeof(double ***));
  for (i=0;i<n1;i++){
    coeffs[i] = (double ***) malloc(m1*sizeof(double **));
    for (j=0;j<m1;j++){
      coeffs[i][j] = (double **) malloc(4*sizeof(double *));
      for (l=0;l<4;l++)
	coeffs[i][j][l] = (double *) malloc(4*sizeof(double));
    }
  }

  /* memory allocation for the derivatives grids */
  deriv_y1 = (double **) malloc((n1+1)*sizeof(double *));
  for (i=0;i<n1+1;i++)
    deriv_y1[i] = (double *) malloc((m1+1)*sizeof(double));

  deriv_T1 = (double **) malloc((n1+1)*sizeof(double *));
  for (i=0;i<n1+1;i++)
    deriv_T1[i] = (double *) malloc((m1+1)*sizeof(double));

  deriv_yT1 = (double **) malloc((n1+1)*sizeof(double *));
  for (i=0;i<n1+1;i++)
    deriv_yT1[i] = (double *) malloc((m1+1)*sizeof(double));



  /* ddl1 => sigmaCoarseGrid1, data */
  for (k=0;k<m1+1;k++)
    for (l=0;l<n1+1;l++)
      sigmaCoarseGrid1[l][k] = ddl1[k*(n1+1)+l];
  
  for (k=(n1+1)*(m1+1);k<(n1+1)*(m1+1)+m1+1;k++)
    data->deriv_y_0[k-(n1+1)*(m1+1)] = ddl1[k];

  for (k=(n1+2)*(m1+1);k<(n1+2)*(m1+1)+m1+1;k++)
    data->deriv_y_n[k-(n1+2)*(m1+1)] = ddl1[k];
  
  for (k=(n1+3)*(m1+1);k<(n1+3)*(m1+1)+n1+1;k++)
    data->deriv_T_0[k-(n1+3)*(m1+1)] = ddl1[k];

  for (k=(n1+3)*(m1+1)+n1+1;k<(n1+3)*(m1+1)+2*(n1+1);k++)
    data->deriv_T_m[k-((n1+3)*(m1+1)+n1+1)] = ddl1[k];
  
  data->deriv_yT_00 = ddl1[(n1+3)*(m1+1)+2*(n1+1)];
  data->deriv_yT_n0 = ddl1[(n1+3)*(m1+1)+2*(n1+1)+1];
  data->deriv_yT_0m = ddl1[(n1+3)*(m1+1)+2*(n1+1)+2];
  data->deriv_yT_nm = ddl1[(n1+3)*(m1+1)+2*(n1+1)+3];
  
  data->n = n1;
  data->m = m1;

  /* computation of the 1st order derivatives and cross derivatives at all the points of the coarse grid*/
  derivativesGrids(deriv_y1,deriv_T1,deriv_yT1,sigmaCoarseGrid1,y_coarseGrid1,T_coarseGrid1,data);
  
  /* computation of the coeffs of the interpolation function */
  interpCoeffs(coeffs,y_coarseGrid1,T_coarseGrid1,deriv_y1,deriv_T1,deriv_yT1,sigmaCoarseGrid1,n1,m1);

  /* for each point of the fine grid, we compute i and j so that sigmaCoarseGrid2[I][J] is in R_ij */
  /* and we evaluate the interpolating function in this point */

  /* d/dy */
  for (J=0;J<m2;J++){
    j = find_index(T_coarseGrid2[J],T_coarseGrid1,m1,0);
    (*pt_ddl2)[(m2+1)*(n2+1)+J] = interpEval2_dy(coeffs,y_coarseGrid2[0],T_coarseGrid2[J],0,j,y_coarseGrid1[0],T_coarseGrid1[j]); /* derivatives wrt y for i=0 */
    (*pt_ddl2)[(m2+1)*(n2+1)+(m2+1)+J] = interpEval2_dy(coeffs,y_coarseGrid2[n2],T_coarseGrid2[J],n1-1,j,y_coarseGrid1[n1-1],T_coarseGrid1[j]); /* derivatives wrt y for i=n2 */
  }
  (*pt_ddl2)[(m2+1)*(n2+1)+m2] = deriv_y1[0][m1];   
  (*pt_ddl2)[(m2+1)*(n2+1)+(m2+1)+m2] = deriv_y1[n1][m1];

  /* d/dT */
  for (I=0;I<n2;I++){
    i = find_index(y_coarseGrid2[I],y_coarseGrid1,n1,0);
    (*pt_ddl2)[(m2+1)*(n2+1)+2*(m2+1)+I] = interpEval2_dT(coeffs,y_coarseGrid2[I],T_coarseGrid2[0],i,0,y_coarseGrid1[i],T_coarseGrid1[0]); /* derivatives wrt T for j=0 */
    (*pt_ddl2)[(m2+1)*(n2+1)+2*(m2+1)+(n2+1)+I] = interpEval2_dT(coeffs,y_coarseGrid2[I],T_coarseGrid2[m2],i,m1-1,y_coarseGrid1[i],T_coarseGrid1[m1-1]); /* derivatives wrt T for j=m2 */
  }
  (*pt_ddl2)[(m2+1)*(n2+1)+2*(m2+1)+n2] = deriv_T1[n1][0];
  (*pt_ddl2)[(m2+1)*(n2+1)+2*(m2+1)+(n2+1)+n2] = deriv_T1[n1][m1];   

  /* d^2/dydT */ 
  for (i=0;i<4;i++)
    (*pt_ddl2)[(n2+3)*(m2+3)-4+i] = ddl1[(n1+3)*(m1+3)-4+i];

  /* function */
  i = 0;
  for (I=0;I<=n2;I++){
    i = find_index(y_coarseGrid2[I],y_coarseGrid1,n1,i);
    j = 0;
    for (J=0;J<=m2;J++){
      j = find_index(T_coarseGrid2[J],T_coarseGrid1,m1,j);
      (*pt_ddl2)[J*(n2+1)+I] = interpEval2(coeffs,y_coarseGrid2[I],T_coarseGrid2[J],i,j,y_coarseGrid1[i],T_coarseGrid1[j]);
    }
  }
  
  /* free */
  for (i=0;i<n1;i++){
    for (j=0;j<m1;j++){
      for (l=0;l<4;l++){
	free(coeffs[i][j][l]);
	coeffs[i][j][l] = NULL;
      }
      free(coeffs[i][j]);
      coeffs[i][j] = NULL;
    }
    free(coeffs[i]);
    coeffs[i] = NULL;
  }
  free(coeffs);
  coeffs = NULL;
  
  for (i=0;i<=n1;i++){
    free(deriv_y1[i]);
    deriv_y1[i] = NULL;
  }
  free(deriv_y1);
  deriv_y1 = NULL;
  
  for (i=0;i<=n1;i++){
    free(deriv_T1[i]);
    deriv_T1[i] = NULL;
  }
  free(deriv_T1);
  deriv_T1 = NULL;
  
  for (i=0;i<=n1;i++){
    free(deriv_yT1[i]);
    deriv_yT1[i] = NULL;
  }
  free(deriv_yT1);
  deriv_yT1 = NULL;
  
}

void printDerivData(struct derivData *d){

  int i;

  printf("DerivData :\n");
  printf("n=%d\n",d->n);
  printf("m=%d\n",d->m);
  printf("deriv_y_0:\n");
  for (i=0;i<d->m+1;i++)
    printf("%lf ",d->deriv_y_0[i]);
  printf("\nderiv_y_n:\n");
  for (i=0;i<d->m+1;i++)
    printf("%lf ",d->deriv_y_n[i]);
  printf("\nderiv_T_0:\n");
  for (i=0;i<d->n+1;i++)
    printf("%lf ",d->deriv_T_0[i]);
  printf("\nderiv_T_m:\n");
  for (i=0;i<d->n+1;i++)
    printf("%lf ",d->deriv_T_m[i]);
  printf("\n deriv_yT_00 = %lf\n",d->deriv_yT_00);
  printf("deriv_yT_0m = %lf\n",d->deriv_yT_0m);
  printf("deriv_yT_n0 = %lf\n",d->deriv_yT_n0);
  printf("deriv_yT_nm = %lf\n",d->deriv_yT_nm);
  


}

