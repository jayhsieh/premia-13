#include "grad.h"
#include "DupirePDE.h"
#include "spline.h"
#include "utilsGrad.h"
#include "sparse.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

void computeGrad_F(double *grad_F, double *sigma_param, double lambda, int N, int M, int n, int m, double *y_fineGrid, double *T_fineGrid, double *y_coarseGrid, double *T_coarseGrid, double r, double q, double S_0, int optionType, struct marketData **marketOptionPrices, int nbMaturites){
/* computes the gradient of F
   OUTPUT:
   - grad_F (gradient evaluated in sigma_param)
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
   - nbMaturites : nb of maturites for which we know market prices 
*/

   
  double **invA_times_B;
  double **invC_times_D;
  int **D_indices;
  int **B_indices;
  struct tridiag *D;
  double **matGrad_F1;
  double *grad_G;
  int i,j,nbParamSigma;

  nbParamSigma = (n+1)*(m+1) + 2*(m+1) + 2*(n+1) + 4;

  grad_G = (double *)malloc(nbParamSigma*sizeof(double));

  /* memory allocation */
  D = (struct tridiag *) malloc(sizeof(struct tridiag));
  D->subdiag = (double *) malloc(m*sizeof(double));
  D->diag = (double *) malloc((m+1)*sizeof(double));
  D->updiag = (double *) malloc(m*sizeof(double));
  D->size = m+1;

  B_indices = (int **) malloc((m+1)*sizeof(int *));
  for (j=0;j<m+1;j++)
    B_indices[j] = (int *) malloc((n+3)*sizeof(int));

  D_indices = (int **) malloc((n+1)*sizeof(int *));
  for (i=0;i<n+1;i++)
    D_indices[i] = (int *) malloc((m+3)*sizeof(int));
  
  invA_times_B = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    invA_times_B[i] = (double *) malloc((n+3)*sizeof(double));
  
  invC_times_D = (double **) malloc((m+1)*sizeof(double *));
  for (j=0;j<m+1;j++)
    invC_times_D[j] = (double *) malloc((m+3)*sizeof(double));
  
  
  matGrad_F1 = (double **) malloc(nbParamSigma*sizeof(double *));
  for (i=0;i<nbParamSigma;i++)
    matGrad_F1[i] = (double *) malloc(nbParamSigma*sizeof(double));
  
  computeGrad_F1(matGrad_F1,invA_times_B,invC_times_D,D_indices,B_indices,D,n,m,y_coarseGrid,T_coarseGrid); /* on peut peut etre l'appeler a l'exterieur, a voir */

  computeGrad_G(grad_G,sigma_param,N,M,n,m,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid,r,q,optionType,S_0,invA_times_B,invC_times_D,D_indices,B_indices,D,marketOptionPrices,nbMaturites);

  /* grad_F = lambda*matGrad_F1*sigma_param + grad_G */

  
  matMultVect(grad_F,matGrad_F1,sigma_param,nbParamSigma,nbParamSigma); 
  
  for (i=0;i<nbParamSigma;i++)
    grad_F[i] = lambda*grad_F[i] + grad_G[i];

  
  /* memory desallocation */

  free(grad_G);

  free(D->subdiag);
  free(D->diag);
  free(D->updiag);
  free(D);

  for (j=0;j<m+1;j++)
    free(B_indices[j]);
  free(B_indices);

  for (i=0;i<n+1;i++)
    free(D_indices[i]);
  free(D_indices);

  for (j=0;j<n+1;j++)
    free(invA_times_B[j]);
  free(invA_times_B);

  for (i=0;i<m+1;i++)
    free(invC_times_D[i]);
  free(invC_times_D);

  for (j=0;j<nbParamSigma;j++)
    free(matGrad_F1[j]);
  free(matGrad_F1);

}



/***************************************************************************************************************************/
/*                             COMPUTATION OF GRAD(F1)                                                                     */
/***************************************************************************************************************************/


void computeGrad_F1(double **matGradF1, double **invA_times_B, double **invC_times_D, int **D_indices, int **B_indices, struct tridiag *D, int n, int m, double *y_coarseGrid, double *T_coarseGrid){
/* computes the matrix representing the gradient of F1 (then, grad_F1(sigma_param) = matGradF1*sigma_param 
   OUTPUT:
   - matGradF1 (matrix of the gradient)
   - invA_times_B (matrix we will need in the computation of gradF2)
   - invC_times_D (matrix we will need in the computation of gradF2) 
   - D_indices(matrix we will need in the computation of gradF2) 
   - B_indices(matrix we will need in the computation of gradF2) 
   - D (tridiag matrix we will need in the computation of gradF2) 
   INPUTS:
   - y_coarseGrid : discretized values of y  (for the coarse grid)                     
   - T_coarseGrid : discretized values of T  (for the coarse grid)        
   - (n,m) = size of the coarse grid
*/

  
  /* declarations */
  int i,j,k,line,col,nbParamSigma;
  double delta1,delta2;
  struct tridiag *A,*B,*C;
  double **invA,**invC,**tempMat,**tempMat2;
  struct sparseMat *N_ij;
  struct sparseMat *M_ij;
  struct cell *pt_cell,*pt_line_Nij;
  struct cell **above_Nij,**sentinelles_lines,**sentinelles_columns;


  nbParamSigma = (n+1)*(m+1)+2*(n+1)+2*(m+1)+4;

  /* memory allocation */

  A = (struct tridiag *) malloc(sizeof(struct tridiag));
  A->subdiag = (double *) malloc(n*sizeof(double));
  A->diag = (double *) malloc((n+1)*sizeof(double));
  A->updiag = (double *) malloc(n*sizeof(double));
  A->size = n+1;

  B = (struct tridiag *) malloc(sizeof(struct tridiag));
  B->subdiag = (double *) malloc(n*sizeof(double));
  B->diag = (double *) malloc((n+1)*sizeof(double));
  B->updiag = (double *) malloc(n*sizeof(double));
  B->size = n+1;

  C = (struct tridiag *) malloc(sizeof(struct tridiag));
  C->subdiag = (double *) malloc(m*sizeof(double));
  C->diag = (double *) malloc((m+1)*sizeof(double));
  C->updiag = (double *) malloc(m*sizeof(double));
  C->size = m+1;

  
  N_ij = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  N_ij->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  N_ij->lines = (struct cell **) malloc(2*sizeof(struct cell *));
  N_ij->nbCol = nbParamSigma;
  N_ij->nbLines = 2;
  for (i=0;i<2;i++)
    N_ij->lines[i] = NULL;
  for (j=0;j<nbParamSigma;j++)
    N_ij->columns[j] = NULL;
  
  above_Nij = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(2*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));


  invA = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    invA[i] = (double *) malloc((n+1)*sizeof(double));

  invC = (double **) malloc((m+1)*sizeof(double *));
  for (j=0;j<m+1;j++)
    invC[j] = (double *) malloc((m+1)*sizeof(double));
  
  tempMat = (double **) malloc((m+1)*sizeof(double *));
  for (j=0;j<m+1;j++)
    tempMat[j] = (double *) malloc((m+1)*sizeof(double));

  tempMat2 = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    tempMat2[i] = (double *) malloc((n+1)*sizeof(double));



  /* initialization of A and B*/
  A->diag[0] = 1;
  A->updiag[0] = 0;
  B->diag[0] = 0;
  B->updiag[0] = 0;
  for (i=0;i<n-1;i++){
    delta1 = y_coarseGrid[i+1]-y_coarseGrid[i];
    delta2 = y_coarseGrid[i+2]-y_coarseGrid[i+1];
    A->updiag[i+1] = delta1;
    A->diag[i+1] = delta1+delta2;
    A->subdiag[i] = delta2;
    B->updiag[i+1] = 3*delta1/delta2; 
    B->diag[i+1] = 3*(delta2/delta1-delta1/delta2);
    B->subdiag[i] = -3*delta2/delta1;
  }
  A->diag[n] = 1;
  A->subdiag[n-1] = 0;
  B->diag[n] = 0;
  B->subdiag[n-1] = 0;

  
 /* initialization of C and D*/
  C->diag[0] = 1;
  C->updiag[0] = 0;
  D->diag[0] = 0;
  D->updiag[0] = 0;
  for (j=0;j<m-1;j++){
    delta1 = T_coarseGrid[j+1]-T_coarseGrid[j];
    delta2 = T_coarseGrid[j+2]-T_coarseGrid[j+1];
    C->updiag[j+1] = delta1;
    C->diag[j+1] = delta1+delta2;
    C->subdiag[j] = delta2;
    D->updiag[j+1] = 3*delta1/delta2; 
    D->diag[j+1] = 3*(delta2/delta1-delta1/delta2);
    D->subdiag[j] = -3*delta2/delta1;
  }
  C->diag[m] = 1;
  C->subdiag[m-1] = 0;
  D->diag[m] = 0;
  D->subdiag[m-1] = 0;
  


  /*initialization of matGradF1 */
  for (line=0;line<nbParamSigma;line++)
    for (col=0;col<nbParamSigma;col++)
      matGradF1[line][col] = 0;




  /* inverses of A and C */
  tridiagMatrixInv(invA,A);

  tridiagMatrixInv(invC,C);


  /*computation of invC*D, stocked in a "full" matrix       */
  /*(we do not stock the blocks of zeros, but we know where */
  /*they are). The invC*D_i only differ by the place of the */
  /*non-empty columns                                       */

  tridiagMatMult(tempMat,D,invC);
  for (j=0;j<m+1;j++){
    for (i=0;i<m+1;i++)
      invC_times_D[j][i] = tempMat[j][i];
    invC_times_D[j][m+1] = invC[j][0];
    invC_times_D[j][m+2] = invC[j][m];
  }
 

  /*computation of invA*B_j, stocked in a "full" matrix  */


  tridiagMatMult(tempMat2,B,invA);
  for (j=0;j<n+1;j++){
    for (i=0;i<n+1;i++)
      invA_times_B[j][i] = tempMat2[j][i];
    invA_times_B[j][n+1] = invA[j][0];
    invA_times_B[j][n+2] = invA[j][n];
  }

  /*we stock once and for all the arrays of indices of non-empty columns               */
  /* of D_i for each i so as not to compute each of them several times.                */
  /*We stock these (n+1) arrays of dimension (m+3) in the (n+1)*(m+3) matrix D_indices */
  for (i=0;i<=n;i++){
    for (j=0;j<=m;j++)
      D_indices[i][j] = i + j*(n+1);
    D_indices[i][m+1] = (m+1)*(n+3) + i;
    D_indices[i][m+2] = (m+1)*(n+3) + (n+1) + i;
  }



  for (j=0;j<=m;j++){
    /* indices of non-empty columns of B_j */
    for (i=0;i<=n;i++)
      B_indices[j][i] = j*(n+1)+i;
    B_indices[j][n+1] = (n+1)*(m+1)+j;
    B_indices[j][n+2] = (n+2)*(m+1)+j;
    

    for (i=0;i<=n;i++){

      /* alloc and init of sentinelles */
      for (col=0;col<nbParamSigma;col++){
	sentinelles_columns[col] = (struct cell *) malloc(sizeof(struct cell));
	sentinelles_columns[col]->down = NULL;
	sentinelles_columns[col]->right = NULL;
      } 
      for (line=0;line<2;line++){
	sentinelles_lines[line] = (struct cell *) malloc(sizeof(struct cell));
	sentinelles_lines[line]->right = NULL;
	sentinelles_lines[line]->down = NULL;
      }  

      /* init of above_Nij */
      for (col=0;col<nbParamSigma;col++)
	above_Nij[col] = sentinelles_columns[col];
      

      // Computation of the first line of N_ij
      pt_line_Nij = sentinelles_lines[0];
      for (k=0;k<n+3;k++)
       	if (invA_times_B[i][k] != 0){
	  /* create a new cell at line 0, column col of Nij*/
	  col = B_indices[j][k];
	  pt_line_Nij->right = (struct cell *) malloc(sizeof(struct cell));
	  pt_line_Nij = pt_line_Nij->right;
	  pt_line_Nij->value = invA_times_B[i][k];
	  pt_line_Nij->line = 0;
	  pt_line_Nij->col = col;
	  pt_line_Nij->right = NULL;
	  pt_line_Nij->down = NULL;
	  /* update above_Nij */
	  above_Nij[col]->down = pt_line_Nij;
	  above_Nij[col] = above_Nij[col]->down;
	}

      /* end of the first line */
      N_ij->lines[0] = sentinelles_lines[0]->right;
      free(sentinelles_lines[0]);
      sentinelles_lines[0] = NULL;


      // Computation of the 2nd line of N_ij
      pt_line_Nij = sentinelles_lines[1];
      for (k=0;k<m+3;k++){
	//printf("j,k=%d,%d\n",j,k);
      	//printf("n,m=%d,%d\n",n,m);
      
       	if (invC_times_D[j][k] != 0){
	  /* create a new cell at line 0, column col of Nij*/
	  col = D_indices[i][k];
	  pt_line_Nij->right = (struct cell *) malloc(sizeof(struct cell));
	  pt_line_Nij = pt_line_Nij->right;
	  pt_line_Nij->value = invC_times_D[j][k];
	  pt_line_Nij->line = 1;
	  pt_line_Nij->col = col;
	  pt_line_Nij->right = NULL;
	  pt_line_Nij->down = NULL;
	  /* update above_Nij */
	  above_Nij[col]->down = pt_line_Nij;
	  above_Nij[col] = above_Nij[col]->down;
	}
      }

      /* end of second line */
       N_ij->lines[1] = sentinelles_lines[1]->right;
      free(sentinelles_lines[1]);
      sentinelles_lines[1] = NULL;

      /* end of the construction of Nij */
      for (col=0;col<nbParamSigma;col++){
	N_ij->columns[col] = sentinelles_columns[col]->down;
	free(sentinelles_columns[col]);
	sentinelles_columns[col] = NULL;
	above_Nij[col] = NULL;
      }

      //computation of M_ij = transpose(N_ij)*N_ij
      sparseMatMult(&M_ij,N_ij);

      //matGradF1 = matGradF1 + transpose(N_ij)*N_ij
      regPlusSparse(matGradF1,M_ij,matGradF1);
      
      desallocSparseMat_cells(N_ij);
      desallocSparseMat_cells(M_ij); 
      free(M_ij->lines);
      M_ij->lines = NULL;
      free(M_ij->columns);
      M_ij->columns = NULL;
      free(M_ij);
      M_ij = NULL;
    }
  }       

  for (i=0;i<nbParamSigma;i++)
    for(j=0;j<nbParamSigma;j++)
      matGradF1[i][j] = 2*matGradF1[i][j];


 /*memory desallocation */

  free(sentinelles_lines);
  free(sentinelles_columns);

  free(above_Nij);

  free(N_ij->lines);
  free(N_ij->columns);
  free(N_ij);

  free(A->subdiag);
  free(A->diag);
  free(A->updiag);
  free(A);

  free(B->subdiag);
  free(B->diag);
  free(B->updiag);
  free(B);

  free(C->subdiag);
  free(C->diag);
  free(C->updiag);
  free(C);

  for (i=0;i<n+1;i++)
    free(invA[i]);
  free(invA);

  for (j=0;j<m+1;j++)
    free(invC[j]);
  free(invC);

  for (j=0;j<m+1;j++)
    free(tempMat[j]);
  free(tempMat);

  for (i=0;i<n+1;i++)
    free(tempMat2[i]);
  free(tempMat2);

}





/***************************************************************************************************************************/
/*                             COMPUTATION OF GRAD(G)                                                                     */
/***************************************************************************************************************************/


void computeGrad_G(double *grad_G, double *sigma_param, int N, int M, int n, int m, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, double r, double q, int optionType, double S_0, double **invA_times_B, double **invC_times_D, int **D_indices, int **B_indices, struct tridiag *D, struct marketData **marketOptionPrices, int nbMaturites){
/* computes the gradient of G 
   OUTPUT:
   - grad_G (gradient evaluated in sigma_param)
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
   - optionType : type of the option (1 for call, 0 for put) 
   - optionPrices : market prices of options
   - nbMaturites : nb of maturites for which we know market prices 
   - invA_times_B (computed in computeGrad_F1)
   - invC_times_D (computed in computeGrad_F1)
   - D_indices(computed in computeGrad_F1)
   - B_indices(computed in computeGrad_F1)
   - D (computed in computeGrad_F1)
*/

  double **U,**prices,**sigmaFineGrid,**sigmaCoarseGrid;
  struct sparseMat *transp_deriv_phi;
  struct derivData *interpolData;
  struct marketData *pt_data;
  double theta; /* param of the PDE resolution */
  int i,j,k,l,line,col;
  int n_d; /* nb of market data */
  double *phi_U,*U_tilde;
  double a,b;
  struct tridiag *A,*transp_deriv_b;
  struct tridiagSystem *tridiagSyst;
  struct bidiagSystem *bidiagSyst;
  double *vect,*vect2;
  struct sparseMat *R;
  struct sparseMat *deriv_Au_j,*prev_deriv_Au_j,*deriv_b_j,*deriv_Hu_j;
  int nbParamSigma;
  double timeStep;
  struct sparseMat *mat_j,*transp_mat_j;
  double *u_bar,*u_bar_prev,*elem_j;
  double **U_bar;


  nbParamSigma = (n+1)*(m+1)+2*(n+1)+2*(m+1)+4;
  theta = .5;
  /* compute n_d */
  n_d = 0;
  for (i=0;i<nbMaturites;i++){
    pt_data = marketOptionPrices[i];
    while (pt_data != NULL){
      n_d++;
      pt_data = pt_data->next;
    }
  }


  u_bar = (double *)malloc((N-1)*sizeof(double));

  elem_j = (double *)malloc(nbParamSigma*sizeof(double));

  phi_U = (double *)malloc(n_d*sizeof(double));
  U_tilde = (double *)malloc(n_d*sizeof(double));

  sigmaCoarseGrid = (double **) malloc((n+1)*sizeof(double *));
  for (i=0;i<n+1;i++)
    sigmaCoarseGrid[i] = (double *) malloc((m+1)*sizeof(double));
 
  sigmaFineGrid = (double **) malloc((N+1)*sizeof(double *));
  for (i=0;i<N+1;i++)
    sigmaFineGrid[i] = (double *) malloc((M+1)*sizeof(double));
  
  interpolData = (struct derivData *) malloc(sizeof(struct derivData));
  interpolData->deriv_y_0 = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_y_n = (double *) malloc((m+1)*sizeof(double));
  interpolData->deriv_T_0 = (double *) malloc((n+1)*sizeof(double));
  interpolData->deriv_T_m = (double *) malloc((n+1)*sizeof(double));

  prices = (double **)malloc((N+1)*sizeof(double *));
  for (i=0;i<N+1;i++)
    prices[i] = (double *)malloc((M+1)*sizeof(double));

  U = (double **)malloc((N-1)*sizeof(double *));
  for (i=0;i<N-1;i++)
    U[i] = (double *)malloc((M+1)*sizeof(double));

  U_bar = (double **)malloc((N-1)*sizeof(double *));
  for (i=0;i<N-1;i++)
    U_bar[i] = (double *)malloc((M+1)*sizeof(double));

  /* memory allocation and init for the matrix (*trans_deriv_phi) */
  transp_deriv_phi = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  transp_deriv_phi->columns = (struct cell **) malloc(n_d*sizeof(struct cell *));
  transp_deriv_phi->lines = (struct cell **) malloc((N-1)*sizeof(struct cell *));
  transp_deriv_phi->nbCol = n_d;
  transp_deriv_phi->nbLines = N-1;
 
  /* memory allocation for the operator A (of type struct tridiag *) */
  A = (struct tridiag *) malloc(sizeof(struct tridiag));
  A->subdiag = (double *) malloc((N-2)*sizeof(double));
  A->diag = (double *) malloc((N-1)*sizeof(double));
  A->updiag = (double *) malloc((N-2)*sizeof(double));
  A->size = N-1;
  
  /* memory allocation for the operator transp_deriv_b (of type struct tridiag *) */
  transp_deriv_b = (struct tridiag *) malloc(sizeof(struct tridiag));
  transp_deriv_b->subdiag = (double *) malloc((N-2)*sizeof(double));
  transp_deriv_b->diag = (double *) malloc((N-1)*sizeof(double));
  transp_deriv_b->updiag = (double *) malloc((N-2)*sizeof(double));
  transp_deriv_b->size = N-1;
  
  /* memory alloc for the tridiag system */
  tridiagSyst = (struct tridiagSystem *)malloc(sizeof(struct tridiagSystem));
  tridiagSyst->T = (struct tridiag *) malloc(sizeof(struct tridiag));
  tridiagSyst->T->subdiag = (double *) malloc((N-2)*sizeof(double));
  tridiagSyst->T->diag = (double *) malloc((N-1)*sizeof(double));
  tridiagSyst->T->updiag = (double *) malloc((N-2)*sizeof(double));
  tridiagSyst->T->size = N-1;
  tridiagSyst->b = (double *) malloc((N-1)*sizeof(double));
  tridiagSyst->size = N-1;

  /* memory alloc for the bidiag system */
  bidiagSyst = (struct bidiagSystem *)malloc(sizeof(struct bidiagSystem));
  bidiagSyst->T = (struct bidiag *) malloc(sizeof(struct bidiag));
  bidiagSyst->T->subdiag = (double *) malloc((N-2)*sizeof(double));
  bidiagSyst->T->diag = (double *) malloc((N-1)*sizeof(double));
  bidiagSyst->T->size = N-1;
  bidiagSyst->b = (double *) malloc((N-1)*sizeof(double));
  bidiagSyst->size = N-1;

  u_bar_prev = (double *)malloc((N-1)*sizeof(double));
  for (i=0;i<N-1;i++)
    u_bar_prev[i] = 0;
 
  vect = (double *)malloc(n_d*sizeof(double));
  vect2 = (double *)malloc((N-1)*sizeof(double));

  deriv_Au_j = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  deriv_Au_j->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  deriv_Au_j->lines = (struct cell **) malloc((N-1)*sizeof(struct cell *));
  deriv_Au_j->nbCol = nbParamSigma;
  deriv_Au_j->nbLines = N-1; 

  prev_deriv_Au_j = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  prev_deriv_Au_j->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  prev_deriv_Au_j->lines = (struct cell **) malloc((N-1)*sizeof(struct cell *));
  prev_deriv_Au_j->nbCol = nbParamSigma;
  prev_deriv_Au_j->nbLines = N-1; 

  deriv_b_j = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  deriv_b_j->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  deriv_b_j->lines = (struct cell **) malloc((N-1)*sizeof(struct cell *));
  deriv_b_j->nbCol = nbParamSigma;
  deriv_b_j->nbLines = N-1; 

  deriv_Hu_j = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  deriv_Hu_j->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  deriv_Hu_j->lines = (struct cell **) malloc((N-1)*sizeof(struct cell *));
  deriv_Hu_j->nbCol = nbParamSigma;
  deriv_Hu_j->nbLines = N-1; 


  /* init grad_G */
  for (i=0;i<nbParamSigma;i++)
    grad_G[i] = 0;
 

 /*interpolation of sigma (coarse grid, defined by sigma_param) to a finer grid (much more points => much more precise) */
     /* 1) sigma_param => sigmaCoarseGrid,interpolData */
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


  /* 2) interpolation of sigmaCoarseGrid,interpolData => computation of sigmaFineGrid */
  interpole(sigmaFineGrid,sigmaCoarseGrid,n,m,N,M,interpolData,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid);

  /* for (i=0;i<N;i++){ */
/*     for (j=0;j<M;j++) */
/* 	printf("%lf ",sigmaFineGrid[i][j]); */
/*     printf("\n"); */
/*   } */

  /* STEP 1 : solve the state equation : from sigma_param (vector of degrees of freedom), we compute U of size (N-1,M+1) */
  solve(optionType,prices,S_0,N,M,r,q,theta,f,sigmaFineGrid,y_fineGrid,T_fineGrid);
  for (i=1;i<N;i++)
    for (j=0;j<M+1;j++)
      U[i-1][j] = prices[i][j];

  /* U[i-1][j] = u_i^j , i=1,...,N and j=0,...,M */


  /* STEP 2 : compute the adjoint state U_bar of size (N-1,M+1) ************************************8*/
  

  /* 1) compute phi_U and U_tilde of size (n_d,1) */
  j = 0;
  l = 0;  
  for (i=0;i<nbMaturites;i++){
    pt_data = marketOptionPrices[i];
    k = 0;
    while (pt_data != NULL){
      k = find_index(log(pt_data->strike),y_fineGrid,N,k); /* y_fineGrid[k] <= log(pt_data->strike) < y_fineGrid[k+1] */
      l = find_index(pt_data->maturity,T_fineGrid,M,l); /* T_fineGrid[l] <= pt_data->maturity < T_fineGrid[l+1] */
      a = (T_fineGrid[l+1]-pt_data->maturity)/(T_fineGrid[l+1]-T_fineGrid[l]);
      b = (y_fineGrid[k+1]-log(pt_data->strike))/(y_fineGrid[k+1]-y_fineGrid[k]);
      phi_U[j] = (1-b)*((1-a)*prices[k+1][l+1] + a*prices[k][l+1]) + b*((1-a)*prices[k+1][l] + a*prices[k][l]);
      U_tilde[j] = pt_data->price;
      j++;
      pt_data = pt_data->next;
    }
  }

  //printf("phi_U   U_tilde\n");
  //for (i=0;i<n_d;i++)
  //printf("%lf  %lf\n",phi_U[i],U_tilde[i]);

  /* 2) j=M : compute u_bar^M */
  /* 2.1) compute A = A^M */ 
  buildOperator(A,r,q,S_0,M,N,sigmaFineGrid,y_fineGrid);
  /* 2.2) compute transp_deriv_phi = [\frac{\partial \varphi}{\partial \u^M}(U)]^T of size (N-1,n_d) */
  buildTranspDerivPhi(&transp_deriv_phi,M,n_d,N,M,nbMaturites,y_fineGrid,T_fineGrid,marketOptionPrices);
  /* 2.3) build tridiag system who will have u_bar^M for solution */
  buildAdjointTridiagSyst(tridiagSyst,M,A,u_bar_prev,theta,T_fineGrid[M]-T_fineGrid[M-1],0,transp_deriv_phi,phi_U,U_tilde,n_d,M);
  /* 2.4) change the tridiag syst to a bidiag system */
  tridiagToBidiagSyst(bidiagSyst,tridiagSyst);
  /* 2.5) solve the bidiag syst (compute u_bar = u_bar^M) */
  solveSyst(u_bar,bidiagSyst);
  /* 2.6) stock u_bar in U_bar and init u_bar_prev*/
  for (line=0;line<N-1;line++){
    U_bar[line][M] = u_bar[line];  
    u_bar_prev[line] = u_bar[line];
  }
  /* 2.7) free the cells of the matrix transp_deriv_phi */
  desallocSparseMat_cells(transp_deriv_phi);


  /* 3) for j=M-1...0, compute u_bar^j */
  for (j=M-1;j>=1;j--){
    /* 3.1) compute A^j */
    buildOperator(A,r,q,S_0,j,N,sigmaFineGrid,y_fineGrid);
    /* 3.2) compute transp_deriv_phi of size (N-1,n_d) */
    buildTranspDerivPhi(&transp_deriv_phi,j,n_d,N,M,nbMaturites,y_fineGrid,T_fineGrid,marketOptionPrices);
    /* 3.3) build tridiag system who will have u_bar^j for solution */
    buildAdjointTridiagSyst(tridiagSyst,j,A,u_bar_prev,theta,T_fineGrid[j]-T_fineGrid[j-1],T_fineGrid[j+1]-T_fineGrid[j],transp_deriv_phi,phi_U,U_tilde,n_d,M);
    /* 3.4) change the tridiag syst to a bidiag system */
    tridiagToBidiagSyst(bidiagSyst,tridiagSyst);
    /* 3.5) solve the bidiag syst (compute u_bar = u_bar^j) */
    solveSyst(u_bar,bidiagSyst);
    /* 3.6) stock u_bar in U_bar and update u_bar_prev*/
    for (line=0;line<N-1;line++){
      U_bar[line][j] = u_bar[line];  
      u_bar_prev[line] = u_bar[line];
    }
    /* 3.7) free the cells of the matrix transp_deriv_phi */
    desallocSparseMat_cells(transp_deriv_phi);
  }

  /* 4) compute u_bar^0 */
  /* 4.1) compute A^0 */
  buildOperator(A,r,q,S_0,0,N,sigmaFineGrid,y_fineGrid);
  /* 4.2) compute transp_deriv_phi of size (N-1,n_d) */
  buildTranspDerivPhi(&transp_deriv_phi,0,n_d,N,M,nbMaturites,y_fineGrid,T_fineGrid,marketOptionPrices);
  /* 4.3) compute u_bar = u_bar^0 */
  for (i=0;i<n_d;i++)
    vect[i] = phi_U[i]-U_tilde[i];
  sparseTimesVect(u_bar,transp_deriv_phi,vect);
  /* we have to add \frac{\partial b^j}{\partial u^j}^T \bar{u}^{j+1} */
  /* build \frac{\partial b^j}{\partial u^j}^T */
  for (i=0;i<N-2;i++){
    transp_deriv_b->subdiag[i] = -theta*(T_fineGrid[1]-T_fineGrid[0])*A->updiag[i];
    transp_deriv_b->diag[i] = 1-theta*(T_fineGrid[1]-T_fineGrid[0])*A->diag[i];
    transp_deriv_b->updiag[i] = -theta*(T_fineGrid[1]-T_fineGrid[0])*A->subdiag[i];
  }   
  transp_deriv_b->diag[N-2] = 1-theta*(T_fineGrid[1]-T_fineGrid[0])*A->diag[N-2];
  /* multiplie by u_bar_prev */
  tridiagTimesVect(vect2,transp_deriv_b,u_bar_prev);
  /* add to tridiagSyst->b */
  for (i=0;i<N-1;i++)
    u_bar[i] = u_bar[i] + vect2[i];
  /* 4.4) stock u_bar in U_bar and update u_bar_prev*/
  for (line=0;line<N-1;line++){
    U_bar[line][0] = u_bar[line];  
    u_bar_prev[line] = u_bar[line];
  }
  /* 4.5) free the cells of the matrix transp_deriv_phi */
  desallocSparseMat_cells(transp_deriv_phi);
  

  /* for (i=0;i<N-1;i++){ */
/*     for (j=0;j<n_d;j++) */
/*       printf("%lf ",U_bar[i][j]); */
/*     printf("\n"); */
/*   } */

  /**************************************************************************************************************************/
  /* STEP 3 : computation of grad G
  /* computation of the sparse matrix R of size 4(n+1)(M+1),nbParamSigma */
  compute_R(&R,n,m,invA_times_B,invC_times_D,B_indices,D_indices,D,T_coarseGrid);
  //printSparseMat(R);



  /* init prev_derv_Au_j */
  compute_deriv_Au_j(&prev_deriv_Au_j,M,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid,N,M,n,m,R,sigmaFineGrid,sigmaCoarseGrid,U);
 
  

  for (j=M-1;j>=0;j--){

    printf("*********************************   j=%d   ****************************************\n",j); 
    /* for each time index j */
    timeStep = T_fineGrid[j+1]-T_fineGrid[j];
    /* 1) compute deriv_Au_j = \frac{\partial A^ju^j}{\partial \Sigma}*/
    compute_deriv_Au_j(&deriv_Au_j,j,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid,N,M,n,m,R,sigmaFineGrid,sigmaCoarseGrid,U);
   
    //printf("deriv_Au_j\n");
    //printSparseMat(deriv_Au_j);
    
    
    /* deriv_Au_j is now computed */
    /* computation of deriv_b_j = \frac{\partial b^j}{\partial \Sigma} */
    affectSparse(&deriv_b_j,deriv_Au_j);
    scalarTimesSparse(-theta*timeStep,&deriv_b_j);
    
    //printSparseMat(deriv_b_j);

    /* computation of deriv_Hu_j = \frac{\partial H^ju^{j+1}}{\partial \Sigma} */
    affectSparse(&deriv_Hu_j,prev_deriv_Au_j);
    scalarTimesSparse((1-theta)*timeStep,&deriv_Hu_j);
    
    //printSparseMat(deriv_Hu_j);

    /* update deriv_G */
    sparsePlusSparse(&mat_j,deriv_b_j,deriv_Hu_j);
    for (i=0;i<N-1;i++)
    u_bar[i] = U_bar[i][j];

    transpSparseMat(&transp_mat_j,mat_j);
    
    sparseTimesVect(elem_j,transp_mat_j,u_bar);
    for (i=0;i<nbParamSigma;i++)
      grad_G[i] = grad_G[i] + elem_j[i];
    
    
    desallocSparseMat_cells(prev_deriv_Au_j);
    affectSparse(&prev_deriv_Au_j,deriv_Au_j);
    desallocSparseMat_cells(deriv_Au_j);
    desallocSparseMat_cells(deriv_b_j);
    desallocSparseMat_cells(deriv_Hu_j); 
    desallocSparseMat_cells(mat_j);
    free(mat_j->lines);
    free(mat_j->columns);
    free(mat_j); 
    desallocSparseMat_cells(transp_mat_j);
    free(transp_mat_j->lines);
    free(transp_mat_j->columns);
    free(transp_mat_j); 
  }

  /* free memory */

  free(deriv_Au_j->lines);
  free(deriv_Au_j->columns);
  free(deriv_Au_j); 

  desallocSparseMat_cells(prev_deriv_Au_j);
  free(prev_deriv_Au_j->lines);
  free(prev_deriv_Au_j->columns);
  free(prev_deriv_Au_j); 
  
  free(deriv_b_j->lines);
  free(deriv_b_j->columns);
  free(deriv_b_j);

  free(deriv_Hu_j->lines);
  free(deriv_Hu_j->columns);
  free(deriv_Hu_j); 

  free(elem_j); 
  
  desallocSparseMat_cells(R);
  free(R->lines);
  free(R->columns);
  free(R);
 
  free(phi_U);
  free(U_tilde);

  for (i=0;i<n+1;i++)
    free(sigmaCoarseGrid[i]);
  free(sigmaCoarseGrid);

  for (i=0;i<N+1;i++)
    free(sigmaFineGrid[i]);
  free(sigmaFineGrid);

  free(interpolData->deriv_y_0);
  free(interpolData->deriv_y_n);
  free(interpolData->deriv_T_0); 
  free(interpolData->deriv_T_m); 
  free(interpolData);

  for (i=0;i<N+1;i++)
    free(prices[i]);
  free(prices);

  for (i=0;i<N-1;i++)
    free(U[i]);
  free(U);

  for (i=0;i<N-1;i++)
    free(U_bar[i]);
  free(U_bar);
  
  free(transp_deriv_phi->columns);
  free(transp_deriv_phi->lines);
  free(transp_deriv_phi);

  free(A->subdiag);
  free(A->diag);
  free(A->updiag);
  free(A);

  free(transp_deriv_b->subdiag);
  free(transp_deriv_b->diag);
  free(transp_deriv_b->updiag);
  free(transp_deriv_b);
  
  free(tridiagSyst->T->subdiag);
  free(tridiagSyst->T->diag);
  free(tridiagSyst->T->updiag);
  free(tridiagSyst->T);
  free(tridiagSyst->b);
  free(tridiagSyst);
  
  free(bidiagSyst->T->subdiag);
  free(bidiagSyst->T->diag);
  free(bidiagSyst->T);
  free(bidiagSyst->b);
  free(bidiagSyst);

  free(u_bar);
  free(u_bar_prev);

  free(vect);
  free(vect2); 

 

}
