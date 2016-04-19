#include "utilsGrad.h"
#include "DupirePDE.h"
#include "spline.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

void buildTranspDerivPhi(struct sparseMat **pt_transp_deriv_phi_j, int j, int n_d, int N, int M, int nbMaturites, double *y_fineGrid, double *T_fineGrid, struct marketData **marketOptionPrices){

  int line,k,l,i,index_maturity;
  struct cell **pt_left;
  struct cell **sentinelles_columns,**sentinelles_lines;
  struct cell *pt_col;
  struct marketData *pt_data;
  double a,b;
  
  /* memory allocation for the sentinelles */
  sentinelles_columns = (struct cell **) malloc(n_d*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc((N-1)*sizeof(struct cell *));

  pt_left = (struct cell **) malloc((N-1)*sizeof(struct cell *));
 
  for (i=0;i<n_d;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;
  } 
  for (line=0;line<N-1;line++){
    sentinelles_lines[line] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[line]->right = NULL;
    sentinelles_lines[line]->down = NULL;
  }  

  for (line=0;line<N-1;line++)
    pt_left[line] = sentinelles_lines[line];


  i = 0;
  l = 0;
  for (index_maturity=0;index_maturity<nbMaturites;index_maturity++){
    pt_data = marketOptionPrices[index_maturity];
    while (pt_data != NULL){
      /* for each column of transp_deriv_phi */
      printf("col = %d\n",i);
      pt_col = sentinelles_columns[i];
      l = find_index(pt_data->maturity,T_fineGrid,M,0);
      k = find_index(log(pt_data->strike),y_fineGrid,N,0);
      a = (T_fineGrid[l+1]-pt_data->maturity)/(T_fineGrid[l+1]-T_fineGrid[l]);
      b = (y_fineGrid[k+1]-log(pt_data->strike))/(y_fineGrid[k+1]-y_fineGrid[k]);
    
      if (j==l){
	/* only the elements on line k-1 and k are not zero (unless k=0 or k=N-1) */

	if (k!=0){
	  pt_col->down = (struct cell *) malloc(sizeof(struct cell));
	  pt_col = pt_col->down;
	  pt_col->line = k-1;
	  pt_col->col = i;
	  pt_col->value = a*b;
	  pt_col->down = NULL;
	  pt_col->right = NULL;
	  pt_left[k-1]->right = pt_col;
	  pt_left[k-1] = pt_left[k-1]->right;
	}
	if (k!=N-1){
	pt_col->down = (struct cell *) malloc(sizeof(struct cell));
	pt_col = pt_col->down;
	pt_col->line = k;
	pt_col->col = i;
	pt_col->value = (1-a)*b;
	pt_col->down = NULL;
	pt_col->right = NULL;
	pt_left[k]->right = pt_col;
	pt_left[k] = pt_left[k]->right;
	}

      }  
      else if (j==l+1){
	/* only the elements on line k-1 and k are not zero (unless k=0 or k=N-1) */

	if (k!=0){
	  pt_col->down = (struct cell *) malloc(sizeof(struct cell));
	  pt_col = pt_col->down;
	  pt_col->line = k-1;
	  pt_col->col = i;
	  pt_col->value = a*(1-b);
	  pt_col->down = NULL;
	  pt_col->right = NULL;
	  pt_left[k-1]->right = pt_col;
	  pt_left[k-1] = pt_left[k-1]->right;
	}
	if (k!=N-1){
	pt_col->down = (struct cell *) malloc(sizeof(struct cell));
	pt_col = pt_col->down;
	pt_col->line = k;
	pt_col->col = i;
	pt_col->value = (1-a)*(1-b);
	pt_col->down = NULL;
	pt_col->right = NULL;
	pt_left[k]->right = pt_col;
	pt_left[k] = pt_left[k]->right;
	}

      }  

      /* end of columns "col" */
      (*pt_transp_deriv_phi_j)->columns[i] = sentinelles_columns[i]->down;
      free(sentinelles_columns[i]);
      i++; 
      pt_data = pt_data->next;
    }

  }

  for (line=0;line<N-1;line++){
    (*pt_transp_deriv_phi_j)->lines[line] = sentinelles_lines[line]->right;
    free(sentinelles_lines[line]);
    pt_left[line] = NULL;
  }

  free(pt_left);
  free(sentinelles_lines);
  free(sentinelles_columns);

}

void buildAdjointTridiagSyst(struct tridiagSystem *tridiagSyst, int j, struct tridiag *A, double *u_prev, double theta, double timeStep, double timeStep2, struct sparseMat *transp_deriv_phi, double *phi_U, double *U_tilde, int n_d, int M){

  int i,N;
  double *vect,*vect2;
  struct tridiag *transp_deriv_b;

  N = A->size+1;

  /* memory allocation for the operator transp_deriv_b (of type struct tridiag *) */
  transp_deriv_b = (struct tridiag *) malloc(sizeof(struct tridiag));
  transp_deriv_b->subdiag = (double *) malloc((N-2)*sizeof(double));
  transp_deriv_b->diag = (double *) malloc((N-1)*sizeof(double));
  transp_deriv_b->updiag = (double *) malloc((N-2)*sizeof(double));
  transp_deriv_b->size = N-1;
  
  vect = (double *)malloc(n_d*sizeof(double));
  vect2 = (double *)malloc((N-1)*sizeof(double));

  for (i=0;i<n_d;i++)
    vect[i] = phi_U[i]-U_tilde[i];

  /* build left hand side tridiag matrix */
  for (i=0;i<N-2;i++){
    tridiagSyst->T->subdiag[i] = (1-theta)*timeStep*A->subdiag[i];
    tridiagSyst->T->diag[i] = 1+(1-theta)*timeStep*A->diag[i];
    tridiagSyst->T->updiag[i] = (1-theta)*timeStep*A->updiag[i];
  }   
  tridiagSyst->T->diag[N-2] = 1+(1-theta)*timeStep*A->diag[N-2];


  /* build the right hand side vector (b) */
  sparseTimesVect(tridiagSyst->b,transp_deriv_phi,vect);

  if (j!=M){
    /* we have to add \frac{\partial b^j}{\partial u^j}^T \bar{u}^{j+1} */
    /* build \frac{\partial b^j}{\partial u^j}^T */
    for (i=0;i<N-2;i++){
      transp_deriv_b->subdiag[i] = -theta*timeStep2*A->updiag[i];
      transp_deriv_b->diag[i] = 1-theta*timeStep2*A->diag[i];
      transp_deriv_b->updiag[i] = -theta*timeStep2*A->subdiag[i];
    }   
    transp_deriv_b->diag[N-2] = 1-theta*timeStep2*A->diag[N-2];
    
    /* multiplie by u_prev */
    tridiagTimesVect(vect2,transp_deriv_b,u_prev);
    
    /* add to tridiagSyst->b */
    for (i=0;i<N-1;i++)
      tridiagSyst->b[i] = tridiagSyst->b[i] + vect2[i];

  }

  free(vect);
  free(vect2);
  free(transp_deriv_b->subdiag);
  free(transp_deriv_b->diag);
  free(transp_deriv_b->updiag);
  free(transp_deriv_b);
}



void findClosestIndices(int *ok, int *k, int *l, double K, double T, double *y_grid, double *t_grid, int N, int M){

  int i;

  i = 0;
  while (i<N+1 && y_grid[i]<log(K))
    i++;
  if (i == 0)
    *ok = 0;
  else{ /* log(K) is between y_grid[i-1] and y_grid[i] */
    *ok = 1;
    if (fabs(y_grid[i-1]-log(K)) <= fabs(y_grid[i]-log(K)))
      /* log(K) closer to y_grid[i=1] */
      *k = i-1;
    else
      *k = i;
  }

  i = 0;
  while (i<M+1 && t_grid[i]<T)
    i++;
  if (i == 0)
    *ok = 0;
  else{ /* T is between t_grid[i-1] and t_grid[i] */
    *ok = 1;
    if (fabs(t_grid[i-1]-T) <= fabs(t_grid[i]-T))
      /* log(K) closer to y_grid[i=1] */
      *l = i-1;
    else
      *l = i;
  }

}



void matMultVect(double *A, double **B, double *C, int nbLinesB, int nbColB){

  int i,j;
  double sum;

  for (i=0;i<nbLinesB;i++){
    sum = 0;
    for (j=0;j<nbColB;j++)
      sum = sum + B[i][j]*C[j];
    A[i] = sum;
  }

}

void vectPlus(double *A, double *B, int nbLines){

  int i;

  for (i=0;i<nbLines;i++)
    A[i] = A[i] + B[i];

}




double g(double y, double **U, struct marketData **data, int nbMaturites, int N, int M, double T_max, double *y_grid,  double *t_grid){

  /* T_max = max_i T_i (market data) */
  struct marketData *pt_data;
  int i,j;
  double res;

  res = 0;
  /* 1) are the any data corresponding to maturity T_max? */
  /* yes because T_max = max_i T_i = data[nbMaturites-1]->maturity*/
  
  /* 2) find (if exists) the data corresponding to strike exp(y) in the last line of the market data */
  pt_data = data[nbMaturites-1];
  while(pt_data != NULL && log(pt_data->strike) < y)
    pt_data = pt_data->next;
  if (pt_data != NULL && log(pt_data->strike) == y){
    /* we found the option of maturity T_max and strike exp(y) */
    /* now we look for the indices (i,j) so that (y,T_max) is in Rij of the fine grid */
    j = find_index(T_max,t_grid,M,0);
    i = find_index(y,y_grid,N,0);
    res = U[i][j]-pt_data->price;
  }
 
  return res;
  
}




void compute_R(struct sparseMat **R, int n, int m, double **invA_times_B, double **invC_times_D, int **B_indices, int **D_indices, struct tridiag *D, double *T_coarseGrid){

  int i,j,k,l,nbParamSigma,index_line_R,line,nbElem;
  int *H0_indices, *Hm_indices;
  struct sparseMat *G,*Li,*Pi,*G_times_Li,*Si;
  struct sparseMat **invC_times_S_tab;
  struct cell *newCell, *ptCell, *ptPrev;
  struct cell **above;
  struct cell **sentinelles_columns, **sentinelles_lines;
  double **invC;
  double delta1,delta2;
  struct tridiag *C;
  
  /* memory allocation and initialization */
  nbParamSigma = (n+1)*(m+1)+2*(n+1)+2*(m+1)+4;

  if (m>1){
    invC = (double **) malloc((m-1)*sizeof(double *));
    for (k=0;k<m-1;k++)
      invC[k] = (double *) malloc((m-1)*sizeof(double));
  }

  H0_indices = (int *)malloc((n+3)*sizeof(int));
  Hm_indices = (int *)malloc((n+3)*sizeof(int));

  if (m>1){

    invC_times_S_tab = (struct sparseMat **)malloc((n+1)*sizeof(struct sparseMat *));

    Li = (struct sparseMat *) malloc(sizeof(struct sparseMat));
    Li->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
    Li->lines = (struct cell **) malloc((m+1)*sizeof(struct cell *));
    Li->nbCol = nbParamSigma;
    Li->nbLines = m+1; 
    
    Pi = (struct sparseMat *) malloc(sizeof(struct sparseMat));
    Pi->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
    Pi->lines = (struct cell **) malloc((m-1)*sizeof(struct cell *));
    Pi->nbCol = nbParamSigma;
    Pi->nbLines = m-1;
  }

  (*R) = (struct sparseMat *)malloc(sizeof(struct sparseMat)); 
  (*R)->nbCol = nbParamSigma;
  (*R)->nbLines = 4*(n+1)*(m+1); 
  (*R)->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  sentinelles_columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  (*R)->lines = (struct cell **) malloc(4*(n+1)*(m+1)*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc(4*(n+1)*(m+1)*sizeof(struct cell *));

  for (i=0;i<nbParamSigma;i++){
    sentinelles_columns[i] = (struct cell *) malloc(sizeof(struct cell)); /* sentinelles */
    sentinelles_columns[i]->down = NULL;
    sentinelles_columns[i]->right = NULL;
  }

  for (i=0;i<4*(m+1)*(n+1);i++){
    sentinelles_lines[i] = (struct cell *) malloc(sizeof(struct cell)); /* sentinelles */
    sentinelles_lines[i]->right = NULL;
    sentinelles_lines[i]->down = NULL;
  }
  

  above = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));

  if (m>1){
    C = (struct tridiag *) malloc(sizeof(struct tridiag));
    C->subdiag = (double *) malloc((m-2)*sizeof(double));
    C->diag = (double *) malloc((m-1)*sizeof(double));
    C->updiag = (double *) malloc((m-2)*sizeof(double));
    C->size = m-1;
  }

  /* end of memory allocation */


  /* initialization of C */

  if (m>1){
    delta1 = T_coarseGrid[1]-T_coarseGrid[0];
    delta2 = T_coarseGrid[2]-T_coarseGrid[1];
    C->diag[0] = delta1+delta2;
    C->updiag[0] = delta1;
    for (k=1;k<m-2;k++){
      delta1 = T_coarseGrid[k+1]-T_coarseGrid[k];
      delta2 = T_coarseGrid[k+2]-T_coarseGrid[k+1];
      C->updiag[k] = delta1;
      C->diag[k] = delta1+delta2;
      C->subdiag[k-1] = delta2;
    }
    delta1 = T_coarseGrid[m-1]-T_coarseGrid[m-2];
    delta2 = T_coarseGrid[m]-T_coarseGrid[m-1];
    C->diag[m-2] = delta1+delta2;
    C->subdiag[m-3] = delta2;
    

    /* computation of invC */
    tridiagMatrixInv(invC,C);
    
  }

  
  for (k=0;k<nbParamSigma;k++)
    above[k] = sentinelles_columns[k];


  //printf("avant R1\n");
  /* R1 */
  for (i=0;i<(n+1)*(m+1);i++){
    (*R)->lines[i] = (struct cell *)malloc(sizeof(struct cell));
    (*R)->lines[i]->value = 1;
    (*R)->lines[i]->col = i;
    (*R)->lines[i]->line = i; 
    (*R)->lines[i]->right = NULL;
    (*R)->columns[i] = (*R)->lines[i];
    above[i]->down = (*R)->lines[i];
    above[i] = above[i]->down;
    free(sentinelles_lines[i]);
    sentinelles_lines[i] = NULL;
  }
  //printf("avant R2\n");
  /* R2 */
  index_line_R = (n+1)*(m+1);
  for (j=0;j<=m;j++){
    /* construction of the matrix R_2^j */
    for (i=0;i<=n;i++){
      /* ith line of R_2^j = jth line of invC_times_Di */
      newCell = sentinelles_lines[index_line_R];
      l=0; /* index of D_indices */
      k=0;/* index of columns of R */
      while (k<nbParamSigma && l<m+3){
	if (D_indices[i][l] == k){
	  if (invC_times_D[j][l] != 0){
	    //k_th column of invC_times_Di not empty
	    newCell->right = (struct cell *)malloc(sizeof(struct cell));
	    newCell = newCell->right;
	    newCell->value = invC_times_D[j][l];
	    newCell->col = k;
	    //printf("newCell->col = %d\n",newCell->col);
	    newCell->line = index_line_R;
	    newCell->down = NULL;
	    newCell->right = NULL;
	    above[k]->down = newCell;
	    above[k] = above[k]->down;
	  }
	  k++;
	  l++;
	  
	}
	else /* if (k < D_indices[i][l]) */
	  //k_th column empty
	  k++;
      }
      /* end of line (of R) */
      
      (*R)->lines[index_line_R] = sentinelles_lines[index_line_R]->right;
      free(sentinelles_lines[index_line_R]);
      sentinelles_lines[index_line_R] = NULL;
      index_line_R++; 
      
    }
    
  }
  //printf("avant R3\n");
  /* R3 */
  
  for (j=0;j<=m;j++){
    /* construction of the matrix R_3^j */
    for (i=0;i<=n;i++){
      /* ith line of R_3^j = ith line of invA_times_Bj */
      newCell = sentinelles_lines[index_line_R];
      l=0; /* index of B_indices */
      k=0;/* index of columns of R */
      while (k<nbParamSigma && l<n+3){
	if (B_indices[j][l] == k){
	  //k_th column of invA_times_Bj not empty
	  newCell->right = (struct cell *)malloc(sizeof(struct cell));
	  newCell = newCell->right;
	  newCell->value = invA_times_B[i][l];
	  newCell->col = k;
	  //printf("newCell->col = %d\n",newCell->col);
	  newCell->line = index_line_R;
	  newCell->down = NULL;
	  newCell->right = NULL;       
	  above[k]->down = newCell;  
	  above[k] = above[k]->down;   
	  k++;     
	  l++;
	}
	else /* if (k < B_indices[j][l])*/
	  //k_th column empty
	  k++; 
      }
      /* end of line (of R) */
      (*R)->lines[index_line_R] = sentinelles_lines[index_line_R]->right;
      free(sentinelles_lines[index_line_R]);
      sentinelles_lines[index_line_R] = NULL;
      index_line_R++; 
      
    }
    
  }
  //printf("avant R4\n");
  /* R4 */
  /* 4.1)    construction of the matrices R_4^0 and R_4^m*/
  /* 4.1.1)  construction of H_0 and H_m*/
  /* H_j (j=0,m) is represented by - a matrix H (in fact = B, cf. gradF1) with the non empty columns of H_j
     that we do not have to explicit
     - an array Hj_indices of the indices of the non empty columns of H_j   */
  /* the 2 last indices of the array Hj_indices are the indices of columns (1 0 ... 0)^T and (0 ... 0 1)^T */
  
  /* 4.1.1.1)  construction of H0_indices */
  for (i=0;i<n+1;i++)
    H0_indices[i] = (m+1)*(n+1) + 2*(m+1) + i;
  H0_indices[n+1] = (m+1)*(n+1) + 2*(m+1) + 2*(n+1);
  H0_indices[n+2] = (m+1)*(n+1) + 2*(m+1) + 2*(n+1) + 1;
  
  /* 4.1.1.2)  construction of Hm_indices */
  for (i=0;i<n+1;i++)
    Hm_indices[i] = (m+1)*(n+1) + 2*(m+1) + (n+1) + i;
  Hm_indices[n+1] = (m+1)*(n+1) + 2*(m+1) + 2*(n+1) + 2;
  Hm_indices[n+2] = (m+1)*(n+1) + 2*(m+1) + 2*(n+1) + 3;
  
  
  /* 4.1.2) construction of invA_times_H0 and invA_times_Hm */
  /* the non empty columns of invA_times_Hj are the same than the non empty columns of Hj */
  /* the indices of these columns are stocked in Hj_indices */
  /* invA_times_Hj (j=0,m) is represented by - a matrix  invA_times_H (in fact invA_times_B) with the non empty columns
     - an array Hj_indices of the indices of the non empty columns  */
  /* => invA_times_H0 and invA_times_Hm are already implicitly built! */
  
  /* 4.2) Fill in R_4^0 */
  /* R_4^0 = invA_times_H0 */
  for (line=0;line<n+1;line++){
    /* for each line of R_4^0 */
    newCell = sentinelles_lines[index_line_R];
    l = 0; /* index of H0_indices */
    k = 0;
    while (k<nbParamSigma && l<n+3){
      /* for each column of R (i.e of R_4^0),while there are still non empty columns left */
      if (k == H0_indices[l]){
	/* if column k of R_4^0 is not empty */	
	if (invA_times_B[line][l] != 0){
	  /* if the value is not zero */
	  newCell->right = (struct cell *)malloc(sizeof(struct cell));	
	  newCell = newCell->right;	 
	  newCell->value = invA_times_B[line][l];      	 
	  newCell->col = k;	  
	  newCell->line = index_line_R;	  
	  newCell->down = NULL; 	
	  newCell->right = NULL;       
	  above[k]->down = newCell;	 
	  above[k] = above[k]->down;   
	}
	l++;
	k++;
      }
      else /* if (k < H0_indices[l])*/ 
	k++;
    }
    /* there are no more non-empty columns in R_4^0 (we are at the end of H0_indices) */
    
    /* end of line */ 
    (*R)->lines[index_line_R] = sentinelles_lines[index_line_R]->right;
    free(sentinelles_lines[index_line_R]);
    index_line_R++; 

  }
  
  //printf("avant R4j\n");  

  if (m>1){
    /* 4.3) Construction of R_4^j, j=1,...,m-1 */
    
    /* 4.3.1) Construction of G */
    /* G = D without its first and its last line */
    tridiagToSparse(D,&G);
    
    /* 4.3.2) we compute and stock once and for all the matrices invC_times_Si, i=0,...,n */
    for (i=0;i<n+1;i++){
      /* Construction of Li */
      compute_Li(&Li,i,nbParamSigma,n,m,invA_times_B,B_indices);
      /* Construction of Pi */
      compute_Pi(&Pi,i,nbParamSigma,n,m,H0_indices,Hm_indices,invA_times_B);
      /* Computation of Si = G*Li-Pi */ 
      sparseTimesSparse(&G_times_Li,G,Li); 
      sparseMinusSparse(&Si,G_times_Li,Pi);   
      /* Computation of invC_times_Si */
      regTimesSparse(&(invC_times_S_tab[i]),invC,Si,m-1);
      /* free memory */
      desallocSparseMat_cells(Li);
      desallocSparseMat_cells(Pi);
      desallocSparseMat_cells(G_times_Li);
      desallocSparseMat_cells(Si);
      
      free(G_times_Li->lines);
      free(G_times_Li->columns);
      free(G_times_Li);
      
      free(Si->lines);
      free(Si->columns);
      free(Si);
    }
  
  }
  
  /* 4.3.3) compute and fill in R_4^j, j=1,...,M-1 */
  for (j=1;j<m;j++){
    /* construction of the matrix R_4^j */
    for (i=0;i<n+1;i++){
      /* ith line of R_4^j = jth line of invC_times_Si */
      newCell = sentinelles_lines[index_line_R];
      /* ptCell points on the jth line of invC_times_Si */
      ptCell = invC_times_S_tab[i]->lines[j-1];
      while (ptCell != NULL){  /* while the line is not finished */
	/* fill in newCell */
	newCell->right = (struct cell *)malloc(sizeof(struct cell));
	newCell = newCell->right;
	newCell->line = index_line_R;
	newCell->col = ptCell->col;
	newCell->value = ptCell->value;
	newCell->down = NULL;
	newCell->right = NULL;
	/* update above */
	above[ptCell->col]->down = newCell;
	above[ptCell->col] =  above[ptCell->col]->down;
	ptCell = ptCell->right;
      }
         /* end of the line */
      (*R)->lines[index_line_R] = sentinelles_lines[index_line_R]->right;
      free(sentinelles_lines[index_line_R]);
      /* goes to next line of R */
      index_line_R++;
    }
  }
 
  //printf("avant R4m\n");
  /* 4.4) Fill in R_4^m */
  /* R_4^m = invA_times_Hm */
  for (line=0;line<n+1;line++){
    /* for each line of R_4^m */
    newCell = sentinelles_lines[index_line_R];
    l = 0; /* index of H0_indices */
    k = 0;
    while (k<nbParamSigma && l<n+3){
      /* for each column of R (i.e of R_4^m) */
      if (k == Hm_indices[l]){
	/* if column k of R_4^m is not empty */	
	if (invA_times_B[line][l] != 0){
	  /* if value is not zero */
	  newCell->right = (struct cell *)malloc(sizeof(struct cell));
	  newCell = newCell->right;
	  newCell->value = invA_times_B[line][l];
	  newCell->col = k;
	  newCell->line = index_line_R;
	  newCell->down = NULL;
	  newCell->right = NULL;
	  above[k]->down = newCell;
	  above[k] = above[k]->down;
	}	 
	k++;
	l++;
      }

      else /* if (k<Hm_indices[l])*/
	l++;
    }
    
    /* end of line */
    (*R)->lines[index_line_R] = sentinelles_lines[index_line_R]->right;
    free(sentinelles_lines[index_line_R]);
    index_line_R++; 
    
  }
  
  for (j=0;j<nbParamSigma;j++){
    (*R)->columns[j] = sentinelles_columns[j]->down;
    free(sentinelles_columns[j]);
  }

  /* memory desallocation */
  free(sentinelles_lines);
  free(sentinelles_columns);

  free(H0_indices);
  free(Hm_indices);

  if (m>1){
    desallocSparseMat_cells(G);
    free(G->lines);
    free(G->columns);
    free(G);

    free(Li->lines);
    free(Li->columns);
    free(Li);
    
    free(Pi->lines);
    free(Pi->columns);
    free(Pi); 
  }

  free(above);

  if (m>1){
    for (i=0;i<n+1;i++){
      desallocSparseMat_cells(invC_times_S_tab[i]);
      free(invC_times_S_tab[i]->lines);
      free(invC_times_S_tab[i]->columns);
      free(invC_times_S_tab[i]);
    }
    free(invC_times_S_tab);

    for (k=0;k<m-1;k++)
      free(invC[k]);
    free(invC);
    
    free(C->subdiag);
    free(C->diag);
    free(C->updiag);
    free(C);
  }    

}



void compute_Li(struct sparseMat **pt_Li, int i, int nbParamSigma, int n, int m, double **invA_times_B, int **B_indices){
  

  struct cell **pt_above,**sentinelles_lines,**sentinelles_columns;
  struct cell *newCell;
  int j,k,l;

  sentinelles_columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc((m+1)*sizeof(struct cell *));
  
  for (j=0;j<nbParamSigma;j++){
    sentinelles_columns[j] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[j]->down = NULL;
    sentinelles_columns[j]->right = NULL;
  }
  for (j=0;j<m+1;j++){
    sentinelles_lines[j] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[j]->right = NULL;
    sentinelles_lines[j]->down = NULL;
  }
  
  /* memory allocation and initialization of pt_above */
  pt_above = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));

  for (j=0;j<nbParamSigma;j++)
    pt_above[j] = sentinelles_columns[j];


  for (j=0;j<m+1;j++){
    /* j_th line of Li = ith line of invA_times_B_j */
    newCell = sentinelles_lines[j];
    l = 0; /* index of B_indices[ ][j] */
    k = 0; /* index of columns of Li */
    while (k<nbParamSigma  && l<n+3){
      if (k == B_indices[j][l]){
	if (invA_times_B[i][l] != 0){
	  /* creation of a new cell */
	  newCell->right = (struct cell *) malloc(sizeof(struct cell));
	  newCell = newCell->right;	  
	  newCell->line = j;
	  newCell->col = k;
	  newCell->right = NULL;
	  newCell->down = NULL;
	  newCell->value = invA_times_B[i][l];
	  /* update of pt_above */
	  pt_above[k]->down = newCell;
	  pt_above[k] = pt_above[k]->down;
	}
	k++;
	l++;
      }
      else //if (k<B_indices[l][j])
	k++;
    }

    /* end of line */ 
    (*pt_Li)->lines[j] = sentinelles_lines[j]->right;
    free(sentinelles_lines[j]);
    sentinelles_lines[j] = NULL;
  }


  for (j=0;j<nbParamSigma;j++){
    (*pt_Li)->columns[j] = sentinelles_columns[j]->down;
    free(sentinelles_columns[j]);
    sentinelles_columns[j] = NULL;
    pt_above[j] = NULL;
  }

  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(pt_above);
  pt_above = NULL;
  
}



void compute_Pi(struct sparseMat **pt_Pi, int i, int nbParamSigma, int n, int m, int *H0_indices, int *Hm_indices, double **invA_times_B){
  
  

  struct cell **pt_above,**sentinelles_lines,**sentinelles_columns;
  struct cell *newCell;
  int j,k,l;

  sentinelles_columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc((m-1)*sizeof(struct cell *));
  
  for (j=0;j<nbParamSigma;j++){
    sentinelles_columns[j] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[j]->down = NULL;
    sentinelles_columns[j]->right = NULL;
  }
  for (j=0;j<m-1;j++){
    sentinelles_lines[j] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[j]->right = NULL;
    sentinelles_lines[j]->down = NULL;
  }
  
  
  /* memory allocation and initialization of pt_above */
  pt_above = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  for (j=0;j<nbParamSigma;j++)
    pt_above[j] = sentinelles_columns[j];
  

  /* 0_th line of Pi = ith line of invA_times_H_0 */
  newCell = sentinelles_lines[0];
  l = 0; /* index of H0_indices[ ] */
  k = 0; /* index of columns of Pi */
  while (k<nbParamSigma  && l<n+3){
    if (k == H0_indices[l]){
      /* if the column k is not empty */
      if (invA_times_B[i][l] != 0){
	/* if the value is not zero */
	/* creation of a new cell */
	newCell->right = (struct cell *) malloc(sizeof(struct cell));
	newCell = newCell->right;	  
	newCell->line = j;
	newCell->col = k;
	newCell->right = NULL;
	newCell->down = NULL;
	newCell->value = invA_times_B[i][l];
	/* update of pt_above */
	pt_above[k]->down = newCell;
	pt_above[k] = pt_above[k]->down;
      }
      k++;
      l++;
    }
    else //if (k<H0_indices[l])
      k++;
  }

  /* end of line */ 
  (*pt_Pi)->lines[0] = sentinelles_lines[0]->right;
  free(sentinelles_lines[0]);
  sentinelles_lines[0] = NULL;

  /* lines 1 to m-3 of Pi are empty */
  for (j=1;j<m-2;j++){
    (*pt_Pi)->lines[j] = NULL;
    free(sentinelles_lines[j]);  
    sentinelles_lines[j] = NULL;
  }

  /* (m-2)_th line of Pi = ith line of invA_times_H_M */  
  newCell = sentinelles_lines[m-2];
  l = 0; /* index of H0_indices[ ] */
  k = 0; /* index of columns of Pi */
  while (k<nbParamSigma  && l<n+3){
    if (k == Hm_indices[l]){
      /* if the column k is not empty */
      if (invA_times_B[i][l] != 0){
	/* if the value is not zero */
	/* creation of a new cell */
	newCell->right = (struct cell *) malloc(sizeof(struct cell));
	newCell = newCell->right;	  
	newCell->line = j;
	newCell->col = k;
	newCell->right = NULL;
	newCell->down = NULL;
	newCell->value = invA_times_B[i][l];
	/* update of pt_above */
	pt_above[k]->down = newCell;
	pt_above[k] = pt_above[k]->down;
      }
      k++;
      l++;
    }
    else //if (k<Hm_indices[l])
      k++;
  }

  /* end of line */ 
  (*pt_Pi)->lines[m-2] = sentinelles_lines[m-2]->right;
  free(sentinelles_lines[m-2]);  
  sentinelles_lines[m-2] = NULL;

  for (j=0;j<nbParamSigma;j++){
    (*pt_Pi)->columns[j] = sentinelles_columns[j]->down;
    free(sentinelles_columns[j]);
    sentinelles_columns[j] = NULL;
    pt_above[j] = NULL;
  }

  free(sentinelles_lines);
  sentinelles_lines = NULL;
  free(sentinelles_columns);
  sentinelles_columns = NULL;
  free(pt_above);
  pt_above = NULL;


}



void compute_deriv_Au_j(struct sparseMat **pt_deriv_Au_j, int j, double *y_coarseGrid, double *T_coarseGrid, double *y_fineGrid, double *T_fineGrid, int N, int M, int n, int m, struct sparseMat *R, double **sigmaFineGrid, double **sigmaCoarseGrid, double **U){

  struct cell **sentinelles_lines,**sentinelles_columns;
  struct cell **pt_above;
  int i,k,l,p,s,line,col;
  struct sparseMat ***transp_Vkl;
  struct sparseMat *sum,*tmp_sum,*deriv_sigma;
  struct sparseMat *deriv_alpha,*deriv_beta,*deriv_gamma,*line_i,*temp;
  double coeff,h_i,h_i_1;
  struct cell *pt_line1, *pt_line2;
  int nbParamSigma;

  
  nbParamSigma = (n+1)*(m+1)+2*(n+1)+2*(m+1)+4;

  /* memory allocation for the sentinelles */
  sentinelles_columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  sentinelles_lines = (struct cell **) malloc((N-1)*sizeof(struct cell *));
  
  /* memory allocation and initialization of pt_above_res */
  pt_above = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  
  /* alloc and init sentinelles, who will help build (*pt_deriv_Au_j) */
  for (col=0;col<nbParamSigma;col++){
    sentinelles_columns[col] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_columns[col]->down = NULL;
    sentinelles_columns[col]->right = NULL;
    pt_above[col] = sentinelles_columns[col];
  }
  for (line=0;line<N-1;line++){
    sentinelles_lines[line] = (struct cell *) malloc(sizeof(struct cell));
    sentinelles_lines[line]->right = NULL;
    sentinelles_lines[line]->down = NULL;
  }  
  
  sum = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  sum->nbLines = 1;
  sum->nbCol = 4*(n+1)*(m+1);
  sum->columns = (struct cell **) malloc(4*(n+1)*(m+1)*sizeof(struct cell *));
  sum->lines = (struct cell **) malloc(sizeof(struct cell *));


/* memory allocation of transp_Vkl */
  transp_Vkl = (struct sparseMat ***) malloc(4*sizeof(struct sparseMat **));
  for (p=0;p<4;p++){
    transp_Vkl[p] = (struct sparseMat **) malloc(4*sizeof(struct sparseMat *));
    for (s=0;s<4;s++){
      transp_Vkl[p][s] = (struct sparseMat *) malloc(sizeof(struct sparseMat));
      transp_Vkl[p][s]->columns = (struct cell **) malloc(4*(n+1)*(m+1)*sizeof(struct cell *));
      for (i=0;i<4*(n+1)*(m+1);i++)
	transp_Vkl[p][s]->columns[i] = NULL;
      transp_Vkl[p][s]->lines = (struct cell **) malloc(sizeof(struct cell *)); // Vkl = column vector => transp_Vkl has 1 line
      transp_Vkl[p][s]->lines[0] = NULL;
      transp_Vkl[p][s]->nbCol = 4*(n+1)*(m+1);
      transp_Vkl[p][s]->nbLines = 1;
    }
  }



  /* memory allocation for the matrix deriv_alpha */
  deriv_alpha = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  deriv_alpha->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  deriv_alpha->lines = (struct cell **) malloc(sizeof(struct cell *));
  deriv_alpha->nbCol = nbParamSigma;
  deriv_alpha->nbLines = 1; 


  /* memory allocation for the matrix deriv_beta */
  deriv_beta = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  deriv_beta->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  deriv_beta->lines = (struct cell **) malloc(sizeof(struct cell *));
  deriv_beta->nbCol = nbParamSigma;
  deriv_beta->nbLines = 1; 

  /* memory allocation for the matrix deriv_gamma */
  deriv_gamma = (struct sparseMat *) malloc(sizeof(struct sparseMat));
  deriv_gamma->columns = (struct cell **) malloc(nbParamSigma*sizeof(struct cell *));
  deriv_gamma->lines = (struct cell **) malloc(sizeof(struct cell *));
  deriv_gamma->nbCol = nbParamSigma;
  deriv_gamma->nbLines = 1; 


  /* initialization of sum */
  sum->lines[0] = NULL;
  for (k=0;k<4*(n+1)*(m+1);k++)
    sum->columns[k] = NULL;
  
  
  for (i=1;i<N;i++){
    printf("i=%d,j=%d\n",i,j);
    /* for each line of (*pt_deriv_Au_j) */
    /* 1) compute deriv_sigma = \frac{\partial \sigma(y_i,T_j)}{\partial \Sigma}, sparse matrix of size 1*nbParamSigma */
    /* 1.1) find k,l so that (y_i,T_j) of the fine grid is in R_kl of the coarse grid */
    k = find_index(y_fineGrid[i],y_coarseGrid,n,0); /* y_coarseGrid[k] <= y_fineGrid[i] < y_coarseGrid[k+1] */
    l = find_index(T_fineGrid[j],T_coarseGrid,m,0); /* T_coarseGrid[l] <= T_fineGrid[j] < T_coarseGrid[l+1] */
    /* 1.2) compute V_kl */
    //printf("k=%d\n",k);
    //printf("l=%d\n",l);
    compute_transp_Vkl(transp_Vkl,k,l,n,m,y_coarseGrid[k+1]-y_coarseGrid[k],T_coarseGrid[l+1]-T_coarseGrid[l]);

  
    /* 1.3) compute the sum */
    for (p=0;p<4;p++)
      for (s=0;s<4;s++){
	//printf("transp_Vkl[%d][%d] :\n",p,s);
	//printSparseMat(transp_Vkl[p][s]);
	scalarTimesSparse((y_fineGrid[i]-y_coarseGrid[k])*(T_fineGrid[j]-T_coarseGrid[l]),&transp_Vkl[p][s]);
	/* tmp_sum = sum + transp_Vkl[p][s] */
	sparsePlusSparse(&tmp_sum,sum,transp_Vkl[p][s]);
	//printf("k=%d,l=%d, transp_Vkl[0][0] = \n",k,l);
	//printSparseMat(transp_Vkl[0][0]);
	//printf("transp_Vkl[%d][%d] :\n",p,s);
	//printSparseMat(transp_Vkl[p][s]);
	desallocSparseMat_cells(sum);
	/* sum = tmp_sum */
	affectSparse(&sum,tmp_sum);

	//if (isEmpty(sum))
	printf("transp_Vkl[%d][%d] :\n",p,s);
	printSparseMat(transp_Vkl[p][s]);
	printf("sum:\n");
	printSparseMat(sum);

	desallocSparseMat_cells(tmp_sum);
	free(tmp_sum->columns);
	free(tmp_sum->lines);
	free(tmp_sum);
	desallocSparseMat_cells(transp_Vkl[p][s]);
      }


    //printf("sum:\n");
    //printSparseMat(sum);

    /*1.4) multiply sum by R => deriv_sigma*/
    sparseTimesSparse(&deriv_sigma,sum,R);
    
    //printf("deriv_sigma:\n");
    //printSparseMat(deriv_sigma);


    /* 2) computation of 
       - deriv_alpha = \frac{\partial \alpha_{i,T_j}}{\partial \Sigma} of size 1,nbParamSigma 
       - deriv_beta = \frac{\partial \beta_{i,T_j}}{\partial \Sigma} of size 1,nbParamSigma 
       - deriv_gamma = \frac{\partial \gamma_{i,T_j}}{\partial \Sigma} of size 1,nbParamSigma  */
    h_i = y_fineGrid[i+1]-y_fineGrid[i];
    h_i_1 = y_fineGrid[i]-y_fineGrid[i-1];
    if (i != 1){
      affectSparse(&deriv_alpha,deriv_sigma);
      coeff = (-2/((h_i+h_i_1)*h_i_1)-1/(2*h_i_1))*sigmaFineGrid[i][j];
      scalarTimesSparse(coeff,&deriv_alpha);
      //printf("deriv_alpha:\n");
      //printSparseMat(deriv_alpha);
    }
    affectSparse(&deriv_beta,deriv_sigma);
    coeff = (2/(h_i+h_i_1)*(1/h_i+1/h_i_1)+0.5*(1/h_i_1-1/h_i))*sigmaFineGrid[i][j];
    scalarTimesSparse(coeff,&deriv_beta);
    //printf("deriv_beta:\n");
    //printSparseMat(deriv_beta);
    if (i != N-1){
      affectSparse(&deriv_gamma,deriv_sigma);
      coeff = (-2/((h_i+h_i_1)*h_i)+1/(2*h_i))*sigmaFineGrid[i][j];
      scalarTimesSparse(coeff,&deriv_gamma);
      //printf("deriv_gamma:\n");
      //printSparseMat(deriv_gamma);
    }
    


    /* 3) computation of line_i = (i-1)_th line of (*pt_deriv_Au_j) */
    if (i == 1){
      scalarTimesSparse(U[0][j],&deriv_beta);
      scalarTimesSparse(U[1][j],&deriv_gamma);
      sparsePlusSparse(&line_i,deriv_beta,deriv_gamma);
    }
    else if (i == N-1){
      scalarTimesSparse(U[N-3][j],&deriv_alpha);
      scalarTimesSparse(U[N-2][j],&deriv_beta);
      sparsePlusSparse(&line_i,deriv_alpha,deriv_beta);
    }
    else{
      scalarTimesSparse(U[i-2][j],&deriv_alpha);
      //printf("deriv_alpha:\n");
      //printSparseMat(deriv_alpha);
      scalarTimesSparse(U[i-1][j],&deriv_beta);
      //printf("deriv_beta:\n");
      //printSparseMat(deriv_beta);
      scalarTimesSparse(U[i][j],&deriv_gamma);
      sparsePlusSparse(&temp,deriv_alpha,deriv_beta);
      //printf("temp:\n");
      //printSparseMat(temp);
      sparsePlusSparse(&line_i,temp,deriv_gamma);
    }
    
    /*4) stock line_i in the (i-1)_th line of (*pt_deriv_Au_j) */
    pt_line1 = line_i->lines[0];
    pt_line2 = sentinelles_lines[i-1];
    while (pt_line1 != NULL){
      pt_line2->right = (struct cell *) malloc(sizeof(struct cell));
      pt_line2 = pt_line2->right;
      pt_line2->value = pt_line1->value;
      pt_line2->line = i-1;
      pt_line2->col = pt_line1->col;
      pt_line2->down = NULL;
      pt_line2->right = NULL;
      pt_above[pt_line1->col]->down = pt_line2;
      pt_above[pt_line1->col] = pt_above[pt_line1->col]->down;
      pt_line1 = pt_line1->right;
    }
    
    /* end of line i-1 */
    (*pt_deriv_Au_j)->lines[i-1] = sentinelles_lines[i-1]->right;
    free(sentinelles_lines[i-1]);
    
    
    /* free memory */
    
    desallocSparseMat_cells(sum);
    
    desallocSparseMat_cells(deriv_sigma);
    free(deriv_sigma->columns);
    free(deriv_sigma->lines);
    free(deriv_sigma);
    
    if (i != 1)
      desallocSparseMat_cells(deriv_alpha);
    desallocSparseMat_cells(deriv_beta);
    
    if (i != N-1)
      desallocSparseMat_cells(deriv_gamma);
    
    desallocSparseMat_cells(line_i);
    free(line_i->lines);
    free(line_i->columns);
    free(line_i);
    
    if (i!=1 && i!=N-1){
      desallocSparseMat_cells(temp);
      free(temp->lines);
      free(temp->columns);
      free(temp);
    }
    
  }

  for (col=0;col<nbParamSigma;col++){
    (*pt_deriv_Au_j)->columns[col] = sentinelles_columns[col]->down;
    free(sentinelles_columns[col]);
  }
  
  free(sentinelles_lines);
  free(sentinelles_columns);
  free(pt_above);
  
  free(deriv_alpha->columns);
  free(deriv_alpha->lines);
  free(deriv_alpha);

  free(deriv_beta->columns);
  free(deriv_beta->lines);
  free(deriv_beta);

  free(deriv_gamma->columns);
  free(deriv_gamma->lines);
  free(deriv_gamma);

  free(sum->lines);
  free(sum->columns);
  free(sum);
  
  for (p=0;p<4;p++){
    for (s=0;s<4;s++){
      free(transp_Vkl[p][s]->lines);
      free(transp_Vkl[p][s]->columns);
      free(transp_Vkl[p][s]);  
    }
    free(transp_Vkl[p]);
  }
  free(transp_Vkl);
  

}



void compute_transp_Vkl(struct sparseMat ***transp_Vkl, int k, int l, int n, int m, double delta_y, double delta_T){
  /* transp_Vkl is a 4*4 array of pointers on V_ij,kl */

  struct cell *newCell;
 

  /************* transp_Vkl[0][0] ******************************/

  transp_Vkl[0][0]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[0][0]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1;
  transp_Vkl[0][0]->columns[l*(n+1)+k] = newCell;
  
  newCell->right = NULL;

  //printf("k=%d,l=%d, transp_Vkl[0][0] = \n",k,l);
  //printSparseMat(transp_Vkl[0][0]);


 /************* transp_Vkl[0][1] ******************************/

  transp_Vkl[0][1]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[0][1]->lines[0];
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1;
  transp_Vkl[0][1]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;
  
  newCell->right = NULL;
  
  /************* transp_Vkl[0][2] ******************************/


  transp_Vkl[0][2]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[0][2]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_T,2));
  transp_Vkl[0][2]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_T,2));
  transp_Vkl[0][2]->columns[(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/delta_T;
  transp_Vkl[0][2]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/delta_T;
  transp_Vkl[0][2]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = NULL;



 /************* transp_Vkl[0][3] ******************************/

  
  transp_Vkl[0][3]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[0][3]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_T,3));
  transp_Vkl[0][3]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_T,3));
  transp_Vkl[0][3]->columns[(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_T,2));
  transp_Vkl[0][3]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_T,2));
  transp_Vkl[0][3]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = NULL;



  /************* transp_Vkl[1][0] ******************************/

  transp_Vkl[1][0]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[1][0]->lines[0];
  newCell->col = 2*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1;
  transp_Vkl[1][0]->columns[2*(n+1)*(m+1)+l*(n+1)+k] = newCell;
  
  newCell->right = NULL;
  


  /************* transp_Vkl[1][1] ******************************/

  transp_Vkl[1][1]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[1][1]->lines[0];
  newCell->col = 3*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1;
  transp_Vkl[1][1]->columns[3*(n+1)*(m+1)+l*(n+1)+k] = newCell;
  
  newCell->right = NULL;



   /************* transp_Vkl[1][2] ******************************/

  transp_Vkl[1][2]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[1][2]->lines[0];
  newCell->col = 2*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_T,2));
  transp_Vkl[1][2]->columns[2*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_T,2));
  transp_Vkl[1][2]->columns[2*(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/delta_T;
  transp_Vkl[1][2]->columns[3*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/delta_T;
  transp_Vkl[1][2]->columns[3*(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = NULL;


 /************* transp_Vkl[1][3] ******************************/

  
  transp_Vkl[1][3]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[1][3]->lines[0];
  newCell->col = 2*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_T,3));
  transp_Vkl[1][3]->columns[2*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_T,3));
  transp_Vkl[1][3]->columns[2*(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_T,2));
  transp_Vkl[1][3]->columns[3*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_T,2));
  transp_Vkl[1][3]->columns[3*(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = NULL;


 /************* transp_Vkl[2][0] ******************************/

  transp_Vkl[2][0]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[2][0]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_y,2));
  transp_Vkl[2][0]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_y,2));
  transp_Vkl[2][0]->columns[l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/delta_y;
  transp_Vkl[2][0]->columns[2*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/delta_y;
  transp_Vkl[2][0]->columns[2*(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = NULL;


/************* transp_Vkl[2][1] ******************************/

  transp_Vkl[2][1]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[2][1]->lines[0];
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_y,2));
  transp_Vkl[2][1]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_y,2));
  transp_Vkl[2][1]->columns[(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/delta_y;
  transp_Vkl[2][1]->columns[3*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/delta_y;
  transp_Vkl[2][1]->columns[3*(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = NULL;


  /************* transp_Vkl[2][2] ******************************/
  
  transp_Vkl[2][2]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[2][2]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 9/(pow(delta_T,2)*pow(delta_y,2));
  transp_Vkl[2][2]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 9/(pow(delta_T,2)*pow(delta_y,2));
  transp_Vkl[2][2]->columns[l*(n+1)+k+1] = newCell;
  
  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -9/(pow(delta_T,2)*pow(delta_y,2));
  transp_Vkl[2][2]->columns[(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 9/(pow(delta_T,2)*pow(delta_y,2));
  transp_Vkl[2][2]->columns[(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 6/(pow(delta_y,2)*delta_T);
  transp_Vkl[2][2]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -6/(pow(delta_y,2)*delta_T);
  transp_Vkl[2][2]->columns[(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_y,2)*delta_T);
  transp_Vkl[2][2]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_y,2)*delta_T);
  transp_Vkl[2][2]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 6/(delta_y*pow(delta_T,2));
  transp_Vkl[2][2]->columns[2*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(delta_y*pow(delta_T,2));
  transp_Vkl[2][2]->columns[2*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -6/(delta_y*pow(delta_T,2));
  transp_Vkl[2][2]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(delta_y*pow(delta_T,2));
  transp_Vkl[2][2]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 4/(delta_T*delta_y);
  transp_Vkl[2][2]->columns[3*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(delta_T*delta_y);
  transp_Vkl[2][2]->columns[3*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(delta_T*delta_y);
  transp_Vkl[2][2]->columns[ 3*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(delta_T*delta_y);
  transp_Vkl[2][2]->columns[3*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = NULL;

  /************* transp_Vkl[2][3] ******************************/

  transp_Vkl[2][3]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[2][3]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -6/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[2][3]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 6/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[2][3]->columns[l*(n+1)+k+1] = newCell;
  
  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 6/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[2][3]->columns[(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -6/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[2][3]->columns[(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[2][3]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[2][3]->columns[(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[2][3]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[2][3]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -4/(delta_y*pow(delta_T,3));
  transp_Vkl[2][3]->columns[2*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(delta_y*pow(delta_T,3));
  transp_Vkl[2][3]->columns[2*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -4/(delta_y*pow(delta_T,3));
  transp_Vkl[2][3]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(delta_y*pow(delta_T,3));
  transp_Vkl[2][3]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(delta_y*pow(delta_T,2));
  transp_Vkl[2][3]->columns[3*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/(delta_y*pow(delta_T,2));
  transp_Vkl[2][3]->columns[3*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(delta_y*pow(delta_T,2));
  transp_Vkl[2][3]->columns[ 3*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/(delta_y*pow(delta_T,2));
  transp_Vkl[2][3]->columns[3*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;
  

  newCell->right = NULL;
  
  /************* transp_Vkl[3][0] ******************************/


 transp_Vkl[3][0]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[3][0]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_y,3));
  transp_Vkl[3][0]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_y,3));
  transp_Vkl[3][0]->columns[l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2));
  transp_Vkl[3][0]->columns[2*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2));
  transp_Vkl[3][0]->columns[2*(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = NULL;




  /************* transp_Vkl[3][1] ******************************/


 transp_Vkl[3][1]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[3][1]->lines[0];
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_y,3));
  transp_Vkl[3][1]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_y,3));
  transp_Vkl[3][1]->columns[(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2));
  transp_Vkl[3][1]->columns[3*(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2));
  transp_Vkl[3][1]->columns[3*(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = NULL;



  /************* transp_Vkl[3][2] ******************************/


  transp_Vkl[3][2]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[3][2]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -6/(pow(delta_T,2)*pow(delta_y,3));
  transp_Vkl[3][2]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 6/(pow(delta_T,2)*pow(delta_y,3));
  transp_Vkl[3][2]->columns[l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 6/(pow(delta_T,2)*pow(delta_y,3));
  transp_Vkl[3][2]->columns[(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -6/(pow(delta_T,2)*pow(delta_y,3));
  transp_Vkl[3][2]->columns[(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -4/(pow(delta_y,3)*delta_T);
  transp_Vkl[3][2]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 4/(pow(delta_y,3)*delta_T);
  transp_Vkl[3][2]->columns[(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_y,3)*delta_T);
  transp_Vkl[3][2]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

 newCell->right = (struct cell *) malloc(sizeof(struct cell)); 
  newCell = newCell->right; 
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k+1; 
  newCell->line = 0;
  newCell->down = NULL; 
  newCell->value = 2/(pow(delta_y,3)*delta_T);  
  transp_Vkl[3][2]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k+1] = newCell;  

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][2]->columns[2*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][2]->columns[2*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][2]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 3/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][2]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(delta_T*pow(delta_y,2));
  transp_Vkl[3][2]->columns[3*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(delta_T*pow(delta_y,2));
  transp_Vkl[3][2]->columns[3*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/(delta_T*pow(delta_y,2));
  transp_Vkl[3][2]->columns[ 3*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -1/(delta_T*pow(delta_y,2));
  transp_Vkl[3][2]->columns[3*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = NULL;

  /************* transp_Vkl[3][3] ******************************/

  transp_Vkl[3][3]->lines[0] = (struct cell *) malloc(sizeof(struct cell));
  newCell = transp_Vkl[3][3]->lines[0];
  newCell->col = l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 4/(pow(delta_y,3)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -4/(pow(delta_y,3)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[l*(n+1)+k+1] = newCell;
  
  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -4/(pow(delta_y,3)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 4/(pow(delta_y,3)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_y,3)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[(n+1)*(m+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_y,3)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[(n+1)*(m+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_y,3)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = (n+1)*(m+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_y,3)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[(n+1)*(m+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[2*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 2/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[2*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 2*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = -2/(pow(delta_y,2)*pow(delta_T,3));
  transp_Vkl[3][3]->columns[2*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+l*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[3*(m+1)*(n+1)+l*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+l*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[3*(m+1)*(n+1)+l*(n+1)+k+1] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col =  3*(m+1)*(n+1)+(l+1)*(n+1)+k;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[3*(m+1)*(n+1)+(l+1)*(n+1)+k] = newCell;

  newCell->right = (struct cell *) malloc(sizeof(struct cell));
  newCell = newCell->right;
  newCell->col = 3*(m+1)*(n+1)+(l+1)*(n+1)+k+1;
  newCell->line = 0;
  newCell->down = NULL;
  newCell->value = 1/(pow(delta_y,2)*pow(delta_T,2));
  transp_Vkl[3][3]->columns[3*(m+1)*(n+1)+(l+1)*(n+1)+k+1] = newCell;

  newCell->right = NULL;

 
}




