#include "solveSystem.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

void tridiagTimesVect(double *res, struct tridiag *A, double *vect){
/* res = A*vect */
  
  int i;
  
  res[0] = A->diag[0]*vect[0] + A->updiag[0]*vect[1];
  for (i=1;i<A->size-1;i++)
    res[i] = A->subdiag[i-1]*vect[i-1] + A->diag[i]*vect[i] + A->updiag[i]*vect[i+1];
  res[A->size-1] = A->subdiag[A->size-2]*vect[A->size-2] + A->diag[A->size-1]*vect[A->size-1];

}

/*changes the tridiagonal system to a triangular system*/
void tridiagToBidiagSyst(struct bidiagSystem *U, struct tridiagSystem *S){

  int size,i;

  size = S->size;
  U->T->diag[size-1] = S->T->diag[size-1];
  U->b[size-1] = S->b[size-1];


  for (i=size-2;i>=0;i--){
    U->T->diag[i] = S->T->diag[i] - S->T->updiag[i]*S->T->subdiag[i]/(U->T->diag[i+1]);
    U->b[i] = S->b[i] - S->T->updiag[i]*U->b[i+1]/(U->T->diag[i+1]);
    U->T->subdiag[i] = S->T->subdiag[i];
  }


}

/*changes a tridiagonal matrix to a bidiagonal one*/
void tridiagToBidiagMat(struct bidiag *B, struct tridiag *T){

  int size,i;

  size = T->size;
  B->diag[size-1] = T->diag[size-1];

  for (i=size-2;i>=0;i--){
    B->diag[i] = T->diag[i] - T->updiag[i]*T->subdiag[i]/(B->diag[i+1]);
    B->subdiag[i] = T->subdiag[i];
  }

}


/*solves the bidiagonal system at each time step*/
void solveSyst(double *sol, struct bidiagSystem *S){

  int size,i;

  size = S->size;
  sol[0] = (S->b)[0]/(S->T->diag[0]);
  
  for (i=1;i<=size-1;i++)
    sol[i] = (S->b[i] - S->T->subdiag[i-1]*sol[i-1])/(S->T->diag[i]);
  
} 


void tridiagMatrixInv(double **inv, struct tridiag *T){

  struct bidiagSystem *U;
  struct tridiagSystem *V;
  double *x_j;
  int i,j,size;

  size = T->size;

  /*memory allocation */
 
  /* memory allocation for the tridiagonal system V*/
  V = (struct tridiagSystem *) malloc(sizeof(struct tridiagSystem));
  V->T = (struct tridiag *) malloc(sizeof(struct tridiag));
  V->T->size = size;
  V->T->subdiag = (double *) malloc((size-1)*sizeof(double));
  V->T->diag = (double *) malloc(size*sizeof(double));
  V->T->updiag = (double *) malloc((size-1)*sizeof(double));
  V->b = (double *) malloc(size*sizeof(double));
  V->size = size;

  /* memory allocation for the bidiagonal system U*/
  U = (struct bidiagSystem *) malloc(sizeof(struct bidiagSystem));
  U->T = (struct bidiag *) malloc(sizeof(struct bidiag));
  U->T->subdiag = (double *) malloc((size-1)*sizeof(double));
  U->T->diag = (double *) malloc(size*sizeof(double));
  U->T->size = size;
  U->b = (double *) malloc(size*sizeof(double));
  U->size = size;

  /*memory alloc for x_j*/
  x_j = (double *) malloc(size*sizeof(double));



  /*initialization of V->T*/
  affectOperator(T,V->T);


  for (j=0;j<size;j++){
    /*initialisation of V->b */
    for(i=0;i<j;i++)
      V->b[i] = 0;
    V->b[j] = 1;
    for(i=j+1;i<size;i++)
      V->b[i] = 0;
    /*computation of the jth column of the inverse */
    tridiagToBidiagSyst(U,V);
    solveSyst(x_j,U);
    for (i=0;i<size;i++){
      inv[i][j] = x_j[i];
    }
  }
    
  /* memory desallocation */
  free(V->T->subdiag);
  V->T->subdiag = NULL;
  free(V->T->diag);  
  V->T->diag = NULL;
  free(V->T->updiag);
  V->T->subdiag = NULL;
  free(V->T);
  V->T = NULL;
  free(V->b);
  V->b = NULL;
  free(V);
  V = NULL;

  free(U->T->subdiag);
  U->T->subdiag = NULL;
  free(U->T->diag);  
  U->T->diag = NULL;
  free(U->T);
  U->T = NULL;
  free(U->b);
  U->b = NULL;
  free(U);
  U = NULL;

  free(x_j);
  x_j = NULL;

}


void tridiagMatMult(double **res, struct tridiag *A, double **B){
/* res = B*A */
  int i,j,size;

  size = A->size;
  for (i=0;i<size;i++){
    res[i][0] = A->diag[0]*B[i][0] + A->subdiag[0]*B[i][1];
    for (j=1;j<size-1;j++)
      res[i][j] = A->updiag[j-1]*B[i][j-1] + A->diag[j]*B[i][j] + A->subdiag[j]*B[i][j+1];
    res[i][size-1] = A->updiag[size-2]*B[i][size-2] + A->diag[size-1]*B[i][size-1];
  }

}

/*equals the tridiag matrix B to the tridiag matrix A when they already have the same size*/
void affectOperator(struct tridiag *A, struct tridiag *B){

  int i;

  for (i=0;i<=B->size-2;i++){
    B->subdiag[i] = A->subdiag[i];
    B->diag[i] = A->diag[i];
    B->updiag[i] = A->updiag[i];
  }
  B->diag[B->size-1] = A->diag[B->size-1];

}


void printBidiagSyst(struct bidiagSystem *B){

  int i;

  printf("\n\nsize = %d\n",B->size);
  printf("subdiag : ");
  for (i=0;i<B->size-1;i++)
    printf("%lf  ",B->T->subdiag[i]);

  printf("\ndiag : ");
  for (i=0;i<B->size;i++)
    printf("%lf  ",B->T->diag[i]);

  printf("\nb : ");
  for (i=0;i<B->size;i++)
    printf("%lf  ",B->b[i]);

}



void printTridiagSyst(struct tridiagSystem *B){

  int i;

  printf("\n\nsize = %d\n",B->size);
  printf("subdiag : ");
  for (i=0;i<B->size-1;i++)
    printf("%lf  ",B->T->subdiag[i]);

  printf("\ndiag : ");
  for (i=0;i<B->size;i++)
    printf("%lf  ",B->T->diag[i]);
 
  printf("\nupdiag : ");
  for (i=0;i<B->size-1;i++)
    printf("%lf  ",B->T->updiag[i]);


  printf("\nb : ");
  for (i=0;i<B->size;i++)
    printf("%lf  ",B->b[i]);

}



