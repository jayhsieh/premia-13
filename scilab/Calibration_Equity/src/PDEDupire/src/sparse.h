#ifndef _SPARSE
#define _SPARSE

#include "solveSystem.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

struct cell{
  int col;
  int line;
  double value;
  struct cell *right;
  struct cell *down;
};
 
struct sparseMat{
  int nbCol;
  int nbLines;
  struct cell **lines;
  struct cell **columns;
};


int isEmpty(struct sparseMat *A);

void transpSparseMat(struct sparseMat **pt_transp_A, struct sparseMat *A);

void sparseTimesSparse(struct sparseMat **pt_res, struct sparseMat *A, struct sparseMat *B);
/* (*pt_res) = A*B */

void sparseTimesVect(double *res, struct sparseMat *A, double *B);
/* res = A*B */

void sparseMinusSparse(struct sparseMat **res, struct sparseMat *A, struct sparseMat *B);
/* (*res) = A-B */

void sparsePlusSparse(struct sparseMat **res, struct sparseMat *A, struct sparseMat *B);
/* (*res) = A+B */

void affectSparse(struct sparseMat **pt_res, struct sparseMat *A);
/* (*pt_res) = A */

void regPlusSparse(double **sum, struct sparseMat *A, double **B);
/* computes sum = A+B */

void regTimesSparse(struct sparseMat **pt_res, double **A, struct sparseMat *B, int nbLinesA);
/* computes (*pt_res) = A*B */

void sparseMatMult(struct sparseMat **res, struct sparseMat *mat);
/* (*res) = transpose(mat)*mat */

void desallocSparseMat_cells(struct sparseMat *mat);
/* frees the cells of the sparse matrix "mat" */

void printSparseMat(struct sparseMat *mat);
/* prints the sparse matrix "mat" */

void testSparseMat(struct sparseMat **pt_mat);
/*builds a test sparse matrix (*pt_mat) */

void testSparseMat2(struct sparseMat **pt_mat);
/*builds a test sparse matrix (*pt_mat) */

void tridiagToSparse(struct tridiag *T, struct sparseMat **pt_res);
/* (*pt_res) = T */

void scalarTimesSparse(double a, struct sparseMat **M);
/* (*M) = a*(*M) */

#endif
