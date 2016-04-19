#ifndef SYST
#define SYST 

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle, September 2002.
*/

/**************type definitions**********************/


/*tridiagonal matrices*/
struct tridiag{
  int size;
  double *subdiag; /*size-1*/
  double *diag;    /*size*/
  double *updiag;  /*size-1*/
}; 

/*bidiagonal matrices*/
struct bidiag{
  int size;
  double *subdiag; /*size-1*/
  double *diag;    /*size*/
}; 

/*tridiagonal system Tx=b */
struct tridiagSystem{
  int size;
  struct tridiag *T;
  double *b;
}; 

/*bidiagonal system Tx=b */
struct bidiagSystem{
  int size;
  struct bidiag *T;
  double *b;
};


/************************functions definitions****************************/

void tridiagTimesVect(double *res, struct tridiag *A, double *vect);
/* res = A*vect */
  
void tridiagToBidiagSyst(struct bidiagSystem *U, struct tridiagSystem *S);
/* AIM : changes the tridiagonal system S  to a triangular (bidiagonal) system U */

void tridiagToBidiagMat(struct bidiag *B, struct tridiag *T);
/* AIM : changes the tridiagonal matrix T to a bidiagonal matrix B */


void solveSyst(double *u_next, struct bidiagSystem *S);
/* AIM : solves the bidiagonal system S           */
/* PARAMETERS :                                   */
/* sol : solution of the system defined by S      */
/* S : bidiagonal system to be solved             */


void tridiagMatrixInv(double **inv, struct tridiag *T);
/* AIM: computes the inverse of the tridiag matrix represented by T */

void tridiagMatMult(double **res, struct tridiag *A, double **B);
/*multiplies the matrix B by the tridiagonal matrix A, both of size N */

void affectOperator(struct tridiag *A, struct tridiag *B);
/* AIM : equals the operator B to the operator A */

void printBidiagSyst(struct bidiagSystem *B);

void printTridiagSyst(struct tridiagSystem *B);


#endif
