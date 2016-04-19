//
// MATHFI Project, Inria Rocquencourt.
// Vincent Barette, June 2002.
//
 
/* 
** FILE
** *****************************************************************
** 
** cholesky.h
** 
** PURPOSE
** 
** Some tools for the cholesky factorization LDL' of a matrix.
** 
*****************************************************************
*/

#ifndef CHOLESKY
#define CHOLESKY


struct chol 
{
  int size;
  double* l; //the length of l is n(n-1)/2 where n=size
  double* d; //the length of d is n where n=size
};

// D=diag( d0, d1, ...,d(n-1) )
//
// L= |1                                |
//    |l0      1                        |
//    |l1   l(n-1)   1                  |
//    |l2      .         .              |
//    |.       .              .         |
//    |.       .                  1     |
//    |l(n-2)  .                 l(K)  1|  where K=n(n-1)/2-1



void init(struct chol *C);//initialize with the identity matrix.


double elementLij(int i,int j,struct chol *C); //Returns the element of the vector l correspoding to the element (i,j) of the matrix L. i>j !!!


double elementDi(int i,struct chol *C);//Returns the ith element of the vector d.


double elementMij(int i,int j,struct chol *C);//returns the element (i,j) of the product LDL' (for i>=j, remember Mij=Mji)


void modify_l(int i,int j,double ll,struct chol *C);//Changes the element (i,j) of the matrix L into the value ll


void modify_d(int i,double dd,struct chol *C);//Changes the ith element of the matrix D into the value dd


void scalarMult(struct chol *C, double a);//Multiplies the matrix by the scalar a


void prodMatVect(struct chol *C,double* v, double* res);//Returns in res the product Cv


void Lsolve(struct chol *C, double* x, double* b);//Changes x into the solution of the linear system Lx=b


void LTsolve(struct chol *C, double* x, double* b);//Changes x into the solution of the linear system L'x=b


void cholSolve(struct chol *C, double* x, double* b);//Changes x into the solution of the linear system Ax=b


void cholSolveOpp(struct chol *C, double* x, double* b);//Changes x into the solution of the linear system Ax=-b


void positiveRankOneModification(struct chol *C, double s, double* z);//C is changed into the cholesky factorization of C + s*z'z (s>0)


void negativeRankOneModification(struct chol *C, double s, double* z, double* v, double epsilon);// C is changed into the cholesky factorization of C - s*z'z (s>0). The algorithm is different from the previous one: 
//- a vector v s.t. Lv=z is needed.
//- the machine precision "epsilon" is needed.


void printChol(struct chol *C);//prints the matrix L and D in the terminal
void printCholM(struct chol *C);//prints the matrix M=LDL' in the terminal.


#endif
