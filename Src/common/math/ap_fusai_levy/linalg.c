/*===================================================================
 linalg.c

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/


#include <math.h>

#include "nrutil.h"


//double determinant(double **a, int n)


// inverse(double **a, int n)

//NUMERICS_EXPORT BOOL linsolve(double **m, double *b, int n, int method)


// void lubksb(double **m, int n, int *indx, double *b)


// ludcmp(double **m, int n, int *indx, double *d)


void matmat(double **a, int nra, int nca, double **b, int ncb, double **prod)
{
    int i, j, k;
    double sum;
    
    for(i = 0; i < nra; i++){
        for(j = 0; j < ncb; j++){
            sum = 0.0;
            for(k = 0; k < nca; k++){
                sum += a[i][k] * b[k][j];
            }
            prod[i][j] = sum;
        }
    }
}

void matvec(double **a, int nra, int nca, double *x, double *b)
{
    int i, j;
    double sum;
    
    for(i = 0; i < nra; i++){
        sum = 0.0;
        for(j = 0; j < nca; j++){
            sum += a[i][j] * x[j];
		//	sum += a[i][j] * x[j];
        }
        b[i] = sum;
    }
}

void transpose(double **a, int nr, int nc, double **at)
{
    int i, j;
    
    for(i = 0; i < nr; ++i){
        for(j = 0; j < nc; ++j){
            at[j][i] = a[i][j];
        }
    }
}

void vecmat(double *x, double **a, int nra, int nca, double *b)
{
    double** t = dmatrix(0, nca - 1, 0, nra - 1);
    
    transpose(a, nra, nca, t);
    matvec(t, nca, nra, x, b);
    
    free_dmatrix(t, 0, nca - 1, 0,nra - 1);
}

double vecvec(double *first1, double* last1, double* first2)
{
    double p = 1.0;
    
    while(first1 < last1){
        p += *first1 * *first2;
        ++first1;
        ++first2;
    }
    
    return p;
}

void pairwdiff(int n, double* x, double* y,   double **dest)
{
    int i;
    int j;

    for(i = 0; i < n; ++i){
        for(j = 0; j < n; ++j){
            dest[i][j] = x[i] - y[j];
        }
    }
}


/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/

