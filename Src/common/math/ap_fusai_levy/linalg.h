
void matmat(double **a, int nra, int nca, double **b, int ncb, double **prod);
/*-------------------------------------------------------------------
 Postmultiplies the matrix a[0..nra-1][0..nca-1] by the matrix
 b[0..nca-1][0..ncb-1] and returns the product in the matrix
 prod[0..nra-1][0..ncb-1].
-------------------------------------------------------------------*/

void matvec(double **a, int nra, int nca, double *x, double *b);
/*-------------------------------------------------------------------
 Postmultiplies the matrix a[0..nra-1][0..nca-1] by the vector
 x[0..nca-1] and returns the product in the vector b[0..nra-1].
-------------------------------------------------------------------*/

void transpose(double **a, int nr, int nc, double **at);
/*-------------------------------------------------------------------
 Returns the transpose of a[0..nr-1][0..nc-1] as
 at[0..nc-1][0..nr-1].
-------------------------------------------------------------------*/

void vecmat(double *x, double **a, int nra, int nca, double *b);
/*-------------------------------------------------------------------
 Premultiplies the matrix a[0..nra-1][0..nca-1] by the vector
 x[0..nra-1] and returns the product in the vector b[0..nca-1].
-------------------------------------------------------------------*/

double vecvec(double *first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns the inner product between the vectors u[0..n-1] and
 v[0..n-1].
-------------------------------------------------------------------*/

//void pairwdiff(double* first1, double* last1, double* first2, double* last2, double* dest);
void pairwdiff(int n, double* x, double* y,   double **dest);

/*-------------------------------------------------------------------
 Computes the pairwise differences between the elements in
 [first1, last1) and the elements in [first2, last2) and places them
 in dest.  dest must be large enough to hold all of the m * n
 differences, where m = last1 - first1 and n = last2 - first2.
-------------------------------------------------------------------*/

