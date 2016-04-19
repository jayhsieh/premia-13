/* zggqrf.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "pnl/pnl_f2c.h"

/* Table of constant values */

static int c__1 = 1;
static int c_n1 = -1;

 int zggqrf_(int *n, int *m, int *p, 
	doublecomplex *a, int *lda, doublecomplex *taua, doublecomplex *b, 
	 int *ldb, doublecomplex *taub, doublecomplex *work, int *
	lwork, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    int nb, nb1, nb2, nb3, lopt;
    extern  int xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *);
    extern  int zgeqrf_(int *, int *, doublecomplex *, 
	     int *, doublecomplex *, doublecomplex *, int *, int *
), zgerqf_(int *, int *, doublecomplex *, int *, 
	    doublecomplex *, doublecomplex *, int *, int *);
    int lwkopt;
    int lquery;
    extern  int zunmqr_(char *, char *, int *, int *, 
	    int *, doublecomplex *, int *, doublecomplex *, 
	    doublecomplex *, int *, doublecomplex *, int *, int *);


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  ZGGQRF computes a generalized QR factorization of an N-by-M matrix A */
/*  and an N-by-P matrix B: */

/*              A = Q*R,        B = Q*T*Z, */

/*  where Q is an N-by-N unitary matrix, Z is a P-by-P unitary matrix, */
/*  and R and T assume one of the forms: */

/*  if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N, */
/*                  (  0  ) N-M                         N   M-N */
/*                     M */

/*  where R11 is upper triangular, and */

/*  if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P, */
/*                   P-N  N                           ( T21 ) P */
/*                                                       P */

/*  where T12 or T21 is upper triangular. */

/*  In particular, if B is square and nonsingular, the GQR factorization */
/*  of A and B implicitly gives the QR factorization of inv(B)*A: */

/*               inv(B)*A = Z'*(inv(T)*R) */

/*  where inv(B) denotes the inverse of the matrix B, and Z' denotes the */
/*  conjugate transpose of matrix Z. */

/*  Arguments */
/*  ========= */

/*  N       (input) INTEGER */
/*          The number of rows of the matrices A and B. N >= 0. */

/*  M       (input) INTEGER */
/*          The number of columns of the matrix A.  M >= 0. */

/*  P       (input) INTEGER */
/*          The number of columns of the matrix B.  P >= 0. */

/*  A       (input/output) COMPLEX*16 array, dimension (LDA,M) */
/*          On entry, the N-by-M matrix A. */
/*          On exit, the elements on and above the diagonal of the array */
/*          contain the MIN(N,M)-by-M upper trapezoidal matrix R (R is */
/*          upper triangular if N >= M); the elements below the diagonal, */
/*          with the array TAUA, represent the unitary matrix Q as a */
/*          product of MIN(N,M) elementary reflectors (see Further */
/*          Details). */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. LDA >= MAX(1,N). */

/*  TAUA    (output) COMPLEX*16 array, dimension (MIN(N,M)) */
/*          The scalar factors of the elementary reflectors which */
/*          represent the unitary matrix Q (see Further Details). */

/*  B       (input/output) COMPLEX*16 array, dimension (LDB,P) */
/*          On entry, the N-by-P matrix B. */
/*          On exit, if N <= P, the upper triangle of the subarray */
/*          B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T; */
/*          if N > P, the elements on and above the (N-P)-th subdiagonal */
/*          contain the N-by-P upper trapezoidal matrix T; the remaining */
/*          elements, with the array TAUB, represent the unitary */
/*          matrix Z as a product of elementary reflectors (see Further */
/*          Details). */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B. LDB >= MAX(1,N). */

/*  TAUB    (output) COMPLEX*16 array, dimension (MIN(N,P)) */
/*          The scalar factors of the elementary reflectors which */
/*          represent the unitary matrix Z (see Further Details). */

/*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. LWORK >= MAX(1,N,M,P). */
/*          For optimum performance LWORK >= MAX(N,M,P)*MAX(NB1,NB2,NB3), */
/*          where NB1 is the optimal blocksize for the QR factorization */
/*          of an N-by-M matrix, NB2 is the optimal blocksize for the */
/*          RQ factorization of an N-by-P matrix, and NB3 is the optimal */
/*          blocksize for a call of ZUNMQR. */

/*          If LWORK = -1, then a workspace query is assumed; the routine */
/*          only calculates the optimal size of the WORK array, returns */
/*          this value as the first entry of the WORK array, and no error */
/*          message related to LWORK is issued by XERBLA. */

/*  INFO    (output) INTEGER */
/*           = 0:  successful exit */
/*           < 0:  if INFO = -i, the i-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(1) H(2) . . . H(k), where k = MIN(n,m). */

/*  Each H(i) has the form */

/*     H(i) = I - taua * v * v' */

/*  where taua is a complex scalar, and v is a complex vector with */
/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i+1:n,i), */
/*  and taua in TAUA(i). */
/*  To form Q explicitly, use LAPACK subroutine ZUNGQR. */
/*  To use Q to update another matrix, use LAPACK subroutine ZUNMQR. */

/*  The matrix Z is represented as a product of elementary reflectors */

/*     Z = H(1) H(2) . . . H(k), where k = MIN(n,p). */

/*  Each H(i) has the form */

/*     H(i) = I - taub * v * v' */

/*  where taub is a complex scalar, and v is a complex vector with */
/*  v(p-k+i+1:p) = 0 and v(p-k+i) = 1; v(1:p-k+i-1) is stored on exit in */
/*  B(n-k+i,1:p-k+i-1), and taub in TAUB(i). */
/*  To form Z explicitly, use LAPACK subroutine ZUNGRQ. */
/*  To use Z to update another matrix, use LAPACK subroutine ZUNMRQ. */

/*  ===================================================================== */

/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --taua;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    --taub;
    --work;

    /* Function Body */
    *info = 0;
    nb1 = ilaenv_(&c__1, "ZGEQRF", " ", n, m, &c_n1, &c_n1);
    nb2 = ilaenv_(&c__1, "ZGERQF", " ", n, p, &c_n1, &c_n1);
    nb3 = ilaenv_(&c__1, "ZUNMQR", " ", n, m, p, &c_n1);
/* Computing MAX */
    i__1 = MAX(nb1,nb2);
    nb = MAX(i__1,nb3);
/* Computing MAX */
    i__1 = MAX(*n,*m);
    lwkopt = MAX(i__1,*p) * nb;
    work[1].r = (double) lwkopt, work[1].i = 0.;
    lquery = *lwork == -1;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = MAX(1,*n), i__1 = MAX(i__1,*m);
	if (*lwork < MAX(i__1,*p) && ! lquery) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("ZGGQRF", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     QR factorization of N-by-M matrix A: A = Q*R */

    zgeqrf_(n, m, &a[a_offset], lda, &taua[1], &work[1], lwork, info);
    lopt = (int) work[1].r;

/*     Update B := Q'*B. */

    i__1 = MIN(*n,*m);
    zunmqr_("Left", "Conjugate Transpose", n, p, &i__1, &a[a_offset], lda, &
	    taua[1], &b[b_offset], ldb, &work[1], lwork, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) work[1].r;
    lopt = MAX(i__1,i__2);

/*     RQ factorization of N-by-P matrix B: B = T*Z. */

    zgerqf_(n, p, &b[b_offset], ldb, &taub[1], &work[1], lwork, info);
/* Computing MAX */
    i__2 = lopt, i__3 = (int) work[1].r;
    i__1 = MAX(i__2,i__3);
    work[1].r = (double) i__1, work[1].i = 0.;

    return 0;

/*     End of ZGGQRF */

} /* zggqrf_ */
