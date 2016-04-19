/* sgeqr2.f -- translated by f2c (version 20061008).
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

 int sgeqr2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    int i__, k;
    float aii;
    extern  int slarf_(char *, int *, int *, float *, 
	    int *, float *, float *, int *, float *), xerbla_(
	    char *, int *), slarfp_(int *, float *, float *, 
	    int *, float *);


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGEQR2 computes a QR factorization of a float m by n matrix A: */
/*  A = Q * R. */

/*  Arguments */
/*  ========= */

/*  M       (input) INTEGER */
/*          The number of rows of the matrix A.  M >= 0. */

/*  N       (input) INTEGER */
/*          The number of columns of the matrix A.  N >= 0. */

/*  A       (input/output) REAL array, dimension (LDA,N) */
/*          On entry, the m by n matrix A. */
/*          On exit, the elements on and above the diagonal of the array */
/*          contain the MIN(m,n) by n upper trapezoidal matrix R (R is */
/*          upper triangular if m >= n); the elements below the diagonal, */
/*          with the array TAU, represent the orthogonal matrix Q as a */
/*          product of elementary reflectors (see Further Details). */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= MAX(1,M). */

/*  TAU     (output) REAL array, dimension (MIN(M,N)) */
/*          The scalar factors of the elementary reflectors (see Further */
/*          Details). */

/*  WORK    (workspace) REAL array, dimension (N) */

/*  INFO    (output) INTEGER */
/*          = 0: successful exit */
/*          < 0: if INFO = -i, the i-th argument had an illegal value */

/*  Further Details */
/*  =============== */

/*  The matrix Q is represented as a product of elementary reflectors */

/*     Q = H(1) H(2) . . . H(k), where k = MIN(m,n). */

/*  Each H(i) has the form */

/*     H(i) = I - tau * v * v' */

/*  where tau is a float scalar, and v is a float vector with */
/*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i), */
/*  and tau in TAU(i). */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("SGEQR2", &i__1);
	return 0;
    }

    k = MIN(*m,*n);

    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i) */

	i__2 = *m - i__ + 1;
/* Computing MIN */
	i__3 = i__ + 1;
	slarfp_(&i__2, &a[i__ + i__ * a_dim1], &a[MIN(i__3, *m)+ i__ * a_dim1]
, &c__1, &tau[i__]);
	if (i__ < *n) {

/*           Apply H(i) to A(i:m,i+1:n) from the left */

	    aii = a[i__ + i__ * a_dim1];
	    a[i__ + i__ * a_dim1] = 1.f;
	    i__2 = *m - i__ + 1;
	    i__3 = *n - i__;
	    slarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &tau[
		    i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1]);
	    a[i__ + i__ * a_dim1] = aii;
	}
/* L10: */
    }
    return 0;

/*     End of SGEQR2 */

} /* sgeqr2_ */