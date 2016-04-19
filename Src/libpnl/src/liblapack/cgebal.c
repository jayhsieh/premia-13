/* cgebal.f -- translated by f2c (version 20061008).
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

 int cgebal_(char *job, int *n, complex *a, int *lda, 
	int *ilo, int *ihi, float *scale, int *info)
{
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;
    float r__1, r__2;

    /* Builtin functions */
    double r_imag(complex *), c_abs(complex *);

    /* Local variables */
    float c__, f, g;
    int i__, j, k, l, m;
    float r__, s, ca, ra;
    int ica, ira, iexc;
    extern int lsame_(char *, char *);
    extern  int cswap_(int *, complex *, int *, 
	    complex *, int *);
    float sfmin1, sfmin2, sfmax1, sfmax2;
    extern int icamax_(int *, complex *, int *);
    extern double slamch_(char *);
    extern  int csscal_(int *, float *, complex *, int 
	    *), xerbla_(char *, int *);
    int noconv;


/*  -- LAPACK routine (version 3.2) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CGEBAL balances a general complex matrix A.  This involves, first, */
/*  permuting A by a similarity transformation to isolate eigenvalues */
/*  in the first 1 to ILO-1 and last IHI+1 to N elements on the */
/*  diagonal; and second, applying a diagonal similarity transformation */
/*  to rows and columns ILO to IHI to make the rows and columns as */
/*  close in norm as possible.  Both steps are optional. */

/*  Balancing may reduce the 1-norm of the matrix, and improve the */
/*  accuracy of the computed eigenvalues and/or eigenvectors. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) CHARACTER*1 */
/*          Specifies the operations to be performed on A: */
/*          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0 */
/*                  for i = 1,...,N; */
/*          = 'P':  permute only; */
/*          = 'S':  scale only; */
/*          = 'B':  both permute and scale. */

/*  N       (input) INTEGER */
/*          The order of the matrix A.  N >= 0. */

/*  A       (input/output) COMPLEX array, dimension (LDA,N) */
/*          On entry, the input matrix A. */
/*          On exit,  A is overwritten by the balanced matrix. */
/*          If JOB = 'N', A is not referenced. */
/*          See Further Details. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A.  LDA >= MAX(1,N). */

/*  ILO     (output) INTEGER */
/*  IHI     (output) INTEGER */
/*          ILO and IHI are set to ints such that on exit */
/*          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N. */
/*          If JOB = 'N' or 'S', ILO = 1 and IHI = N. */

/*  SCALE   (output) REAL array, dimension (N) */
/*          Details of the permutations and scaling factors applied to */
/*          A.  If P(j) is the index of the row and column interchanged */
/*          with row and column j and D(j) is the scaling factor */
/*          applied to row and column j, then */
/*          SCALE(j) = P(j)    for j = 1,...,ILO-1 */
/*                   = D(j)    for j = ILO,...,IHI */
/*                   = P(j)    for j = IHI+1,...,N. */
/*          The order in which the interchanges are made is N to IHI+1, */
/*          then 1 to ILO-1. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit. */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*  The permutations consist of row and column interchanges which put */
/*  the matrix in the form */

/*             ( T1   X   Y  ) */
/*     P A P = (  0   B   Z  ) */
/*             (  0   0   T2 ) */

/*  where T1 and T2 are upper triangular matrices whose eigenvalues lie */
/*  along the diagonal.  The column indices ILO and IHI mark the starting */
/*  and ending columns of the submatrix B. Balancing consists of applying */
/*  a diagonal similarity transformation inv(D) * B * D to make the */
/*  1-norms of each row of B and its corresponding column nearly equal. */
/*  The output matrix is */

/*     ( T1     X*D          Y    ) */
/*     (  0  inv(D)*B*D  inv(D)*Z ). */
/*     (  0      0           T2   ) */

/*  Information about the permutations P and the diagonal matrix D is */
/*  returned in the vector SCALE. */

/*  This subroutine is based on the EISPACK routine CBAL. */

/*  Modified by Tzu-Yi Chen, Computer Science Division, University of */
/*    California at Berkeley, USA */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --scale;

    /* Function Body */
    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(job, "S") 
	    && ! lsame_(job, "B")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CGEBAL", &i__1);
	return 0;
    }

    k = 1;
    l = *n;

    if (*n == 0) {
	goto L210;
    }

    if (lsame_(job, "N")) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    scale[i__] = 1.f;
/* L10: */
	}
	goto L210;
    }

    if (lsame_(job, "S")) {
	goto L120;
    }

/*     Permutation to isolate eigenvalues if possible */

    goto L50;

/*     Row and column exchange. */

L20:
    scale[m] = (float) j;
    if (j == m) {
	goto L30;
    }

    cswap_(&l, &a[j * a_dim1 + 1], &c__1, &a[m * a_dim1 + 1], &c__1);
    i__1 = *n - k + 1;
    cswap_(&i__1, &a[j + k * a_dim1], lda, &a[m + k * a_dim1], lda);

L30:
    switch (iexc) {
	case 1:  goto L40;
	case 2:  goto L80;
    }

/*     Search for rows isolating an eigenvalue and push them down. */

L40:
    if (l == 1) {
	goto L210;
    }
    --l;

L50:
    for (j = l; j >= 1; --j) {

	i__1 = l;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (i__ == j) {
		goto L60;
	    }
	    i__2 = j + i__ * a_dim1;
	    if (a[i__2].r != 0.f || r_imag(&a[j + i__ * a_dim1]) != 0.f) {
		goto L70;
	    }
L60:
	    ;
	}

	m = l;
	iexc = 1;
	goto L20;
L70:
	;
    }

    goto L90;

/*     Search for columns isolating an eigenvalue and push them left. */

L80:
    ++k;

L90:
    i__1 = l;
    for (j = k; j <= i__1; ++j) {

	i__2 = l;
	for (i__ = k; i__ <= i__2; ++i__) {
	    if (i__ == j) {
		goto L100;
	    }
	    i__3 = i__ + j * a_dim1;
	    if (a[i__3].r != 0.f || r_imag(&a[i__ + j * a_dim1]) != 0.f) {
		goto L110;
	    }
L100:
	    ;
	}

	m = k;
	iexc = 2;
	goto L20;
L110:
	;
    }

L120:
    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) {
	scale[i__] = 1.f;
/* L130: */
    }

    if (lsame_(job, "P")) {
	goto L210;
    }

/*     Balance the submatrix in rows K to L. */

/*     Iterative loop for norm reduction */

    sfmin1 = slamch_("S") / slamch_("P");
    sfmax1 = 1.f / sfmin1;
    sfmin2 = sfmin1 * 2.f;
    sfmax2 = 1.f / sfmin2;
L140:
    noconv = FALSE;

    i__1 = l;
    for (i__ = k; i__ <= i__1; ++i__) {
	c__ = 0.f;
	r__ = 0.f;

	i__2 = l;
	for (j = k; j <= i__2; ++j) {
	    if (j == i__) {
		goto L150;
	    }
	    i__3 = j + i__ * a_dim1;
	    c__ += (r__1 = a[i__3].r, ABS(r__1)) + (r__2 = r_imag(&a[j + i__ 
		    * a_dim1]), ABS(r__2));
	    i__3 = i__ + j * a_dim1;
	    r__ += (r__1 = a[i__3].r, ABS(r__1)) + (r__2 = r_imag(&a[i__ + j 
		    * a_dim1]), ABS(r__2));
L150:
	    ;
	}
	ica = icamax_(&l, &a[i__ * a_dim1 + 1], &c__1);
	ca = c_abs(&a[ica + i__ * a_dim1]);
	i__2 = *n - k + 1;
	ira = icamax_(&i__2, &a[i__ + k * a_dim1], lda);
	ra = c_abs(&a[i__ + (ira + k - 1) * a_dim1]);

/*        Guard against zero C or R due to underflow. */

	if (c__ == 0.f || r__ == 0.f) {
	    goto L200;
	}
	g = r__ / 2.f;
	f = 1.f;
	s = c__ + r__;
L160:
/* Computing MAX */
	r__1 = MAX(f,c__);
/* Computing MIN */
	r__2 = MIN(r__,g);
	if (c__ >= g || MAX(r__1,ca) >= sfmax2 || MIN(r__2,ra) <= sfmin2) {
	    goto L170;
	}
	f *= 2.f;
	c__ *= 2.f;
	ca *= 2.f;
	r__ /= 2.f;
	g /= 2.f;
	ra /= 2.f;
	goto L160;

L170:
	g = c__ / 2.f;
L180:
/* Computing MIN */
	r__1 = MIN(f,c__), r__1 = MIN(r__1,g);
	if (g < r__ || MAX(r__,ra) >= sfmax2 || MIN(r__1,ca) <= sfmin2) {
	    goto L190;
	}
	f /= 2.f;
	c__ /= 2.f;
	g /= 2.f;
	ca /= 2.f;
	r__ *= 2.f;
	ra *= 2.f;
	goto L180;

/*        Now balance. */

L190:
	if (c__ + r__ >= s * .95f) {
	    goto L200;
	}
	if (f < 1.f && scale[i__] < 1.f) {
	    if (f * scale[i__] <= sfmin1) {
		goto L200;
	    }
	}
	if (f > 1.f && scale[i__] > 1.f) {
	    if (scale[i__] >= sfmax1 / f) {
		goto L200;
	    }
	}
	g = 1.f / f;
	scale[i__] *= f;
	noconv = TRUE;

	i__2 = *n - k + 1;
	csscal_(&i__2, &g, &a[i__ + k * a_dim1], lda);
	csscal_(&l, &f, &a[i__ * a_dim1 + 1], &c__1);

L200:
	;
    }

    if (noconv) {
	goto L140;
    }

L210:
    *ilo = k;
    *ihi = l;

    return 0;

/*     End of CGEBAL */

} /* cgebal_ */
