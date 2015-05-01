/* qrfac.f -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* Table of constant values */

static integer c__4 = 4;

/* DECK QRFAC */
/* Subroutine */ int qrfac_(integer *m, integer *n, real *a, integer *lda, 
	logical *pivot, integer *ipvt, integer *lipvt, real *sigma, real *
	acnorm, real *wa)
{
    /* Initialized data */

    static real one = 1.f;
    static real p05 = .05f;
    static real zero = 0.f;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, j, k, jp1;
    static real sum;
    static integer kmax;
    static real temp;
    static integer minmn;
    extern doublereal enorm_(integer *, real *), r1mach_(integer *);
    static real epsmch, ajnorm;

/* ***BEGIN PROLOGUE  QRFAC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SNLS1, SNLS1E, SNSQ and SNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (QRFAC-S, DQRFAC-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine uses Householder transformations with column */
/*     pivoting (optional) to compute a QR factorization of the */
/*     M by N matrix A. That is, QRFAC determines an orthogonal */
/*     matrix Q, a permutation matrix P, and an upper trapezoidal */
/*     matrix R with diagonal elements of nonincreasing magnitude, */
/*     such that A*P = Q*R. The Householder transformation for */
/*     column K, K = 1,2,...,MIN(M,N), is of the form */

/*                           T */
/*           I - (1/U(K))*U*U */

/*     where U has zeros in the first K-1 positions. The form of */
/*     this transformation and the method of pivoting first */
/*     appeared in the corresponding LINPACK subroutine. */

/*     The subroutine statement is */

/*       SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA) */

/*     where */

/*       M is a positive integer input variable set to the number */
/*         of rows of A. */

/*       N is a positive integer input variable set to the number */
/*         of columns of A. */

/*       A is an M by N array. On input A contains the matrix for */
/*         which the QR factorization is to be computed. On output */
/*         the strict upper trapezoidal part of A contains the strict */
/*         upper trapezoidal part of R, and the lower trapezoidal */
/*         part of A contains a factored form of Q (the non-trivial */
/*         elements of the U vectors described above). */

/*       LDA is a positive integer input variable not less than M */
/*         which specifies the leading dimension of the array A. */

/*       PIVOT is a logical input variable. If pivot is set .TRUE., */
/*         then column pivoting is enforced. If pivot is set .FALSE., */
/*         then no column pivoting is done. */

/*       IPVT is an integer output array of length LIPVT. IPVT */
/*         defines the permutation matrix P such that A*P = Q*R. */
/*         Column J of P is column IPVT(J) of the identity matrix. */
/*         If pivot is .FALSE., IPVT is not referenced. */

/*       LIPVT is a positive integer input variable. If PIVOT is */
/*             .FALSE., then LIPVT may be as small as 1. If PIVOT is */
/*             .TRUE., then LIPVT must be at least N. */

/*       SIGMA is an output array of length N which contains the */
/*         diagonal elements of R. */

/*       ACNORM is an output array of length N which contains the */
/*         norms of the corresponding columns of the input matrix A. */
/*         If this information is not needed, then ACNORM can coincide */
/*         with SIGMA. */

/*       WA is a work array of length N. If pivot is .FALSE., then WA */
/*         can coincide with SIGMA. */

/* ***SEE ALSO  SNLS1, SNLS1E, SNSQ, SNSQE */
/* ***ROUTINES CALLED  ENORM, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  QRFAC */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --sigma;
    --acnorm;
    --wa;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  QRFAC */
    epsmch = r1mach_(&c__4);

/*     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	acnorm[j] = enorm_(m, &a[j * a_dim1 + 1]);
	sigma[j] = acnorm[j];
	wa[j] = sigma[j];
	if (*pivot) {
	    ipvt[j] = j;
	}
/* L10: */
    }

/*     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS. */

    minmn = min(*m,*n);
    i__1 = minmn;
    for (j = 1; j <= i__1; ++j) {
	if (! (*pivot)) {
	    goto L40;
	}

/*        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION. */

	kmax = j;
	i__2 = *n;
	for (k = j; k <= i__2; ++k) {
	    if (sigma[k] > sigma[kmax]) {
		kmax = k;
	    }
/* L20: */
	}
	if (kmax == j) {
	    goto L40;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[i__ + kmax * a_dim1];
	    a[i__ + kmax * a_dim1] = temp;
/* L30: */
	}
	sigma[kmax] = sigma[j];
	wa[kmax] = wa[j];
	k = ipvt[j];
	ipvt[j] = ipvt[kmax];
	ipvt[kmax] = k;
L40:

/*        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE */
/*        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR. */

	i__2 = *m - j + 1;
	ajnorm = enorm_(&i__2, &a[j + j * a_dim1]);
	if (ajnorm == zero) {
	    goto L100;
	}
	if (a[j + j * a_dim1] < zero) {
	    ajnorm = -ajnorm;
	}
	i__2 = *m;
	for (i__ = j; i__ <= i__2; ++i__) {
	    a[i__ + j * a_dim1] /= ajnorm;
/* L50: */
	}
	a[j + j * a_dim1] += one;

/*        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS */
/*        AND UPDATE THE NORMS. */

	jp1 = j + 1;
	if (*n < jp1) {
	    goto L100;
	}
	i__2 = *n;
	for (k = jp1; k <= i__2; ++k) {
	    sum = zero;
	    i__3 = *m;
	    for (i__ = j; i__ <= i__3; ++i__) {
		sum += a[i__ + j * a_dim1] * a[i__ + k * a_dim1];
/* L60: */
	    }
	    temp = sum / a[j + j * a_dim1];
	    i__3 = *m;
	    for (i__ = j; i__ <= i__3; ++i__) {
		a[i__ + k * a_dim1] -= temp * a[i__ + j * a_dim1];
/* L70: */
	    }
	    if (! (*pivot) || sigma[k] == zero) {
		goto L80;
	    }
	    temp = a[j + k * a_dim1] / sigma[k];
/* Computing MAX */
/* Computing 2nd power */
	    r__3 = temp;
	    r__1 = zero, r__2 = one - r__3 * r__3;
	    sigma[k] *= sqrt((dmax(r__1,r__2)));
/* Computing 2nd power */
	    r__1 = sigma[k] / wa[k];
	    if (p05 * (r__1 * r__1) > epsmch) {
		goto L80;
	    }
	    i__3 = *m - j;
	    sigma[k] = enorm_(&i__3, &a[jp1 + k * a_dim1]);
	    wa[k] = sigma[k];
L80:
/* L90: */
	    ;
	}
L100:
	sigma[j] = -ajnorm;
/* L110: */
    }
    return 0;

/*     LAST CARD OF SUBROUTINE QRFAC. */

} /* qrfac_ */

