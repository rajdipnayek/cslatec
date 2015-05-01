/* qform.f -- translated by f2c (version 12.02.01).
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

/* DECK QFORM */
/* Subroutine */ int qform_(integer *m, integer *n, real *q, integer *ldq, 
	real *wa)
{
    /* Initialized data */

    static real one = 1.f;
    static real zero = 0.f;

    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, jm1, np1;
    static real sum, temp;
    static integer minmn;

/* ***BEGIN PROLOGUE  QFORM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SNSQ and SNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (QFORM-S, DQFORM-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine proceeds from the computed QR factorization of */
/*     an M by N matrix A to accumulate the M by M orthogonal matrix */
/*     Q from its factored form. */

/*     The subroutine statement is */

/*       SUBROUTINE QFORM(M,N,Q,LDQ,WA) */

/*     where */

/*       M is a positive integer input variable set to the number */
/*         of rows of A and the order of Q. */

/*       N is a positive integer input variable set to the number */
/*         of columns of A. */

/*       Q is an M by M array. On input the full lower trapezoid in */
/*         the first min(M,N) columns of Q contains the factored form. */
/*         On output Q has been accumulated into a square matrix. */

/*       LDQ is a positive integer input variable not less than M */
/*         which specifies the leading dimension of the array Q. */

/*       WA is a work array of length M. */

/* ***SEE ALSO  SNSQ, SNSQE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  QFORM */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --wa;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  QFORM */
    minmn = min(*m,*n);
    if (minmn < 2) {
	goto L30;
    }
    i__1 = minmn;
    for (j = 2; j <= i__1; ++j) {
	jm1 = j - 1;
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + j * q_dim1] = zero;
/* L10: */
	}
/* L20: */
    }
L30:

/*     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX. */

    np1 = *n + 1;
    if (*m < np1) {
	goto L60;
    }
    i__1 = *m;
    for (j = np1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    q[i__ + j * q_dim1] = zero;
/* L40: */
	}
	q[j + j * q_dim1] = one;
/* L50: */
    }
L60:

/*     ACCUMULATE Q FROM ITS FACTORED FORM. */

    i__1 = minmn;
    for (l = 1; l <= i__1; ++l) {
	k = minmn - l + 1;
	i__2 = *m;
	for (i__ = k; i__ <= i__2; ++i__) {
	    wa[i__] = q[i__ + k * q_dim1];
	    q[i__ + k * q_dim1] = zero;
/* L70: */
	}
	q[k + k * q_dim1] = one;
	if (wa[k] == zero) {
	    goto L110;
	}
	i__2 = *m;
	for (j = k; j <= i__2; ++j) {
	    sum = zero;
	    i__3 = *m;
	    for (i__ = k; i__ <= i__3; ++i__) {
		sum += q[i__ + j * q_dim1] * wa[i__];
/* L80: */
	    }
	    temp = sum / wa[k];
	    i__3 = *m;
	    for (i__ = k; i__ <= i__3; ++i__) {
		q[i__ + j * q_dim1] -= temp * wa[i__];
/* L90: */
	    }
/* L100: */
	}
L110:
/* L120: */
	;
    }
    return 0;

/*     LAST CARD OF SUBROUTINE QFORM. */

} /* qform_ */

