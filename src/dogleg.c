/* dogleg.f -- translated by f2c (version 12.02.01).
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

/* DECK DOGLEG */
/* Subroutine */ int dogleg_(integer *n, real *r__, integer *lr, real *diag, 
	real *qtb, real *delta, real *x, real *wa1, real *wa2)
{
    /* Initialized data */

    static real one = 1.f;
    static real zero = 0.f;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2, r__3, r__4;

    /* Local variables */
    static integer i__, j, k, l, jj, jp1;
    static real sum, temp, alpha, bnorm;
    extern doublereal enorm_(integer *, real *);
    static real gnorm, qnorm;
    extern doublereal r1mach_(integer *);
    static real epsmch, sgnorm;

/* ***BEGIN PROLOGUE  DOGLEG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SNSQ and SNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (DOGLEG-S, DDOGLG-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an M by N matrix A, an N by N nonsingular DIAGONAL */
/*     matrix D, an M-vector B, and a positive number DELTA, the */
/*     problem is to determine the convex combination X of the */
/*     Gauss-Newton and scaled gradient directions that minimizes */
/*     (A*X - B) in the least squares sense, subject to the */
/*     restriction that the Euclidean norm of D*X be at most DELTA. */

/*     This subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     QR factorization of A. That is, if A = Q*R, where Q has */
/*     orthogonal columns and R is an upper triangular matrix, */
/*     then DOGLEG expects the full upper triangle of R and */
/*     the first N components of (Q TRANSPOSE)*B. */

/*     The subroutine statement is */

/*       SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an input array of length LR which must contain the upper */
/*         triangular matrix R stored by rows. */

/*       LR is a positive integer input variable not less than */
/*         (N*(N+1))/2. */

/*       DIAG is an input array of length N which must contain the */
/*         diagonal elements of the matrix D. */

/*       QTB is an input array of length N which must contain the first */
/*         N elements of the vector (Q TRANSPOSE)*B. */

/*       DELTA is a positive input variable which specifies an upper */
/*         bound on the Euclidean norm of D*X. */

/*       X is an output array of length N which contains the desired */
/*         convex combination of the Gauss-Newton direction and the */
/*         scaled gradient direction. */

/*       WA1 and WA2 are work arrays of length N. */

/* ***SEE ALSO  SNSQ, SNSQE */
/* ***ROUTINES CALLED  ENORM, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DOGLEG */
    /* Parameter adjustments */
    --r__;
    --diag;
    --qtb;
    --x;
    --wa1;
    --wa2;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DOGLEG */
    epsmch = r1mach_(&c__4);

/*     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION. */

    jj = *n * (*n + 1) / 2 + 1;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	j = *n - k + 1;
	jp1 = j + 1;
	jj -= k;
	l = jj + 1;
	sum = zero;
	if (*n < jp1) {
	    goto L20;
	}
	i__2 = *n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	    sum += r__[l] * x[i__];
	    ++l;
/* L10: */
	}
L20:
	temp = r__[jj];
	if (temp != zero) {
	    goto L40;
	}
	l = j;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
	    r__2 = temp, r__3 = (r__1 = r__[l], dabs(r__1));
	    temp = dmax(r__2,r__3);
	    l = l + *n - i__;
/* L30: */
	}
	temp = epsmch * temp;
	if (temp == zero) {
	    temp = epsmch;
	}
L40:
	x[j] = (qtb[j] - sum) / temp;
/* L50: */
    }

/*     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = zero;
	wa2[j] = diag[j] * x[j];
/* L60: */
    }
    qnorm = enorm_(n, &wa2[1]);
    if (qnorm <= *delta) {
	goto L140;
    }

/*     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE. */
/*     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION. */

    l = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	temp = qtb[j];
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    wa1[i__] += r__[l] * temp;
	    ++l;
/* L70: */
	}
	wa1[j] /= diag[j];
/* L80: */
    }

/*     CALCULATE THE NORM OF THE SCALED GRADIENT DIRECTION, */
/*     NORMALIZE, AND RESCALE THE GRADIENT. */

    gnorm = enorm_(n, &wa1[1]);
    sgnorm = zero;
    alpha = *delta / qnorm;
    if (gnorm == zero) {
	goto L120;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = wa1[j] / gnorm / diag[j];
/* L90: */
    }

/*     CALCULATE THE POINT ALONG THE SCALED GRADIENT */
/*     AT WHICH THE QUADRATIC IS MINIMIZED. */

    l = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum = zero;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum += r__[l] * wa1[i__];
	    ++l;
/* L100: */
	}
	wa2[j] = sum;
/* L110: */
    }
    temp = enorm_(n, &wa2[1]);
    sgnorm = gnorm / temp / temp;

/*     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE. */

    alpha = zero;
    if (sgnorm >= *delta) {
	goto L120;
    }

/*     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE. */
/*     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG */
/*     AT WHICH THE QUADRATIC IS MINIMIZED. */

    bnorm = enorm_(n, &qtb[1]);
    temp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / *delta);
/* Computing 2nd power */
    r__1 = sgnorm / *delta;
/* Computing 2nd power */
    r__2 = temp - *delta / qnorm;
/* Computing 2nd power */
    r__3 = *delta / qnorm;
/* Computing 2nd power */
    r__4 = sgnorm / *delta;
    temp = temp - *delta / qnorm * (r__1 * r__1) + sqrt(r__2 * r__2 + (one - 
	    r__3 * r__3) * (one - r__4 * r__4));
/* Computing 2nd power */
    r__1 = sgnorm / *delta;
    alpha = *delta / qnorm * (one - r__1 * r__1) / temp;
L120:

/*     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON */
/*     DIRECTION AND THE SCALED GRADIENT DIRECTION. */

    temp = (one - alpha) * dmin(sgnorm,*delta);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = temp * wa1[j] + alpha * x[j];
/* L130: */
    }
L140:
    return 0;

/*     LAST CARD OF SUBROUTINE DOGLEG. */

} /* dogleg_ */

