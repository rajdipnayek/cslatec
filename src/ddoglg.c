/* ddoglg.f -- translated by f2c (version 12.02.01).
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

/* DECK DDOGLG */
/* Subroutine */ int ddoglg_(integer *n, doublereal *r__, integer *lr, 
	doublereal *diag, doublereal *qtb, doublereal *delta, doublereal *x, 
	doublereal *wa1, doublereal *wa2)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i__, j, k, l, jj, jp1;
    static doublereal sum, temp, alpha, bnorm, gnorm, qnorm;
    extern doublereal d1mach_(integer *);
    static doublereal epsmch;
    extern doublereal denorm_(integer *, doublereal *);
    static doublereal sgnorm;

/* ***BEGIN PROLOGUE  DDOGLG */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DNSQ and DNSQE */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (DOGLEG-S, DDOGLG-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Given an M by N matrix A, an N by N nonsingular diagonal */
/*     matrix D, an M-vector B, and a positive number DELTA, the */
/*     problem is to determine the convex combination X of the */
/*     Gauss-Newton and scaled gradient directions that minimizes */
/*     (A*X - B) in the least squares sense, subject to the */
/*     restriction that the Euclidean norm of D*X be at most DELTA. */

/*     This subroutine completes the solution of the problem */
/*     if it is provided with the necessary information from the */
/*     QR factorization of A. That is, if A = Q*R, where Q has */
/*     orthogonal columns and R is an upper triangular matrix, */
/*     then DDOGLG expects the full upper triangle of R and */
/*     the first N components of (Q transpose)*B. */

/*     The subroutine statement is */

/*       SUBROUTINE DDOGLG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2) */

/*     where */

/*       N is a positive integer input variable set to the order of R. */

/*       R is an input array of length LR which must contain the upper */
/*         triangular matrix R stored by rows. */

/*       LR is a positive integer input variable not less than */
/*         (N*(N+1))/2. */

/*       DIAG is an input array of length N which must contain the */
/*         diagonal elements of the matrix D. */

/*       QTB is an input array of length N which must contain the first */
/*         N elements of the vector (Q transpose)*B. */

/*       DELTA is a positive input variable which specifies an upper */
/*         bound on the Euclidean norm of D*X. */

/*       X is an output array of length N which contains the desired */
/*         convex combination of the Gauss-Newton direction and the */
/*         scaled gradient direction. */

/*       WA1 and WA2 are work arrays of length N. */

/* ***SEE ALSO  DNSQ, DNSQE */
/* ***ROUTINES CALLED  D1MACH, DENORM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DDOGLG */
    /* Parameter adjustments */
    --wa2;
    --wa1;
    --x;
    --qtb;
    --diag;
    --r__;

    /* Function Body */

/*     EPSMCH IS THE MACHINE PRECISION. */

/* ***FIRST EXECUTABLE STATEMENT  DDOGLG */
    epsmch = d1mach_(&c__4);

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
	    d__2 = temp, d__3 = (d__1 = r__[l], abs(d__1));
	    temp = max(d__2,d__3);
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
    qnorm = denorm_(n, &wa2[1]);
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

/*     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR */
/*     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO. */

    gnorm = denorm_(n, &wa1[1]);
    sgnorm = zero;
    alpha = *delta / qnorm;
    if (gnorm == zero) {
	goto L120;
    }

/*     CALCULATE THE POINT ALONG THE SCALED GRADIENT */
/*     AT WHICH THE QUADRATIC IS MINIMIZED. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = wa1[j] / gnorm / diag[j];
/* L90: */
    }
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
    temp = denorm_(n, &wa2[1]);
    sgnorm = gnorm / temp / temp;

/*     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE. */

    alpha = zero;
    if (sgnorm >= *delta) {
	goto L120;
    }

/*     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE. */
/*     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG */
/*     AT WHICH THE QUADRATIC IS MINIMIZED. */

    bnorm = denorm_(n, &qtb[1]);
    temp = bnorm / gnorm * (bnorm / qnorm) * (sgnorm / *delta);
/* Computing 2nd power */
    d__1 = sgnorm / *delta;
/* Computing 2nd power */
    d__2 = temp - *delta / qnorm;
/* Computing 2nd power */
    d__3 = *delta / qnorm;
/* Computing 2nd power */
    d__4 = sgnorm / *delta;
    temp = temp - *delta / qnorm * (d__1 * d__1) + sqrt(d__2 * d__2 + (one - 
	    d__3 * d__3) * (one - d__4 * d__4));
/* Computing 2nd power */
    d__1 = sgnorm / *delta;
    alpha = *delta / qnorm * (one - d__1 * d__1) / temp;
L120:

/*     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON */
/*     DIRECTION AND THE SCALED GRADIENT DIRECTION. */

    temp = (one - alpha) * min(sgnorm,*delta);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = temp * wa1[j] + alpha * x[j];
/* L130: */
    }
L140:
    return 0;

/*     LAST CARD OF SUBROUTINE DDOGLG. */

} /* ddoglg_ */

