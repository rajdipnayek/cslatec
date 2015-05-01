/* ddntp.f -- translated by f2c (version 12.02.01).
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

/* DECK DDNTP */
/* Subroutine */ int ddntp_(doublereal *h__, integer *k, integer *n, integer *
	nq, doublereal *t, doublereal *tout, doublereal *yh, doublereal *y)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal r__;
    static integer jj, kk, kused;
    static doublereal factor;

/* ***BEGIN PROLOGUE  DDNTP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine DDNTP interpolates the K-th derivative of Y at */
/*            TOUT, using the data in the YH array.  If K has a value */
/*            greater than NQ, the NQ-th derivative is calculated. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      DOUBLE PRECISION (SDNTP-S, DDNTP-D, CDNTP-C) */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  DDNTP */
/* ***FIRST EXECUTABLE STATEMENT  DDNTP */
    /* Parameter adjustments */
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --y;

    /* Function Body */
    if (*k == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	    y[i__] = yh[i__ + (*nq + 1) * yh_dim1];
	}
	r__ = (*tout - *t) / *h__;
	i__1 = *nq;
	for (jj = 1; jj <= i__1; ++jj) {
	    j = *nq + 1 - jj;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
		y[i__] = yh[i__ + j * yh_dim1] + r__ * y[i__];
	    }
	}
    } else {
	kused = min(*k,*nq);
	factor = 1.;
	i__2 = kused;
	for (kk = 1; kk <= i__2; ++kk) {
/* L40: */
	    factor *= *nq + 1 - kk;
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L50: */
	    y[i__] = factor * yh[i__ + (*nq + 1) * yh_dim1];
	}
	r__ = (*tout - *t) / *h__;
	i__2 = *nq;
	for (jj = kused + 1; jj <= i__2; ++jj) {
	    j = kused + 1 + *nq - jj;
	    factor = 1.;
	    i__1 = kused;
	    for (kk = 1; kk <= i__1; ++kk) {
/* L60: */
		factor *= j - kk;
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
		y[i__] = factor * yh[i__ + j * yh_dim1] + r__ * y[i__];
	    }
/* L80: */
	}
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L100: */
	    i__1 = -kused;
	    y[i__] *= pow_di(h__, &i__1);
	}
    }
    return 0;
} /* ddntp_ */

