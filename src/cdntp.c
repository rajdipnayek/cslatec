/* cdntp.f -- translated by f2c (version 12.02.01).
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

/* DECK CDNTP */
/* Subroutine */ int cdntp_(real *h__, integer *k, integer *n, integer *nq, 
	real *t, real *tout, complex *yh, complex *y)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__, j;
    static real r__;
    static integer jj, kk, kused;
    static real factor;

/* ***BEGIN PROLOGUE  CDNTP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine CDNTP interpolates the K-th derivative of Y at */
/*            TOUT, using the data in the YH array.  If K has a value */
/*            greater than NQ, the NQ-th derivative is calculated. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      COMPLEX (SDNTP-S, DDNTP-D, CDNTP-C) */
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
/* ***END PROLOGUE  CDNTP */
/* ***FIRST EXECUTABLE STATEMENT  CDNTP */
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
	    i__2 = i__;
	    i__3 = i__ + (*nq + 1) * yh_dim1;
	    y[i__2].r = yh[i__3].r, y[i__2].i = yh[i__3].i;
	}
	r__ = (*tout - *t) / *h__;
	i__2 = *nq;
	for (jj = 1; jj <= i__2; ++jj) {
	    j = *nq + 1 - jj;
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L20: */
		i__1 = i__;
		i__4 = i__ + j * yh_dim1;
		i__5 = i__;
		q__2.r = r__ * y[i__5].r, q__2.i = r__ * y[i__5].i;
		q__1.r = yh[i__4].r + q__2.r, q__1.i = yh[i__4].i + q__2.i;
		y[i__1].r = q__1.r, y[i__1].i = q__1.i;
	    }
	}
    } else {
	kused = min(*k,*nq);
	factor = 1.f;
	i__1 = kused;
	for (kk = 1; kk <= i__1; ++kk) {
/* L40: */
	    factor *= *nq + 1 - kk;
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	    i__4 = i__;
	    i__5 = i__ + (*nq + 1) * yh_dim1;
	    q__1.r = factor * yh[i__5].r, q__1.i = factor * yh[i__5].i;
	    y[i__4].r = q__1.r, y[i__4].i = q__1.i;
	}
	r__ = (*tout - *t) / *h__;
	i__4 = *nq;
	for (jj = kused + 1; jj <= i__4; ++jj) {
	    j = kused + 1 + *nq - jj;
	    factor = 1.f;
	    i__5 = kused;
	    for (kk = 1; kk <= i__5; ++kk) {
/* L60: */
		factor *= j - kk;
	    }
	    i__5 = *n;
	    for (i__ = 1; i__ <= i__5; ++i__) {
/* L70: */
		i__1 = i__;
		i__3 = i__ + j * yh_dim1;
		q__2.r = factor * yh[i__3].r, q__2.i = factor * yh[i__3].i;
		i__2 = i__;
		q__3.r = r__ * y[i__2].r, q__3.i = r__ * y[i__2].i;
		q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
		y[i__1].r = q__1.r, y[i__1].i = q__1.i;
	    }
/* L80: */
	}
	i__4 = *n;
	for (i__ = 1; i__ <= i__4; ++i__) {
/* L100: */
	    i__1 = i__;
	    i__3 = i__;
	    i__2 = -kused;
	    r__1 = pow_ri(h__, &i__2);
	    q__1.r = r__1 * y[i__3].r, q__1.i = r__1 * y[i__3].i;
	    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
	}
    }
    return 0;
} /* cdntp_ */

