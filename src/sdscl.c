/* sdscl.f -- translated by f2c (version 12.02.01).
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

/* DECK SDSCL */
/* Subroutine */ int sdscl_(real *hmax, integer *n, integer *nq, real *rmax, 
	real *h__, real *rc, real *rh, real *yh)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j;
    static real r1;

/* ***BEGIN PROLOGUE  SDSCL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine SDSCL rescales the YH array whenever the step */
/*            size is changed. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      SINGLE PRECISION (SDSCL-S, DDSCL-D, CDSCL-C) */
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
/* ***END PROLOGUE  SDSCL */
/* ***FIRST EXECUTABLE STATEMENT  SDSCL */
    /* Parameter adjustments */
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;

    /* Function Body */
    if (*h__ < 1.f) {
/* Computing MIN */
	r__1 = dabs(*h__) * *rh, r__2 = dabs(*h__) * *rmax, r__1 = min(r__1,
		r__2);
	*rh = dmin(r__1,*hmax) / dabs(*h__);
    } else {
/* Computing MIN */
	r__1 = min(*rh,*rmax), r__2 = *hmax / dabs(*h__);
	*rh = dmin(r__1,r__2);
    }
    r1 = 1.f;
    i__1 = *nq;
    for (j = 1; j <= i__1; ++j) {
	r1 *= *rh;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L10: */
	    yh[i__ + (j + 1) * yh_dim1] *= r1;
	}
    }
    *h__ *= *rh;
    *rc *= *rh;
    return 0;
} /* sdscl_ */

