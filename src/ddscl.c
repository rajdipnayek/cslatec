/* ddscl.f -- translated by f2c (version 12.02.01).
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

/* DECK DDSCL */
/* Subroutine */ int ddscl_(doublereal *hmax, integer *n, integer *nq, 
	doublereal *rmax, doublereal *h__, doublereal *rc, doublereal *rh, 
	doublereal *yh)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal r1;

/* ***BEGIN PROLOGUE  DDSCL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine DDSCL rescales the YH array whenever the step */
/*            size is changed. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      DOUBLE PRECISION (SDSCL-S, DDSCL-D, CDSCL-C) */
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
/* ***END PROLOGUE  DDSCL */
/* ***FIRST EXECUTABLE STATEMENT  DDSCL */
    /* Parameter adjustments */
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;

    /* Function Body */
    if (*h__ < 1.) {
/* Computing MIN */
	d__1 = abs(*h__) * *rh, d__2 = abs(*h__) * *rmax, d__1 = min(d__1,
		d__2);
	*rh = min(d__1,*hmax) / abs(*h__);
    } else {
/* Computing MIN */
	d__1 = min(*rh,*rmax), d__2 = *hmax / abs(*h__);
	*rh = min(d__1,d__2);
    }
    r1 = 1.;
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
} /* ddscl_ */

