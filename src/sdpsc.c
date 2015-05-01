/* sdpsc.f -- translated by f2c (version 12.02.01).
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

/* DECK SDPSC */
/* Subroutine */ int sdpsc_(integer *ksgn, integer *n, integer *nq, real *yh)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, j1, j2;

/* ***BEGIN PROLOGUE  SDPSC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine SDPSC computes the predicted YH values by */
/*            effectively multiplying the YH array by the Pascal triangle */
/*            matrix when KSGN is +1, and performs the inverse function */
/*            when KSGN is -1. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      SINGLE PRECISION (SDPSC-S, DDPSC-D, CDPSC-C) */
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
/* ***END PROLOGUE  SDPSC */
/* ***FIRST EXECUTABLE STATEMENT  SDPSC */
    /* Parameter adjustments */
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;

    /* Function Body */
    if (*ksgn > 0) {
	i__1 = *nq;
	for (j1 = 1; j1 <= i__1; ++j1) {
	    i__2 = *nq;
	    for (j2 = j1; j2 <= i__2; ++j2) {
		j = *nq - j2 + j1;
		i__3 = *n;
		for (i__ = 1; i__ <= i__3; ++i__) {
/* L10: */
		    yh[i__ + j * yh_dim1] += yh[i__ + (j + 1) * yh_dim1];
		}
	    }
	}
    } else {
	i__3 = *nq;
	for (j1 = 1; j1 <= i__3; ++j1) {
	    i__2 = *nq;
	    for (j2 = j1; j2 <= i__2; ++j2) {
		j = *nq - j2 + j1;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
		    yh[i__ + j * yh_dim1] -= yh[i__ + (j + 1) * yh_dim1];
		}
	    }
	}
    }
    return 0;
} /* sdpsc_ */

