/* cdpsc.f -- translated by f2c (version 12.02.01).
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

/* DECK CDPSC */
/* Subroutine */ int cdpsc_(integer *ksgn, integer *n, integer *nq, complex *
	yh)
{
    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1;

    /* Local variables */
    static integer i__, j, j1, j2;

/* ***BEGIN PROLOGUE  CDPSC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine CDPSC computes the predicted YH values by */
/*            effectively multiplying the YH array by the Pascal triangle */
/*            matrix when KSGN is +1, and performs the inverse function */
/*            when KSGN is -1. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      COMPLEX (SDPSC-S, DDPSC-D, CDPSC-C) */
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
/* ***END PROLOGUE  CDPSC */
/* ***FIRST EXECUTABLE STATEMENT  CDPSC */
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
		    i__4 = i__ + j * yh_dim1;
		    i__5 = i__ + j * yh_dim1;
		    i__6 = i__ + (j + 1) * yh_dim1;
		    q__1.r = yh[i__5].r + yh[i__6].r, q__1.i = yh[i__5].i + 
			    yh[i__6].i;
		    yh[i__4].r = q__1.r, yh[i__4].i = q__1.i;
		}
	    }
	}
    } else {
	i__4 = *nq;
	for (j1 = 1; j1 <= i__4; ++j1) {
	    i__5 = *nq;
	    for (j2 = j1; j2 <= i__5; ++j2) {
		j = *nq - j2 + j1;
		i__6 = *n;
		for (i__ = 1; i__ <= i__6; ++i__) {
/* L30: */
		    i__3 = i__ + j * yh_dim1;
		    i__2 = i__ + j * yh_dim1;
		    i__1 = i__ + (j + 1) * yh_dim1;
		    q__1.r = yh[i__2].r - yh[i__1].r, q__1.i = yh[i__2].i - 
			    yh[i__1].i;
		    yh[i__3].r = q__1.r, yh[i__3].i = q__1.i;
		}
	    }
	}
    }
    return 0;
} /* cdpsc_ */

