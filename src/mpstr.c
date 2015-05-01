/* mpstr.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    integer b, t, m, lun, mxr, r__[30];
} mpcom_;

#define mpcom_1 mpcom_

/* DECK MPSTR */
/* Subroutine */ int mpstr_(integer *x, integer *y)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  MPSTR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPSTR-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Sets Y = X for 'mp' X and Y. */

/*  The arguments X(*) and Y(*) are INTEGER arrays of size 30.  See the */
/*  comments in the routine MPBLAS for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPSTR */
/* ***FIRST EXECUTABLE STATEMENT  MPSTR */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    i__1 = mpcom_1.t + 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = x[i__];
/* L10: */
    }
    return 0;
} /* mpstr_ */

