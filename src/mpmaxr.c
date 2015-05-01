/* mpmaxr.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c__4 = 4;

/* DECK MPMAXR */
/* Subroutine */ int mpmaxr_(integer *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, it;
    extern /* Subroutine */ int mpchk_(integer *, integer *);

/* ***BEGIN PROLOGUE  MPMAXR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPMAXR-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Sets X to the largest possible positive 'mp' number. */

/*  The argument X(*) is an INTEGER arrays of size 30.  See the comments */
/*  in the routine MPBLAS for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPCHK */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPMAXR */
/* ***FIRST EXECUTABLE STATEMENT  MPMAXR */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    mpchk_(&c__1, &c__4);
    it = mpcom_1.b - 1;
/* SET FRACTION DIGITS TO B-1 */
    i__1 = mpcom_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	x[i__ + 2] = it;
    }
/* SET SIGN AND EXPONENT */
    x[1] = 1;
    x[2] = mpcom_1.m;
    return 0;
} /* mpmaxr_ */

