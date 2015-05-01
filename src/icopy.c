/* icopy.f -- translated by f2c (version 12.02.01).
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

/* DECK ICOPY */
/* Subroutine */ int icopy_(integer *n, integer *ix, integer *incx, integer *
	iy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, ns, mp1, iix, iiy;

/* ***BEGIN PROLOGUE  ICOPY */
/* ***PURPOSE  Copy a vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A5 */
/* ***TYPE      INTEGER (ICOPY-S, DCOPY-D, CCOPY-C, ICOPY-I) */
/* ***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Boland, W. Robert, (LANL) */
/*           Clemens, Reginald, (PLK) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       IX  integer vector with N elements */
/*     INCX  storage spacing between elements of IX */
/*       IY  integer vector with N elements */
/*     INCY  storage spacing between elements of IY */

/*     --Output-- */
/*       IY  copy of vector IX (unchanged if N .LE. 0) */

/*     Copy integer IX to integer IY. */
/*     For I = 0 to N-1, copy  IX(LX+I*INCX) to IY(LY+I*INCY), */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   930201  DATE WRITTEN */
/* ***END PROLOGUE  ICOPY */
/* ***FIRST EXECUTABLE STATEMENT  ICOPY */
    /* Parameter adjustments */
    --iy;
    --ix;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == *incy) {
	if ((i__1 = *incx - 1) < 0) {
	    goto L5;
	} else if (i__1 == 0) {
	    goto L20;
	} else {
	    goto L60;
	}
    }

/*     Code for unequal or nonpositive increments. */

L5:
    iix = 1;
    iiy = 1;
    if (*incx < 0) {
	iix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iiy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iy[iiy] = ix[iix];
	iix += *incx;
	iiy += *incy;
/* L10: */
    }
    return 0;

/*     Code for both increments equal to 1. */

/*     Clean-up loop so remaining vector length is a multiple of 7. */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iy[i__] = ix[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	iy[i__] = ix[i__];
	iy[i__ + 1] = ix[i__ + 1];
	iy[i__ + 2] = ix[i__ + 2];
	iy[i__ + 3] = ix[i__ + 3];
	iy[i__ + 4] = ix[i__ + 4];
	iy[i__ + 5] = ix[i__ + 5];
	iy[i__ + 6] = ix[i__ + 6];
/* L50: */
    }
    return 0;

/*     Code for equal, positive, non-unit increments. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	iy[i__] = ix[i__];
/* L70: */
    }
    return 0;
} /* icopy_ */

