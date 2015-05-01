/* dcopym.f -- translated by f2c (version 12.02.01).
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

/* DECK DCOPYM */
/* Subroutine */ int dcopym_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, ix, iy, ns, mp1;

/* ***BEGIN PROLOGUE  DCOPYM */
/* ***PURPOSE  Copy the negative of a vector to a vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A5 */
/* ***TYPE      DOUBLE PRECISION (SCOPYM-S, DCOPYM-D) */
/* ***KEYWORDS  BLAS, COPY, VECTOR */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/* ***DESCRIPTION */

/*       Description of Parameters */
/*           The * Flags Output Variables */

/*       N   Number of elements in vector(s) */
/*      DX   Double precision vector with N elements */
/*    INCX   Storage spacing between elements of DX */
/*      DY*  Double precision negative copy of DX */
/*    INCY   Storage spacing between elements of DY */

/*      ***  Note that DY = -DX  *** */

/*     Copy negative of d.p. DX to d.p. DY.  For I=0 to N-1, */
/*     copy  -DX(LX+I*INCX) to DY(LY+I*INCY), where LX=1 if */
/*     INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is defined */
/*     in a similar way using INCY. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/* ***END PROLOGUE  DCOPYM */
/* ***FIRST EXECUTABLE STATEMENT  DCOPYM */
    /* Parameter adjustments */
    --dy;
    --dx;

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
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = -dx[ix];
	ix += *incx;
	iy += *incy;
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
	dy[i__] = -dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = -dx[i__];
	dy[i__ + 1] = -dx[i__ + 1];
	dy[i__ + 2] = -dx[i__ + 2];
	dy[i__ + 3] = -dx[i__ + 3];
	dy[i__ + 4] = -dx[i__ + 4];
	dy[i__ + 5] = -dx[i__ + 5];
	dy[i__ + 6] = -dx[i__ + 6];
/* L50: */
    }
    return 0;

/*     Code for equal, positive, non-unit increments. */

L60:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dy[i__] = -dx[i__];
/* L70: */
    }
    return 0;
} /* dcopym_ */

