/* ccopy.f -- translated by f2c (version 12.02.01).
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

/* DECK CCOPY */
/* Subroutine */ int ccopy_(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, ns, kx, ky;

/* ***BEGIN PROLOGUE  CCOPY */
/* ***PURPOSE  Copy a vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A5 */
/* ***TYPE      COMPLEX (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I) */
/* ***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       CX  complex vector with N elements */
/*     INCX  storage spacing between elements of CX */
/*       CY  complex vector with N elements */
/*     INCY  storage spacing between elements of CY */

/*     --Output-- */
/*       CY  copy of vector CX (unchanged if N .LE. 0) */

/*     Copy complex CX to complex CY. */
/*     For I = 0 to N-1, copy CX(LX+I*INCX) to CY(LY+I*INCY), */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CCOPY */
/* ***FIRST EXECUTABLE STATEMENT  CCOPY */
    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == *incy && *incx > 0) {
	goto L20;
    }

/*     Code for unequal or nonpositive increments. */

    kx = 1;
    ky = 1;
    if (*incx < 0) {
	kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
	ky = (1 - *n) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ky;
	i__3 = kx;
	cy[i__2].r = cx[i__3].r, cy[i__2].i = cx[i__3].i;
	kx += *incx;
	ky += *incy;
/* L10: */
    }
    return 0;

/*     Code for equal, positive increments. */

L20:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	i__3 = i__;
	i__4 = i__;
	cy[i__3].r = cx[i__4].r, cy[i__3].i = cx[i__4].i;
/* L30: */
    }
    return 0;
} /* ccopy_ */

