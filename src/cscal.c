/* cscal.f -- translated by f2c (version 12.02.01).
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

/* DECK CSCAL */
/* Subroutine */ int cscal_(integer *n, complex *ca, complex *cx, integer *
	incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    complex q__1;

    /* Local variables */
    static integer i__, ix;

/* ***BEGIN PROLOGUE  CSCAL */
/* ***PURPOSE  Multiply a vector by a constant. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A6 */
/* ***TYPE      COMPLEX (SSCAL-S, DSCAL-D, CSCAL-C) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       CA  complex scale factor */
/*       CX  complex vector with N elements */
/*     INCX  storage spacing between elements of CX */

/*     --Output-- */
/*       CX  complex result (unchanged if N .LE. 0) */

/*     Replace complex CX by complex CA*CX. */
/*     For I = 0 to N-1, replace CX(IX+I*INCX) with CA*CX(IX+I*INCX), */
/*     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */

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
/*   900821  Modified to correct problem with a negative increment. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSCAL */
/* ***FIRST EXECUTABLE STATEMENT  CSCAL */
    /* Parameter adjustments */
    --cx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }

    if (*incx == 1) {
	goto L20;
    }

/*     Code for increment not equal to 1. */

    ix = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	i__3 = ix;
	q__1.r = ca->r * cx[i__3].r - ca->i * cx[i__3].i, q__1.i = ca->r * cx[
		i__3].i + ca->i * cx[i__3].r;
	cx[i__2].r = q__1.r, cx[i__2].i = q__1.i;
	ix += *incx;
/* L10: */
    }
    return 0;

/*     Code for increment equal to 1. */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	q__1.r = ca->r * cx[i__3].r - ca->i * cx[i__3].i, q__1.i = ca->r * cx[
		i__3].i + ca->i * cx[i__3].r;
	cx[i__2].r = q__1.r, cx[i__2].i = q__1.i;
/* L30: */
    }
    return 0;
} /* cscal_ */

