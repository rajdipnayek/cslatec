/* icamax.f -- translated by f2c (version 12.02.01).
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

/* DECK ICAMAX */
integer icamax_(integer *n, complex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__, ix;
    static real xmag, smax;

/* ***BEGIN PROLOGUE  ICAMAX */
/* ***PURPOSE  Find the smallest index of the component of a complex */
/*            vector having the maximum sum of magnitudes of real */
/*            and imaginary parts. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A2 */
/* ***TYPE      COMPLEX (ISAMAX-S, IDAMAX-D, ICAMAX-C) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR */
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

/*     --Output-- */
/*   ICAMAX  smallest index (zero if N .LE. 0) */

/*     Returns the smallest index of the component of CX having the */
/*     largest sum of magnitudes of real and imaginary parts. */
/*     ICAMAX = first I, I = 1 to N, to maximize */
/*     ABS(REAL(CX(IX+(I-1)*INCX))) + ABS(IMAG(CX(IX+(I-1)*INCX))), */
/*     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900821  Modified to correct problem with a negative increment. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  ICAMAX */
/* ***FIRST EXECUTABLE STATEMENT  ICAMAX */
    /* Parameter adjustments */
    --cx;

    /* Function Body */
    ret_val = 0;
    if (*n <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }

    if (*incx == 1) {
	goto L20;
    }

/*     Code for increment not equal to 1. */

    ix = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = ix;
    smax = (r__1 = cx[i__1].r, dabs(r__1)) + (r__2 = r_imag(&cx[ix]), dabs(
	    r__2));
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = ix;
	xmag = (r__1 = cx[i__2].r, dabs(r__1)) + (r__2 = r_imag(&cx[ix]), 
		dabs(r__2));
	if (xmag > smax) {
	    ret_val = i__;
	    smax = xmag;
	}
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*     Code for increment equal to 1. */

L20:
    smax = (r__1 = cx[1].r, dabs(r__1)) + (r__2 = r_imag(&cx[1]), dabs(r__2));
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__;
	xmag = (r__1 = cx[i__2].r, dabs(r__1)) + (r__2 = r_imag(&cx[i__]), 
		dabs(r__2));
	if (xmag > smax) {
	    ret_val = i__;
	    smax = xmag;
	}
/* L30: */
    }
    return ret_val;
} /* icamax_ */

