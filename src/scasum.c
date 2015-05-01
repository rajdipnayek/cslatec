/* scasum.f -- translated by f2c (version 12.02.01).
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

/* DECK SCASUM */
doublereal scasum_(integer *n, complex *cx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2;

    /* Local variables */
    static integer i__, ix;

/* ***BEGIN PROLOGUE  SCASUM */
/* ***PURPOSE  Compute the sum of the magnitudes of the real and */
/*            imaginary elements of a complex vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A3A */
/* ***TYPE      COMPLEX (SASUM-S, DASUM-D, SCASUM-C) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR */
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
/*   SCASUM  single precision result (zero if N .LE. 0) */

/*     Returns sums of magnitudes of real and imaginary parts of */
/*     components of CX.  Note that this is not the L1 norm of CX. */
/*     CASUM = sum from 0 to N-1 of ABS(REAL(CX(IX+I*INCX))) + */
/*             ABS(IMAG(CX(IX+I*INCX))), */
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
/* ***END PROLOGUE  SCASUM */
/* ***FIRST EXECUTABLE STATEMENT  SCASUM */
    /* Parameter adjustments */
    --cx;

    /* Function Body */
    ret_val = 0.f;
    if (*n <= 0) {
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
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	ret_val = ret_val + (r__1 = cx[i__2].r, dabs(r__1)) + (r__2 = r_imag(&
		cx[ix]), dabs(r__2));
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*     Code for increment equal to 1. */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	ret_val = ret_val + (r__1 = cx[i__2].r, dabs(r__1)) + (r__2 = r_imag(&
		cx[i__]), dabs(r__2));
/* L30: */
    }
    return ret_val;
} /* scasum_ */

