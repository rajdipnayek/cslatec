/* cdcdot.f -- translated by f2c (version 12.02.01).
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

/* DECK CDCDOT */
/* Complex */ void cdcdot_(complex * ret_val, integer *n, complex *cb, 
	complex *cx, integer *incx, complex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    complex q__1;

    /* Local variables */
    static integer i__, kx, ky;
    static doublereal dt1, dt2, dt3, dt4, dsdoti, dsdotr;

/* ***BEGIN PROLOGUE  CDCDOT */
/* ***PURPOSE  Compute the inner product of two vectors with extended */
/*            precision accumulation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A4 */
/* ***TYPE      COMPLEX (SDSDOT-S, CDCDOT-C) */
/* ***KEYWORDS  BLAS, DOT PRODUCT, INNER PRODUCT, LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       CB  complex scalar to be added to inner product */
/*       CX  complex vector with N elements */
/*     INCX  storage spacing between elements of CX */
/*       CY  complex vector with N elements */
/*     INCY  storage spacing between elements of CY */

/*     --Output-- */
/*   CDCDOT  complex dot product (CB if N .LE. 0) */

/*     Returns complex result with dot product accumulated in D.P. */
/*     CDCDOT = CB + sum for I = 0 to N-1 of CX(LX+I*INCY)*CY(LY+I*INCY) */
/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CDCDOT */
/* ***FIRST EXECUTABLE STATEMENT  CDCDOT */
    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    dsdotr = (doublereal) cb->r;
    dsdoti = (doublereal) r_imag(cb);
    if (*n <= 0) {
	goto L10;
    }
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
	i__2 = kx;
	dt1 = (doublereal) cx[i__2].r;
	i__2 = ky;
	dt2 = (doublereal) cy[i__2].r;
	dt3 = (doublereal) r_imag(&cx[kx]);
	dt4 = (doublereal) r_imag(&cy[ky]);
	dsdotr = dsdotr + dt1 * dt2 - dt3 * dt4;
	dsdoti = dsdoti + dt1 * dt4 + dt3 * dt2;
	kx += *incx;
	ky += *incy;
/* L5: */
    }
L10:
    r__1 = (real) dsdotr;
    r__2 = (real) dsdoti;
    q__1.r = r__1, q__1.i = r__2;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;
} /* cdcdot_ */

