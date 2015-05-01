/* dcdot.f -- translated by f2c (version 12.02.01).
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

/* DECK DCDOT */
/* Subroutine */ int dcdot_(integer *n, doublereal *fm, complex *cx, integer *
	incx, complex *cy, integer *incy, doublereal *dcr, doublereal *dci)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, kx, ky;
    static doublereal dt1, dt2, dt3, dt4;

/* ***BEGIN PROLOGUE  DCDOT */
/* ***PURPOSE  Compute the inner product of two vectors with extended */
/*            precision accumulation and result. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A4 */
/* ***TYPE      COMPLEX (DSDOT-D, DCDOT-C) */
/* ***KEYWORDS  BLAS, COMPLEX VECTORS, DOT PRODUCT, INNER PRODUCT, */
/*             LINEAR ALGEBRA, VECTOR */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*    Compute the dot product of 2 complex vectors, CX and CY, e.g. */
/*    CX DOT CY, or, CXconjugate DOT CY.  The real and imaginary */
/*    parts of CX and CY are converted to double precision, the dot */
/*    product accumulation is done in double precision and the output */
/*    is given as 2 double precision numbers, corresponding to the real */
/*    and imaginary part of the result. */
/*     Input */
/*      N:  Number of complex components of CX and CY. */
/*      FM: =+1.0   compute CX DOT CY. */
/*          =-1.0   compute CXconjugate DOT CY. */
/*      CX(N): */
/*      CY(N):  Complex arrays of length N. */
/*      INCX:(Integer)   Spacing of elements of CX to use */
/*      INCY:(Integer)   Spacing of elements of CY to use. */
/*     Output */
/*      DCR:(Double Precision) Real part of dot product. */
/*      DCI:(Double Precision) Imaginary part of dot product. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DCDOT */
/* ***FIRST EXECUTABLE STATEMENT  DCDOT */
    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    *dcr = 0.;
    *dci = 0.;
    if (*n <= 0) {
	goto L20;
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
	*dcr = *dcr + dt1 * dt2 - *fm * (dt3 * dt4);
	*dci = *dci + dt1 * dt4 + *fm * (dt3 * dt2);
	kx += *incx;
	ky += *incy;
/* L10: */
    }
L20:
    return 0;
} /* dcdot_ */

