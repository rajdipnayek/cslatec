/* cdotc.f -- translated by f2c (version 12.02.01).
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

/* DECK CDOTC */
/* Complex */ void cdotc_(complex * ret_val, integer *n, complex *cx, integer 
	*incx, complex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__, ns, kx, ky;

/* ***BEGIN PROLOGUE  CDOTC */
/* ***PURPOSE  Dot product of two complex vectors using the complex */
/*            conjugate of the first vector. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A4 */
/* ***TYPE      COMPLEX (CDOTC-C) */
/* ***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR */
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
/*    CDOTC  complex result (zero if N .LE. 0) */

/*     Returns the dot product of complex CX and CY, using CONJUGATE(CX) */
/*     CDOTC = SUM for I = 0 to N-1 of CONJ(CX(LX+I*INCX))*CY(LY+I*INCY), */
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
/* ***END PROLOGUE  CDOTC */
/* ***FIRST EXECUTABLE STATEMENT  CDOTC */
    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
     ret_val->r = 0.f,  ret_val->i = 0.f;
    if (*n <= 0) {
	return ;
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
	r_cnjg(&q__3, &cx[kx]);
	i__2 = ky;
	q__2.r = q__3.r * cy[i__2].r - q__3.i * cy[i__2].i, q__2.i = q__3.r * 
		cy[i__2].i + q__3.i * cy[i__2].r;
	q__1.r =  ret_val->r + q__2.r, q__1.i =  ret_val->i + q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
	kx += *incx;
	ky += *incy;
/* L10: */
    }
    return ;

/*     Code for equal, positive increments. */

L20:
    ns = *n * *incx;
    i__1 = ns;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	r_cnjg(&q__3, &cx[i__]);
	i__3 = i__;
	q__2.r = q__3.r * cy[i__3].r - q__3.i * cy[i__3].i, q__2.i = q__3.r * 
		cy[i__3].i + q__3.i * cy[i__3].r;
	q__1.r =  ret_val->r + q__2.r, q__1.i =  ret_val->i + q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
/* L30: */
    }
    return ;
} /* cdotc_ */

