/* csrot.f -- translated by f2c (version 12.02.01).
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

/* DECK CSROT */
/* Subroutine */ int csrot_(integer *n, complex *cx, integer *incx, complex *
	cy, integer *incy, real *c__, real *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer i__, ix, iy;
    static complex ctemp;

/* ***BEGIN PROLOGUE  CSROT */
/* ***PURPOSE  Apply a plane Givens rotation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B10 */
/* ***TYPE      COMPLEX (SROT-S, DROT-D, CSROT-C) */
/* ***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION, */
/*             LINEAR ALGEBRA, PLANE ROTATION, VECTOR */
/* ***AUTHOR  Dongarra, J., (ANL) */
/* ***DESCRIPTION */

/*     CSROT applies the complex Givens rotation */

/*          (X)   ( C S)(X) */
/*          (Y) = (-S C)(Y) */

/*     N times where for I = 0,...,N-1 */

/*          X = CX(LX+I*INCX) */
/*          Y = CY(LY+I*INCY), */

/*     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is */
/*     defined in a similar way using INCY. */

/*     Argument Description */

/*        N      (integer)  number of elements in each vector */

/*        CX     (complex array)  beginning of one vector */

/*        INCX   (integer)  memory spacing of successive elements */
/*               of vector CX */

/*        CY     (complex array)  beginning of the other vector */

/*        INCY   (integer)  memory spacing of successive elements */
/*               of vector CY */

/*        C      (real)  cosine term of the rotation */

/*        S      (real)  sine term of the rotation. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810223  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CSROT */
/* ***FIRST EXECUTABLE STATEMENT  CSROT */
    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*     Code for unequal increments or equal increments not equal to 1. */

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
	i__2 = ix;
	q__2.r = *c__ * cx[i__2].r, q__2.i = *c__ * cx[i__2].i;
	i__3 = iy;
	q__3.r = *s * cy[i__3].r, q__3.i = *s * cy[i__3].i;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
	i__2 = iy;
	i__3 = iy;
	q__2.r = *c__ * cy[i__3].r, q__2.i = *c__ * cy[i__3].i;
	i__4 = ix;
	q__3.r = *s * cx[i__4].r, q__3.i = *s * cx[i__4].i;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	cy[i__2].r = q__1.r, cy[i__2].i = q__1.i;
	i__2 = ix;
	cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*     Code for both increments equal to 1. */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	q__2.r = *c__ * cx[i__2].r, q__2.i = *c__ * cx[i__2].i;
	i__3 = i__;
	q__3.r = *s * cy[i__3].r, q__3.i = *s * cy[i__3].i;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	ctemp.r = q__1.r, ctemp.i = q__1.i;
	i__2 = i__;
	i__3 = i__;
	q__2.r = *c__ * cy[i__3].r, q__2.i = *c__ * cy[i__3].i;
	i__4 = i__;
	q__3.r = *s * cx[i__4].r, q__3.i = *s * cx[i__4].i;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	cy[i__2].r = q__1.r, cy[i__2].i = q__1.i;
	i__2 = i__;
	cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
/* L30: */
    }
    return 0;
} /* csrot_ */

