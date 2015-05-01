/* drot.f -- translated by f2c (version 12.02.01).
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

/* DECK DROT */
/* Subroutine */ int drot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *dc, doublereal *ds)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal w, z__;
    static integer kx, ky, nsteps;

/* ***BEGIN PROLOGUE  DROT */
/* ***PURPOSE  Apply a plane Givens rotation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A8 */
/* ***TYPE      DOUBLE PRECISION (SROT-S, DROT-D, CSROT-C) */
/* ***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION, */
/*             LINEAR ALGEBRA, PLANE ROTATION, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*        N  number of elements in input vector(s) */
/*       DX  double precision vector with N elements */
/*     INCX  storage spacing between elements of DX */
/*       DY  double precision vector with N elements */
/*     INCY  storage spacing between elements of DY */
/*       DC  D.P. element of rotation matrix */
/*       DS  D.P. element of rotation matrix */

/*     --Output-- */
/*       DX  rotated vector DX (unchanged if N .LE. 0) */
/*       DY  rotated vector DY (unchanged if N .LE. 0) */

/*     Multiply the 2 x 2 matrix  ( DC DS) times the 2 x N matrix (DX**T) */
/*                                (-DS DC)                        (DY**T) */
/*     where **T indicates transpose.  The elements of DX are in */
/*     DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else */
/*     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY. */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920310  Corrected definition of LX in DESCRIPTION.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DROT */
    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DROT */
    if (*n <= 0 || *ds == zero && *dc == one) {
	goto L40;
    }
    if (! (*incx == *incy && *incx > 0)) {
	goto L20;
    }

/*          Code for equal and positive increments. */

    nsteps = *incx * *n;
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = dx[i__];
	z__ = dy[i__];
	dx[i__] = *dc * w + *ds * z__;
	dy[i__] = -(*ds) * w + *dc * z__;
/* L10: */
    }
    goto L40;

/*     Code for unequal or nonpositive increments. */

L20:
    kx = 1;
    ky = 1;

    if (*incx < 0) {
	kx = 1 - (*n - 1) * *incx;
    }
    if (*incy < 0) {
	ky = 1 - (*n - 1) * *incy;
    }

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = dx[kx];
	z__ = dy[ky];
	dx[kx] = *dc * w + *ds * z__;
	dy[ky] = -(*ds) * w + *dc * z__;
	kx += *incx;
	ky += *incy;
/* L30: */
    }
L40:

    return 0;
} /* drot_ */

