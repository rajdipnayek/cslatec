/* drotm.f -- translated by f2c (version 12.02.01).
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

/* DECK DROTM */
/* Subroutine */ int drotm_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal w, z__;
    static integer kx, ky;
    static doublereal dh11, dh12, dh22, dh21, dflag;
    static integer nsteps;

/* ***BEGIN PROLOGUE  DROTM */
/* ***PURPOSE  Apply a modified Givens transformation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A8 */
/* ***TYPE      DOUBLE PRECISION (SROTM-S, DROTM-D) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR */
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
/*   DPARAM  5-element D.P. vector.  DPARAM(1) is DFLAG described below. */
/*           Locations 2-5 of SPARAM contain elements of the */
/*           transformation matrix H described below. */

/*     --Output-- */
/*       DX  rotated vector (unchanged if N .LE. 0) */
/*       DY  rotated vector (unchanged if N .LE. 0) */

/*     Apply the modified Givens transformation, H, to the 2 by N matrix */
/*     (DX**T) */
/*     (DY**T) , where **T indicates transpose.  The elements of DX are */
/*     in DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else */
/*     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY. */

/*     With DPARAM(1)=DFLAG, H has one of the following forms: */

/*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */

/*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */

/*     See DROTMG for a description of data storage in DPARAM. */

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
/* ***END PROLOGUE  DROTM */
    /* Parameter adjustments */
    --dparam;
    --dy;
    --dx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DROTM */
    dflag = dparam[1];
    if (*n <= 0 || dflag + two == zero) {
	goto L140;
    }
    if (! (*incx == *incy && *incx > 0)) {
	goto L70;
    }

    nsteps = *n * *incx;
    if (dflag < 0.) {
	goto L50;
    } else if (dflag == 0) {
	goto L10;
    } else {
	goto L30;
    }
L10:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = dx[i__];
	z__ = dy[i__];
	dx[i__] = w + z__ * dh12;
	dy[i__] = w * dh21 + z__;
/* L20: */
    }
    goto L140;
L30:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = nsteps;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
	w = dx[i__];
	z__ = dy[i__];
	dx[i__] = w * dh11 + z__;
	dy[i__] = -w + dh22 * z__;
/* L40: */
    }
    goto L140;
L50:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = dx[i__];
	z__ = dy[i__];
	dx[i__] = w * dh11 + z__ * dh12;
	dy[i__] = w * dh21 + z__ * dh22;
/* L60: */
    }
    goto L140;
L70:
    kx = 1;
    ky = 1;
    if (*incx < 0) {
	kx = (1 - *n) * *incx + 1;
    }
    if (*incy < 0) {
	ky = (1 - *n) * *incy + 1;
    }

    if (dflag < 0.) {
	goto L120;
    } else if (dflag == 0) {
	goto L80;
    } else {
	goto L100;
    }
L80:
    dh12 = dparam[4];
    dh21 = dparam[3];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = dx[kx];
	z__ = dy[ky];
	dx[kx] = w + z__ * dh12;
	dy[ky] = w * dh21 + z__;
	kx += *incx;
	ky += *incy;
/* L90: */
    }
    goto L140;
L100:
    dh11 = dparam[2];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = dx[kx];
	z__ = dy[ky];
	dx[kx] = w * dh11 + z__;
	dy[ky] = -w + dh22 * z__;
	kx += *incx;
	ky += *incy;
/* L110: */
    }
    goto L140;
L120:
    dh11 = dparam[2];
    dh12 = dparam[4];
    dh21 = dparam[3];
    dh22 = dparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = dx[kx];
	z__ = dy[ky];
	dx[kx] = w * dh11 + z__ * dh12;
	dy[ky] = w * dh21 + z__ * dh22;
	kx += *incx;
	ky += *incy;
/* L130: */
    }
L140:
    return 0;
} /* drotm_ */

