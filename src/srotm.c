/* srotm.f -- translated by f2c (version 12.02.01).
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

/* DECK SROTM */
/* Subroutine */ int srotm_(integer *n, real *sx, integer *incx, real *sy, 
	integer *incy, real *sparam)
{
    /* Initialized data */

    static real zero = 0.f;
    static real two = 2.f;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__;
    static real w, z__;
    static integer kx, ky;
    static real sh11, sh12, sh21, sh22, sflag;
    static integer nsteps;

/* ***BEGIN PROLOGUE  SROTM */
/* ***PURPOSE  Apply a modified Givens transformation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1A8 */
/* ***TYPE      SINGLE PRECISION (SROTM-S, DROTM-D) */
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
/*       SX  single precision vector with N elements */
/*     INCX  storage spacing between elements of SX */
/*       SY  single precision vector with N elements */
/*     INCY  storage spacing between elements of SY */
/*   SPARAM  5-element vector. SPARAM(1) is SFLAG described below. */
/*           Locations 2-5 of SPARAM contain elements of the */
/*           transformation matrix H described below. */

/*     --Output-- */
/*       SX  rotated vector (unchanged if N .LE. 0) */
/*       SY  rotated vector (unchanged if N .LE. 0) */

/*     Apply the modified Givens transformation, H, to the 2 by N matrix */
/*     (SX**T) */
/*     (SY**T) , where **T indicates transpose.  The elements of SX are */
/*     in SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX .GE. 0, else */
/*     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY. */

/*     With SPARAM(1)=SFLAG, H has one of the following forms: */

/*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0 */

/*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0). */

/*     See SROTMG for a description of data storage in SPARAM. */

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
/* ***END PROLOGUE  SROTM */
    /* Parameter adjustments */
    --sparam;
    --sy;
    --sx;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  SROTM */
    sflag = sparam[1];
    if (*n <= 0 || sflag + two == zero) {
	goto L140;
    }
    if (! (*incx == *incy && *incx > 0)) {
	goto L70;
    }

    nsteps = *n * *incx;
    if (sflag < 0.f) {
	goto L50;
    } else if (sflag == 0) {
	goto L10;
    } else {
	goto L30;
    }
L10:
    sh12 = sparam[4];
    sh21 = sparam[3];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = sx[i__];
	z__ = sy[i__];
	sx[i__] = w + z__ * sh12;
	sy[i__] = w * sh21 + z__;
/* L20: */
    }
    goto L140;
L30:
    sh11 = sparam[2];
    sh22 = sparam[5];
    i__2 = nsteps;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
	w = sx[i__];
	z__ = sy[i__];
	sx[i__] = w * sh11 + z__;
	sy[i__] = -w + sh22 * z__;
/* L40: */
    }
    goto L140;
L50:
    sh11 = sparam[2];
    sh12 = sparam[4];
    sh21 = sparam[3];
    sh22 = sparam[5];
    i__1 = nsteps;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	w = sx[i__];
	z__ = sy[i__];
	sx[i__] = w * sh11 + z__ * sh12;
	sy[i__] = w * sh21 + z__ * sh22;
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

    if (sflag < 0.f) {
	goto L120;
    } else if (sflag == 0) {
	goto L80;
    } else {
	goto L100;
    }
L80:
    sh12 = sparam[4];
    sh21 = sparam[3];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = sx[kx];
	z__ = sy[ky];
	sx[kx] = w + z__ * sh12;
	sy[ky] = w * sh21 + z__;
	kx += *incx;
	ky += *incy;
/* L90: */
    }
    goto L140;
L100:
    sh11 = sparam[2];
    sh22 = sparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = sx[kx];
	z__ = sy[ky];
	sx[kx] = w * sh11 + z__;
	sy[ky] = -w + sh22 * z__;
	kx += *incx;
	ky += *incy;
/* L110: */
    }
    goto L140;
L120:
    sh11 = sparam[2];
    sh12 = sparam[4];
    sh21 = sparam[3];
    sh22 = sparam[5];
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	w = sx[kx];
	z__ = sy[ky];
	sx[kx] = w * sh11 + z__ * sh12;
	sy[ky] = w * sh21 + z__ * sh22;
	kx += *incx;
	ky += *incy;
/* L130: */
    }
L140:
    return 0;
} /* srotm_ */

