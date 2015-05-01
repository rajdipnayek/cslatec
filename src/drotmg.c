/* drotmg.f -- translated by f2c (version 12.02.01).
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

/* DECK DROTMG */
/* Subroutine */ int drotmg_(doublereal *dd1, doublereal *dd2, doublereal *
	dx1, doublereal *dy1, doublereal *dparam)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal gam = 4096.;
    static doublereal gamsq = 16777216.;
    static doublereal rgamsq = 5.9604645e-8;

    /* Format strings */
    static char fmt_120[] = "";
    static char fmt_150[] = "";
    static char fmt_180[] = "";
    static char fmt_210[] = "";

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal du, dp1, dp2, dq1, dq2, dh11, dh12, dh21, dh22;
    static integer igo;
    static doublereal dflag, dtemp;

    /* Assigned format variables */
    static char *igo_fmt;

/* ***BEGIN PROLOGUE  DROTMG */
/* ***PURPOSE  Construct a modified Givens transformation. */
/* ***LIBRARY   SLATEC (BLAS) */
/* ***CATEGORY  D1B10 */
/* ***TYPE      DOUBLE PRECISION (SROTMG-S, DROTMG-D) */
/* ***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR */
/* ***AUTHOR  Lawson, C. L., (JPL) */
/*           Hanson, R. J., (SNLA) */
/*           Kincaid, D. R., (U. of Texas) */
/*           Krogh, F. T., (JPL) */
/* ***DESCRIPTION */

/*                B L A S  Subprogram */
/*    Description of Parameters */

/*     --Input-- */
/*      DD1  double precision scalar */
/*      DD2  double precision scalar */
/*      DX1  double precision scalar */
/*      DX2  double precision scalar */
/*   DPARAM  D.P. 5-vector. DPARAM(1)=DFLAG defined below. */
/*           Locations 2-5 contain the rotation matrix. */

/*     --Output-- */
/*      DD1  changed to represent the effect of the transformation */
/*      DD2  changed to represent the effect of the transformation */
/*      DX1  changed to represent the effect of the transformation */
/*      DX2  unchanged */

/*     Construct the modified Givens transformation matrix H which zeros */
/*     the second component of the 2-vector  (SQRT(DD1)*DX1,SQRT(DD2)* */
/*     DY2)**T. */
/*     With DPARAM(1)=DFLAG, H has one of the following forms: */

/*     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0 */

/*       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0) */
/*     H=(          )    (          )    (          )    (          ) */
/*       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0). */

/*     Locations 2-5 of DPARAM contain DH11, DH21, DH12, and DH22, */
/*     respectively.  (Values of 1.D0, -1.D0, or 0.D0 implied by the */
/*     value of DPARAM(1) are not stored in DPARAM.) */

/* ***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T. */
/*                 Krogh, Basic linear algebra subprograms for Fortran */
/*                 usage, Algorithm No. 539, Transactions on Mathematical */
/*                 Software 5, 3 (September 1979), pp. 308-323. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920316  Prologue corrected.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DROTMG */
    /* Parameter adjustments */
    --dparam;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DROTMG */
    if (! (*dd1 < zero)) {
	goto L10;
    }
/*       GO ZERO-H-D-AND-DX1.. */
    goto L60;
L10:
/*     CASE-DD1-NONNEGATIVE */
    dp2 = *dd2 * *dy1;
    if (! (dp2 == zero)) {
	goto L20;
    }
    dflag = -two;
    goto L260;
/*     REGULAR-CASE.. */
L20:
    dp1 = *dd1 * *dx1;
    dq2 = dp2 * *dy1;
    dq1 = dp1 * *dx1;

    if (! (abs(dq1) > abs(dq2))) {
	goto L40;
    }
    dh21 = -(*dy1) / *dx1;
    dh12 = dp2 / dp1;

    du = one - dh12 * dh21;

    if (! (du <= zero)) {
	goto L30;
    }
/*         GO ZERO-H-D-AND-DX1.. */
    goto L60;
L30:
    dflag = zero;
    *dd1 /= du;
    *dd2 /= du;
    *dx1 *= du;
/*         GO SCALE-CHECK.. */
    goto L100;
L40:
    if (! (dq2 < zero)) {
	goto L50;
    }
/*         GO ZERO-H-D-AND-DX1.. */
    goto L60;
L50:
    dflag = one;
    dh11 = dp1 / dp2;
    dh22 = *dx1 / *dy1;
    du = one + dh11 * dh22;
    dtemp = *dd2 / du;
    *dd2 = *dd1 / du;
    *dd1 = dtemp;
    *dx1 = *dy1 * du;
/*         GO SCALE-CHECK */
    goto L100;
/*     PROCEDURE..ZERO-H-D-AND-DX1.. */
L60:
    dflag = -one;
    dh11 = zero;
    dh12 = zero;
    dh21 = zero;
    dh22 = zero;

    *dd1 = zero;
    *dd2 = zero;
    *dx1 = zero;
/*         RETURN.. */
    goto L220;
/*     PROCEDURE..FIX-H.. */
L70:
    if (! (dflag >= zero)) {
	goto L90;
    }

    if (! (dflag == zero)) {
	goto L80;
    }
    dh11 = one;
    dh22 = one;
    dflag = -one;
    goto L90;
L80:
    dh21 = -one;
    dh12 = one;
    dflag = -one;
L90:
    switch (igo) {
	case 0: goto L120;
	case 1: goto L150;
	case 2: goto L180;
	case 3: goto L210;
    }
/*     PROCEDURE..SCALE-CHECK */
L100:
L110:
    if (! (*dd1 <= rgamsq)) {
	goto L130;
    }
    if (*dd1 == zero) {
	goto L160;
    }
    igo = 0;
    igo_fmt = fmt_120;
/*              FIX-H.. */
    goto L70;
L120:
/* Computing 2nd power */
    d__1 = gam;
    *dd1 *= d__1 * d__1;
    *dx1 /= gam;
    dh11 /= gam;
    dh12 /= gam;
    goto L110;
L130:
L140:
    if (! (*dd1 >= gamsq)) {
	goto L160;
    }
    igo = 1;
    igo_fmt = fmt_150;
/*              FIX-H.. */
    goto L70;
L150:
/* Computing 2nd power */
    d__1 = gam;
    *dd1 /= d__1 * d__1;
    *dx1 *= gam;
    dh11 *= gam;
    dh12 *= gam;
    goto L140;
L160:
L170:
    if (! (abs(*dd2) <= rgamsq)) {
	goto L190;
    }
    if (*dd2 == zero) {
	goto L220;
    }
    igo = 2;
    igo_fmt = fmt_180;
/*              FIX-H.. */
    goto L70;
L180:
/* Computing 2nd power */
    d__1 = gam;
    *dd2 *= d__1 * d__1;
    dh21 /= gam;
    dh22 /= gam;
    goto L170;
L190:
L200:
    if (! (abs(*dd2) >= gamsq)) {
	goto L220;
    }
    igo = 3;
    igo_fmt = fmt_210;
/*              FIX-H.. */
    goto L70;
L210:
/* Computing 2nd power */
    d__1 = gam;
    *dd2 /= d__1 * d__1;
    dh21 *= gam;
    dh22 *= gam;
    goto L200;
L220:
    if (dflag < 0.) {
	goto L250;
    } else if (dflag == 0) {
	goto L230;
    } else {
	goto L240;
    }
L230:
    dparam[3] = dh21;
    dparam[4] = dh12;
    goto L260;
L240:
    dparam[2] = dh11;
    dparam[5] = dh22;
    goto L260;
L250:
    dparam[2] = dh11;
    dparam[3] = dh21;
    dparam[4] = dh12;
    dparam[5] = dh22;
L260:
    dparam[1] = dflag;
    return 0;
} /* drotmg_ */

