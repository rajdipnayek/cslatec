/* xred.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    real radix, radixl, rad2l, dlg10r;
    integer l, l2, kmax;
} xblk2_;

#define xblk2_1 xblk2_

/* DECK XRED */
/* Subroutine */ int xred_(real *x, integer *ix, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;
    static real xa;
    static integer ixa, ixa1, ixa2;

/* ***BEGIN PROLOGUE  XRED */
/* ***PURPOSE  To provide single-precision floating-point arithmetic */
/*            with an extended exponent range. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  A3D */
/* ***TYPE      SINGLE PRECISION (XRED-S, DXRED-D) */
/* ***KEYWORDS  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */
/*     REAL X */
/*     INTEGER IX */

/*                  IF */
/*                  RADIX**(-2L) .LE. (ABS(X),IX) .LE. RADIX**(2L) */
/*                  THEN XRED TRANSFORMS (X,IX) SO THAT IX=0. */
/*                  IF (X,IX) IS OUTSIDE THE ABOVE RANGE, */
/*                  THEN XRED TAKES NO ACTION. */
/*                  THIS SUBROUTINE IS USEFUL IF THE */
/*                  RESULTS OF EXTENDED-RANGE CALCULATIONS */
/*                  ARE TO BE USED IN SUBSEQUENT ORDINARY */
/*                  SINGLE-PRECISION CALCULATIONS. */

/* ***SEE ALSO  XSET */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    XBLK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  XRED */

/* ***FIRST EXECUTABLE STATEMENT  XRED */
    *ierror = 0;
    if (*x == 0.f) {
	goto L90;
    }
    xa = dabs(*x);
    if (*ix == 0) {
	goto L70;
    }
    ixa = abs(*ix);
    ixa1 = ixa / xblk2_1.l2;
    ixa2 = ixa % xblk2_1.l2;
    if (*ix > 0) {
	goto L40;
    }
L10:
    if (xa > 1.f) {
	goto L20;
    }
    xa *= xblk2_1.rad2l;
    ++ixa1;
    goto L10;
L20:
    xa /= pow_ri(&xblk2_1.radix, &ixa2);
    if (ixa1 == 0) {
	goto L70;
    }
    i__1 = ixa1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (xa < 1.f) {
	    goto L100;
	}
	xa /= xblk2_1.rad2l;
/* L30: */
    }
    goto L70;

L40:
    if (xa < 1.f) {
	goto L50;
    }
    xa /= xblk2_1.rad2l;
    ++ixa1;
    goto L40;
L50:
    xa *= pow_ri(&xblk2_1.radix, &ixa2);
    if (ixa1 == 0) {
	goto L70;
    }
    i__1 = ixa1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (xa > 1.f) {
	    goto L100;
	}
	xa *= xblk2_1.rad2l;
/* L60: */
    }
L70:
    if (xa > xblk2_1.rad2l) {
	goto L100;
    }
    if (xa > 1.f) {
	goto L80;
    }
    if (xblk2_1.rad2l * xa < 1.f) {
	goto L100;
    }
L80:
    *x = r_sign(&xa, x);
L90:
    *ix = 0;
L100:
    return 0;
} /* xred_ */

