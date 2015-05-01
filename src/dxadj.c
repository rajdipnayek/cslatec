/* dxadj.f -- translated by f2c (version 12.02.01).
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
    doublereal radix, radixl, rad2l, dlg10r;
    integer l, l2, kmax;
} dxblk2_;

#define dxblk2_1 dxblk2_

/* Table of constant values */

static integer c__207 = 207;
static integer c__1 = 1;

/* DECK DXADJ */
/* Subroutine */ int dxadj_(doublereal *x, integer *ix, integer *ierror)
{
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DXADJ */
/* ***PURPOSE  To provide double-precision floating-point arithmetic */
/*            with an extended exponent range. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  A3D */
/* ***TYPE      DOUBLE PRECISION (XADJ-S, DXADJ-D) */
/* ***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */
/*     DOUBLE PRECISION X */
/*     INTEGER IX */

/*                  TRANSFORMS (X,IX) SO THAT */
/*                  RADIX**(-L) .LE. ABS(X) .LT. RADIX**L. */
/*                  ON MOST COMPUTERS THIS TRANSFORMATION DOES */
/*                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS */
/*                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC. */

/* ***SEE ALSO  DXSET */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  XERMSG */
/* ***COMMON BLOCKS    DXBLK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*           CALLs to XERROR changed to CALLs to XERMSG.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXADJ */

/*   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE */
/* IS */
/*     2*L .LE. KMAX */

/* THIS CONDITION MUST BE MET BY APPROPRIATE CODING */
/* IN SUBROUTINE DXSET. */

/* ***FIRST EXECUTABLE STATEMENT  DXADJ */
    *ierror = 0;
    if (*x == 0.) {
	goto L50;
    }
    if (abs(*x) >= 1.) {
	goto L20;
    }
    if (dxblk2_1.radixl * abs(*x) >= 1.) {
	goto L60;
    }
    *x *= dxblk2_1.rad2l;
    if (*ix < 0) {
	goto L10;
    }
    *ix -= dxblk2_1.l2;
    goto L70;
L10:
    if (*ix < -dxblk2_1.kmax + dxblk2_1.l2) {
	goto L40;
    }
    *ix -= dxblk2_1.l2;
    goto L70;
L20:
    if (abs(*x) < dxblk2_1.radixl) {
	goto L60;
    }
    *x /= dxblk2_1.rad2l;
    if (*ix > 0) {
	goto L30;
    }
    *ix += dxblk2_1.l2;
    goto L70;
L30:
    if (*ix > dxblk2_1.kmax - dxblk2_1.l2) {
	goto L40;
    }
    *ix += dxblk2_1.l2;
    goto L70;
L40:
    xermsg_("SLATEC", "DXADJ", "overflow in auxiliary index", &c__207, &c__1, 
	    (ftnlen)6, (ftnlen)5, (ftnlen)27);
    *ierror = 207;
    return 0;
L50:
    *ix = 0;
L60:
    if (abs(*ix) > dxblk2_1.kmax) {
	goto L40;
    }
L70:
    return 0;
} /* dxadj_ */

