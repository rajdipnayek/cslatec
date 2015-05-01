/* dxadd.f -- translated by f2c (version 12.02.01).
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

/* DECK DXADD */
/* Subroutine */ int dxadd_(doublereal *x, integer *ix, doublereal *y, 
	integer *iy, doublereal *z__, integer *iz, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static doublereal s, t;
    static integer i1, i2, is;
    extern /* Subroutine */ int dxadj_(doublereal *, integer *, integer *);

/* ***BEGIN PROLOGUE  DXADD */
/* ***PURPOSE  To provide double-precision floating-point arithmetic */
/*            with an extended exponent range. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  A3D */
/* ***TYPE      DOUBLE PRECISION (XADD-S, DXADD-D) */
/* ***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */
/*     DOUBLE PRECISION X, Y, Z */
/*     INTEGER IX, IY, IZ */

/*                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) = */
/*                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED */
/*                  BEFORE RETURNING. THE INPUT OPERANDS */
/*                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR */
/*                  PRINCIPAL PARTS MUST SATISFY */
/*                  RADIX**(-2L).LE.ABS(X).LE.RADIX**(2L), */
/*                  RADIX**(-2L).LE.ABS(Y).LE.RADIX**(2L). */

/* ***SEE ALSO  DXSET */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DXADJ */
/* ***COMMON BLOCKS    DXBLK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXADD */

/*   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE */
/* ARE */
/*     (1) 1 .LT. L .LE. 0.5D0*LOGR(0.5D0*DZERO) */

/*     (2) NRADPL .LT. L .LE. KMAX/6 */

/*     (3) KMAX .LE. (2**NBITS - 4*L - 1)/2 */

/* THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING */
/* IN SUBROUTINE DXSET. */

/* ***FIRST EXECUTABLE STATEMENT  DXADD */
    *ierror = 0;
    if (*x != 0.) {
	goto L10;
    }
    *z__ = *y;
    *iz = *iy;
    goto L220;
L10:
    if (*y != 0.) {
	goto L20;
    }
    *z__ = *x;
    *iz = *ix;
    goto L220;
L20:
    if (*ix >= 0 && *iy >= 0) {
	goto L40;
    }
    if (*ix < 0 && *iy < 0) {
	goto L40;
    }
    if (abs(*ix) <= dxblk2_1.l * 6 && abs(*iy) <= dxblk2_1.l * 6) {
	goto L40;
    }
    if (*ix >= 0) {
	goto L30;
    }
    *z__ = *y;
    *iz = *iy;
    goto L220;
L30:
    *z__ = *x;
    *iz = *ix;
    goto L220;
L40:
    i__ = *ix - *iy;
    if (i__ < 0) {
	goto L80;
    } else if (i__ == 0) {
	goto L50;
    } else {
	goto L90;
    }
L50:
    if (abs(*x) > 1. && abs(*y) > 1.) {
	goto L60;
    }
    if (abs(*x) < 1. && abs(*y) < 1.) {
	goto L70;
    }
    *z__ = *x + *y;
    *iz = *ix;
    goto L220;
L60:
    s = *x / dxblk2_1.radixl;
    t = *y / dxblk2_1.radixl;
    *z__ = s + t;
    *iz = *ix + dxblk2_1.l;
    goto L220;
L70:
    s = *x * dxblk2_1.radixl;
    t = *y * dxblk2_1.radixl;
    *z__ = s + t;
    *iz = *ix - dxblk2_1.l;
    goto L220;
L80:
    s = *y;
    is = *iy;
    t = *x;
    goto L100;
L90:
    s = *x;
    is = *ix;
    t = *y;
L100:

/*  AT THIS POINT, THE ONE OF (X,IX) OR (Y,IY) THAT HAS THE */
/* LARGER AUXILIARY INDEX IS STORED IN (S,IS). THE PRINCIPAL */
/* PART OF THE OTHER INPUT IS STORED IN T. */

    i1 = abs(i__) / dxblk2_1.l;
    i2 = abs(i__) % dxblk2_1.l;
    if (abs(t) >= dxblk2_1.radixl) {
	goto L130;
    }
    if (abs(t) >= 1.) {
	goto L120;
    }
    if (dxblk2_1.radixl * abs(t) >= 1.) {
	goto L110;
    }
    j = i1 + 1;
    i__1 = dxblk2_1.l - i2;
    t *= pow_di(&dxblk2_1.radix, &i__1);
    goto L140;
L110:
    j = i1;
    i__1 = -i2;
    t *= pow_di(&dxblk2_1.radix, &i__1);
    goto L140;
L120:
    j = i1 - 1;
    if (j < 0) {
	goto L110;
    }
    i__1 = -i2;
    t = t * pow_di(&dxblk2_1.radix, &i__1) / dxblk2_1.radixl;
    goto L140;
L130:
    j = i1 - 2;
    if (j < 0) {
	goto L120;
    }
    i__1 = -i2;
    t = t * pow_di(&dxblk2_1.radix, &i__1) / dxblk2_1.rad2l;
L140:

/*  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE */
/* AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT */
/* OF T.  THE SHIFTED VALUE OF T SATISFIES */

/*       RADIX**(-2*L) .LE. ABS(T) .LE. 1.0D0 */

/* AND, IF J=0, NO FURTHER SHIFTING REMAINS TO BE DONE. */

    if (j == 0) {
	goto L190;
    }
    if (abs(s) >= dxblk2_1.radixl || j > 3) {
	goto L150;
    }
    if (abs(s) >= 1.) {
	switch (j) {
	    case 1:  goto L180;
	    case 2:  goto L150;
	    case 3:  goto L150;
	}
    }
    if (dxblk2_1.radixl * abs(s) >= 1.) {
	switch (j) {
	    case 1:  goto L180;
	    case 2:  goto L170;
	    case 3:  goto L150;
	}
    }
    switch (j) {
	case 1:  goto L180;
	case 2:  goto L170;
	case 3:  goto L160;
    }
L150:
    *z__ = s;
    *iz = is;
    goto L220;
L160:
    s *= dxblk2_1.radixl;
L170:
    s *= dxblk2_1.radixl;
L180:
    s *= dxblk2_1.radixl;
L190:

/*   AT THIS POINT, THE REMAINING DIFFERENCE IN THE */
/* AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT */
/* OF S.  IF THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED */
/* RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE */
/* SUM. */

    if (abs(s) > 1. && abs(t) > 1.) {
	goto L200;
    }
    if (abs(s) < 1. && abs(t) < 1.) {
	goto L210;
    }
    *z__ = s + t;
    *iz = is - j * dxblk2_1.l;
    goto L220;
L200:
    s /= dxblk2_1.radixl;
    t /= dxblk2_1.radixl;
    *z__ = s + t;
    *iz = is - j * dxblk2_1.l + dxblk2_1.l;
    goto L220;
L210:
    s *= dxblk2_1.radixl;
    t *= dxblk2_1.radixl;
    *z__ = s + t;
    *iz = is - j * dxblk2_1.l - dxblk2_1.l;
L220:
    dxadj_(z__, iz, ierror);
    return 0;
} /* dxadd_ */

