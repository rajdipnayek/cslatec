/* xadd.f -- translated by f2c (version 12.02.01).
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

/* DECK XADD */
/* Subroutine */ int xadd_(real *x, integer *ix, real *y, integer *iy, real *
	z__, integer *iz, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static real s, t;
    static integer i1, i2, is;
    extern /* Subroutine */ int xadj_(real *, integer *, integer *);

/* ***BEGIN PROLOGUE  XADD */
/* ***PURPOSE  To provide single-precision floating-point arithmetic */
/*            with an extended exponent range. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  A3D */
/* ***TYPE      SINGLE PRECISION (XADD-S, DXADD-D) */
/* ***KEYWORDS  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */
/*     REAL X, Y, Z */
/*     INTEGER IX, IY, IZ */

/*                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) = */
/*                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED */
/*                  BEFORE RETURNING. THE INPUT OPERANDS */
/*                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR */
/*                  PRINCIPAL PARTS MUST SATISFY */
/*                  RADIX**(-2L).LE.ABS(X).LE.RADIX**(2L), */
/*                  RADIX**(-2L).LE.ABS(Y).LE.RADIX**(2L). */

/* ***SEE ALSO  XSET */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  XADJ */
/* ***COMMON BLOCKS    XBLK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  XADD */


/*   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE */
/* ARE */
/*     (1) 1 .LT. L .LE. 0.5*LOGR(0.5*DZERO) */

/*     (2) NRADPL .LT. L .LE. KMAX/6 */

/*     (3) KMAX .LE. (2**NBITS - 4*L - 1)/2 */

/* THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING */
/* IN SUBROUTINE XSET. */

/* ***FIRST EXECUTABLE STATEMENT  XADD */
    *ierror = 0;
    if (*x != 0.f) {
	goto L10;
    }
    *z__ = *y;
    *iz = *iy;
    goto L220;
L10:
    if (*y != 0.f) {
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
    if (abs(*ix) <= xblk2_1.l * 6 && abs(*iy) <= xblk2_1.l * 6) {
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
    if (dabs(*x) > 1.f && dabs(*y) > 1.f) {
	goto L60;
    }
    if (dabs(*x) < 1.f && dabs(*y) < 1.f) {
	goto L70;
    }
    *z__ = *x + *y;
    *iz = *ix;
    goto L220;
L60:
    s = *x / xblk2_1.radixl;
    t = *y / xblk2_1.radixl;
    *z__ = s + t;
    *iz = *ix + xblk2_1.l;
    goto L220;
L70:
    s = *x * xblk2_1.radixl;
    t = *y * xblk2_1.radixl;
    *z__ = s + t;
    *iz = *ix - xblk2_1.l;
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

    i1 = abs(i__) / xblk2_1.l;
    i2 = abs(i__) % xblk2_1.l;
    if (dabs(t) >= xblk2_1.radixl) {
	goto L130;
    }
    if (dabs(t) >= 1.f) {
	goto L120;
    }
    if (xblk2_1.radixl * dabs(t) >= 1.f) {
	goto L110;
    }
    j = i1 + 1;
    i__1 = xblk2_1.l - i2;
    t *= pow_ri(&xblk2_1.radix, &i__1);
    goto L140;
L110:
    j = i1;
    i__1 = -i2;
    t *= pow_ri(&xblk2_1.radix, &i__1);
    goto L140;
L120:
    j = i1 - 1;
    if (j < 0) {
	goto L110;
    }
    i__1 = -i2;
    t = t * pow_ri(&xblk2_1.radix, &i__1) / xblk2_1.radixl;
    goto L140;
L130:
    j = i1 - 2;
    if (j < 0) {
	goto L120;
    }
    i__1 = -i2;
    t = t * pow_ri(&xblk2_1.radix, &i__1) / xblk2_1.rad2l;
L140:

/*  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE */
/* AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT */
/* OF T.  THE SHIFTED VALUE OF T SATISFIES */

/*       RADIX**(-2*L) .LE. ABS(T) .LE. 1.0 */

/* AND, IF J=0, NO FURTHER SHIFTING REMAINS TO BE DONE. */

    if (j == 0) {
	goto L190;
    }
    if (dabs(s) >= xblk2_1.radixl || j > 3) {
	goto L150;
    }
    if (dabs(s) >= 1.f) {
	switch (j) {
	    case 1:  goto L180;
	    case 2:  goto L150;
	    case 3:  goto L150;
	}
    }
    if (xblk2_1.radixl * dabs(s) >= 1.f) {
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
    s *= xblk2_1.radixl;
L170:
    s *= xblk2_1.radixl;
L180:
    s *= xblk2_1.radixl;
L190:

/*   AT THIS POINT, THE REMAINING DIFFERENCE IN THE */
/* AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT */
/* OF S.  IF THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED */
/* RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE */
/* SUM. */

    if (dabs(s) > 1.f && dabs(t) > 1.f) {
	goto L200;
    }
    if (dabs(s) < 1.f && dabs(t) < 1.f) {
	goto L210;
    }
    *z__ = s + t;
    *iz = is - j * xblk2_1.l;
    goto L220;
L200:
    s /= xblk2_1.radixl;
    t /= xblk2_1.radixl;
    *z__ = s + t;
    *iz = is - j * xblk2_1.l + xblk2_1.l;
    goto L220;
L210:
    s *= xblk2_1.radixl;
    t *= xblk2_1.radixl;
    *z__ = s + t;
    *iz = is - j * xblk2_1.l - xblk2_1.l;
L220:
    xadj_(z__, iz, ierror);
    return 0;
} /* xadd_ */

