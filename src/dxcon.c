/* dxcon.f -- translated by f2c (version 12.02.01).
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

static doublereal c_b11 = 10.;

/* DECK DXCON */
/* Subroutine */ int dxcon_(doublereal *x, integer *ix, integer *ierror)
{
    /* Initialized data */

    static integer ispace = 1;

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal a, b;
    static integer i__, j;
    static doublereal z__;
    static integer i1, j1, j2;
    extern /* Subroutine */ int dxc210_(integer *, doublereal *, integer *, 
	    integer *);
    static integer icase;
    extern /* Subroutine */ int dxadj_(doublereal *, integer *, integer *), 
	    dxred_(doublereal *, integer *, integer *);
    static integer itemp;

/* ***BEGIN PROLOGUE  DXCON */
/* ***PURPOSE  To provide double-precision floating-point arithmetic */
/*            with an extended exponent range. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  A3D */
/* ***TYPE      DOUBLE PRECISION (XCON-S, DXCON-D) */
/* ***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC */
/* ***AUTHOR  Lozier, Daniel W., (National Bureau of Standards) */
/*           Smith, John M., (NBS and George Mason University) */
/* ***DESCRIPTION */
/*     DOUBLE PRECISION X */
/*     INTEGER IX */

/*                  CONVERTS (X,IX) = X*RADIX**IX */
/*                  TO DECIMAL FORM IN PREPARATION FOR */
/*                  PRINTING, SO THAT (X,IX) = X*10**IX */
/*                  WHERE 1/10 .LE. ABS(X) .LT. 1 */
/*                  IS RETURNED, EXCEPT THAT IF */
/*                  (ABS(X),IX) IS BETWEEN RADIX**(-2L) */
/*                  AND RADIX**(2L) THEN THE REDUCED */
/*                  FORM WITH IX = 0 IS RETURNED. */

/* ***SEE ALSO  DXSET */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DXADJ, DXC210, DXRED */
/* ***COMMON BLOCKS    DXBLK2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820712  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXCON */

/*   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE */
/* ARE */
/*    (1) 4 .LE. L .LE. 2**NBITS - 1 - KMAX */

/*    (2) KMAX .LE. ((2**NBITS)-2)/LOG10R - L */

/* THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING */
/* IN SUBROUTINE DXSET. */



/*   THE PARAMETER ISPACE IS THE INCREMENT USED IN FORM- */
/* ING THE AUXILIARY INDEX OF THE DECIMAL EXTENDED-RANGE */
/* FORM. THE RETURNED VALUE OF IX WILL BE AN INTEGER MULT- */
/* IPLE OF ISPACE. ISPACE MUST SATISFY 1 .LE. ISPACE .LE. */
/* L/2. IF A VALUE GREATER THAN 1 IS TAKEN, THE RETURNED */
/* VALUE OF X WILL SATISFY 10**(-ISPACE) .LE. ABS(X) .LE. 1 */
/* WHEN (ABS(X),IX) .LT. RADIX**(-2L), AND 1/10 .LE. ABS(X) */
/* .LT. 10**(ISPACE-1) WHEN (ABS(X),IX) .GT. RADIX**(2L). */

/* ***FIRST EXECUTABLE STATEMENT  DXCON */
    *ierror = 0;
    dxred_(x, ix, ierror);
    if (*ierror != 0) {
	return 0;
    }
    if (*ix == 0) {
	goto L150;
    }
    dxadj_(x, ix, ierror);
    if (*ierror != 0) {
	return 0;
    }

/* CASE 1 IS WHEN (X,IX) IS LESS THAN RADIX**(-2L) IN MAGNITUDE, */
/* CASE 2 IS WHEN (X,IX) IS GREATER THAN RADIX**(2L) IN MAGNITUDE. */
    itemp = 1;
    icase = (i_sign(&itemp, ix) + 3) / 2;
    switch (icase) {
	case 1:  goto L10;
	case 2:  goto L20;
    }
L10:
    if (abs(*x) < 1.) {
	goto L30;
    }
    *x /= dxblk2_1.radixl;
    *ix += dxblk2_1.l;
    goto L30;
L20:
    if (abs(*x) >= 1.) {
	goto L30;
    }
    *x *= dxblk2_1.radixl;
    *ix -= dxblk2_1.l;
L30:

/* AT THIS POINT, RADIX**(-L) .LE. ABS(X) .LT. 1.0D0     IN CASE 1, */
/*                      1.0D0 .LE. ABS(X) .LT. RADIX**L  IN CASE 2. */
    d__1 = abs(*x);
    i__ = (integer) (d_lg10(&d__1) / dxblk2_1.dlg10r);
    a = pow_di(&dxblk2_1.radix, &i__);
    switch (icase) {
	case 1:  goto L40;
	case 2:  goto L60;
    }
L40:
    if (a <= dxblk2_1.radix * abs(*x)) {
	goto L50;
    }
    --i__;
    a /= dxblk2_1.radix;
    goto L40;
L50:
    if (abs(*x) < a) {
	goto L80;
    }
    ++i__;
    a *= dxblk2_1.radix;
    goto L50;
L60:
    if (a <= abs(*x)) {
	goto L70;
    }
    --i__;
    a /= dxblk2_1.radix;
    goto L60;
L70:
    if (abs(*x) < dxblk2_1.radix * a) {
	goto L80;
    }
    ++i__;
    a *= dxblk2_1.radix;
    goto L70;
L80:

/* AT THIS POINT I IS SUCH THAT */
/* RADIX**(I-1) .LE. ABS(X) .LT. RADIX**I      IN CASE 1, */
/*     RADIX**I .LE. ABS(X) .LT. RADIX**(I+1)  IN CASE 2. */
    itemp = (integer) (ispace / dxblk2_1.dlg10r);
    a = pow_di(&dxblk2_1.radix, &itemp);
    b = pow_di(&c_b11, &ispace);
L90:
    if (a <= b) {
	goto L100;
    }
    --itemp;
    a /= dxblk2_1.radix;
    goto L90;
L100:
    if (b < a * dxblk2_1.radix) {
	goto L110;
    }
    ++itemp;
    a *= dxblk2_1.radix;
    goto L100;
L110:

/* AT THIS POINT ITEMP IS SUCH THAT */
/* RADIX**ITEMP .LE. 10**ISPACE .LT. RADIX**(ITEMP+1). */
    if (itemp > 0) {
	goto L120;
    }
/* ITEMP = 0 IF, AND ONLY IF, ISPACE = 1 AND RADIX = 16.0D0 */
    i__1 = -i__;
    *x *= pow_di(&dxblk2_1.radix, &i__1);
    *ix += i__;
    dxc210_(ix, &z__, &j, ierror);
    if (*ierror != 0) {
	return 0;
    }
    *x *= z__;
    *ix = j;
    switch (icase) {
	case 1:  goto L130;
	case 2:  goto L140;
    }
L120:
    i1 = i__ / itemp;
    i__1 = -i1 * itemp;
    *x *= pow_di(&dxblk2_1.radix, &i__1);
    *ix += i1 * itemp;

/* AT THIS POINT, */
/* RADIX**(-ITEMP) .LE. ABS(X) .LT. 1.0D0        IN CASE 1, */
/*           1.0D0 .LE. ABS(X) .LT. RADIX**ITEMP IN CASE 2. */
    dxc210_(ix, &z__, &j, ierror);
    if (*ierror != 0) {
	return 0;
    }
    j1 = j / ispace;
    j2 = j - j1 * ispace;
    *x = *x * z__ * pow_di(&c_b11, &j2);
    *ix = j1 * ispace;

/* AT THIS POINT, */
/*  10.0D0**(-2*ISPACE) .LE. ABS(X) .LT. 1.0D0                IN CASE 1, */
/*           10.0D0**-1 .LE. ABS(X) .LT. 10.0D0**(2*ISPACE-1) IN CASE 2. */
    switch (icase) {
	case 1:  goto L130;
	case 2:  goto L140;
    }
L130:
    if (b * abs(*x) >= 1.) {
	goto L150;
    }
    *x *= b;
    *ix -= ispace;
    goto L130;
L140:
    if (abs(*x) * 10. < b) {
	goto L150;
    }
    *x /= b;
    *ix += ispace;
    goto L140;
L150:
    return 0;
} /* dxcon_ */

