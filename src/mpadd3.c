/* mpadd3.f -- translated by f2c (version 12.02.01).
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
    integer b, t, m, lun, mxr, r__[30];
} mpcom_;

#define mpcom_1 mpcom_

/* DECK MPADD3 */
/* Subroutine */ int mpadd3_(integer *x, integer *y, integer *s, integer *med,
	 integer *re)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer c__, i__, j, i2, i2p, ted;

/* ***BEGIN PROLOGUE  MPADD3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPADD3-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   Called by MPADD2; does inner loops of addition */

/*   The arguments X(*) and Y(*) and the variable R in COMMON are all */
/*   INTEGER arrays of size 30.  See the comments in the routine MPBLAS */
/*   for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPADD3 */
/* ***FIRST EXECUTABLE STATEMENT  MPADD3 */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    ted = mpcom_1.t + *med;
    i2 = mpcom_1.t + 4;
    i__ = i2;
    c__ = 0;
/* CLEAR GUARD DIGITS TO RIGHT OF X DIGITS */
L10:
    if (i__ <= ted) {
	goto L20;
    }
    mpcom_1.r__[i__ - 1] = 0;
    --i__;
    goto L10;
L20:
    if (*s < 0) {
	goto L130;
    }
/* HERE DO ADDITION, EXPONENT(Y) .GE. EXPONENT(X) */
    if (i__ < mpcom_1.t) {
	goto L40;
    }
L30:
    j = i__ - *med;
    mpcom_1.r__[i__ - 1] = x[j + 2];
    --i__;
    if (i__ > mpcom_1.t) {
	goto L30;
    }
L40:
    if (i__ <= *med) {
	goto L60;
    }
    j = i__ - *med;
    c__ = y[i__ + 2] + x[j + 2] + c__;
    if (c__ < mpcom_1.b) {
	goto L50;
    }
/* CARRY GENERATED HERE */
    mpcom_1.r__[i__ - 1] = c__ - mpcom_1.b;
    c__ = 1;
    --i__;
    goto L40;
/* NO CARRY GENERATED HERE */
L50:
    mpcom_1.r__[i__ - 1] = c__;
    c__ = 0;
    --i__;
    goto L40;
L60:
    if (i__ <= 0) {
	goto L90;
    }
    c__ = y[i__ + 2] + c__;
    if (c__ < mpcom_1.b) {
	goto L70;
    }
    mpcom_1.r__[i__ - 1] = 0;
    c__ = 1;
    --i__;
    goto L60;
L70:
    mpcom_1.r__[i__ - 1] = c__;
    --i__;
/* NO CARRY POSSIBLE HERE */
L80:
    if (i__ <= 0) {
	return 0;
    }
    mpcom_1.r__[i__ - 1] = y[i__ + 2];
    --i__;
    goto L80;
L90:
    if (c__ == 0) {
	return 0;
    }
/* MUST SHIFT RIGHT HERE AS CARRY OFF END */
    i2p = i2 + 1;
    i__1 = i2;
    for (j = 2; j <= i__1; ++j) {
	i__ = i2p - j;
/* L100: */
	mpcom_1.r__[i__] = mpcom_1.r__[i__ - 1];
    }
    mpcom_1.r__[0] = 1;
    ++(*re);
    return 0;
/* HERE DO SUBTRACTION, ABS(Y) .GT. ABS(X) */
L110:
    j = i__ - *med;
    mpcom_1.r__[i__ - 1] = c__ - x[j + 2];
    c__ = 0;
    if (mpcom_1.r__[i__ - 1] >= 0) {
	goto L120;
    }
/* BORROW GENERATED HERE */
    c__ = -1;
    mpcom_1.r__[i__ - 1] += mpcom_1.b;
L120:
    --i__;
L130:
    if (i__ > mpcom_1.t) {
	goto L110;
    }
L140:
    if (i__ <= *med) {
	goto L160;
    }
    j = i__ - *med;
    c__ = y[i__ + 2] + c__ - x[j + 2];
    if (c__ >= 0) {
	goto L150;
    }
/* BORROW GENERATED HERE */
    mpcom_1.r__[i__ - 1] = c__ + mpcom_1.b;
    c__ = -1;
    --i__;
    goto L140;
/* NO BORROW GENERATED HERE */
L150:
    mpcom_1.r__[i__ - 1] = c__;
    c__ = 0;
    --i__;
    goto L140;
L160:
    if (i__ <= 0) {
	return 0;
    }
    c__ = y[i__ + 2] + c__;
    if (c__ >= 0) {
	goto L70;
    }
    mpcom_1.r__[i__ - 1] = c__ + mpcom_1.b;
    c__ = -1;
    --i__;
    goto L160;
} /* mpadd3_ */

