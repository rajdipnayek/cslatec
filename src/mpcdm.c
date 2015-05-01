/* mpcdm.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__1 = 1;
static integer c__4 = 4;
static integer c__0 = 0;

/* DECK MPCDM */
/* Subroutine */ int mpcdm_(doublereal *dx, integer *z__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, k, i2;
    static doublereal db;
    static integer ib;
    static doublereal dj;
    static integer ie, re, tp, rs;
    extern /* Subroutine */ int mpchk_(integer *, integer *), mpnzr_(integer *
	    , integer *, integer *, integer *), mpdivi_(integer *, integer *, 
	    integer *), mpmuli_(integer *, integer *, integer *);

/* ***BEGIN PROLOGUE  MPCDM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPCDM-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* Converts double-precision number DX to multiple-precision Z. */
/* Some numbers will not convert exactly on machines with base */
/* other than two, four or sixteen. This routine is not called */
/* by any other routine in 'mp', so may be omitted if double- */
/* precision is not available. */

/* The argument Z(*) and the variable R in COMMON are both INTEGER */
/* arrays of size 30.  See the comments in the routine MPBLAS for the */
/* for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPCHK, MPDIVI, MPMULI, MPNZR */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPCDM */
/* ***FIRST EXECUTABLE STATEMENT  MPCDM */
    /* Parameter adjustments */
    --z__;

    /* Function Body */
    mpchk_(&c__1, &c__4);
    i2 = mpcom_1.t + 4;
/* CHECK SIGN */
    if (*dx < 0.) {
	goto L20;
    } else if (*dx == 0) {
	goto L10;
    } else {
	goto L30;
    }
/* IF DX = 0D0 RETURN 0 */
L10:
    z__[1] = 0;
    return 0;
/* DX .LT. 0D0 */
L20:
    rs = -1;
    dj = -(*dx);
    goto L40;
/* DX .GT. 0D0 */
L30:
    rs = 1;
    dj = *dx;
L40:
    ie = 0;
L50:
    if (dj < 1.) {
	goto L60;
    }
/* INCREASE IE AND DIVIDE DJ BY 16. */
    ++ie;
    dj *= .0625;
    goto L50;
L60:
    if (dj >= .0625) {
	goto L70;
    }
    --ie;
    dj *= 16.;
    goto L60;
/* NOW DJ IS DY DIVIDED BY SUITABLE POWER OF 16 */
/* SET EXPONENT TO 0 */
L70:
    re = 0;
    db = (doublereal) mpcom_1.b;
/* CONVERSION LOOP (ASSUME DOUBLE-PRECISION OPS. EXACT) */
    i__1 = i2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dj = db * dj;
	mpcom_1.r__[i__ - 1] = (integer) dj;
/* L80: */
	dj -= (doublereal) mpcom_1.r__[i__ - 1];
    }
/* NORMALIZE RESULT */
    mpnzr_(&rs, &re, &z__[1], &c__0);
/* Computing MAX */
    i__1 = mpcom_1.b * 7 * mpcom_1.b;
    ib = max(i__1,32767) / 16;
    tp = 1;
/* NOW MULTIPLY BY 16**IE */
    if (ie < 0) {
	goto L90;
    } else if (ie == 0) {
	goto L130;
    } else {
	goto L110;
    }
L90:
    k = -ie;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tp <<= 4;
	if (tp <= ib && tp != mpcom_1.b && i__ < k) {
	    goto L100;
	}
	mpdivi_(&z__[1], &tp, &z__[1]);
	tp = 1;
L100:
	;
    }
    return 0;
L110:
    i__1 = ie;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tp <<= 4;
	if (tp <= ib && tp != mpcom_1.b && i__ < ie) {
	    goto L120;
	}
	mpmuli_(&z__[1], &tp, &z__[1]);
	tp = 1;
L120:
	;
    }
L130:
    return 0;
} /* mpcdm_ */

