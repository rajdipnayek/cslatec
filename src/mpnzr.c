/* mpnzr.f -- translated by f2c (version 12.02.01).
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

/* DECK MPNZR */
/* Subroutine */ int mpnzr_(integer *rs, integer *re, integer *z__, integer *
	trunc)
{
    /* Format strings */
    static char fmt_30[] = "(\002 *** SIGN NOT 0, +1 OR -1 IN CALL TO MPNZR"
	    ",\002,\002 POSSIBLE OVERWRITING PROBLEM ***\002)";
    static char fmt_160[] = "(\002 *** OVERFLOW OCCURRED IN MPNZR ***\002)";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, b2, i2, is, it, i2m, i2p;
    extern /* Subroutine */ int mperr_(void), mpunfl_(integer *), mpovfl_(
	    integer *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_160, 0 };


/* ***BEGIN PROLOGUE  MPNZR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPNZR-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Modified for use with BLAS.  Blank COMMON changed to named COMMON. */
/*  Assumes long (i.e. (t+4)-DIGIT) fraction in R, sign = RS, exponent */
/*  = RE.  Normalizes, and returns 'mp' result in Z. Integer arguments */
/*  RS and RE are not preserved. R*-rounding is used if TRUNC.EQ.0 */

/*  The argument Z(*) and the variable R in COMMON are INTEGER arrays */
/*  of size 30.  See the comments in the routine MPBLAS for the reason */
/*  for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPERR, MPOVFL, MPUNFL */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPNZR */
/* ***FIRST EXECUTABLE STATEMENT  MPNZR */
    /* Parameter adjustments */
    --z__;

    /* Function Body */
    i2 = mpcom_1.t + 4;
    if (*rs != 0) {
	goto L20;
    }
/* STORE ZERO IN Z */
L10:
    z__[1] = 0;
    return 0;
/* CHECK THAT SIGN = +-1 */
L20:
    if (abs(*rs) <= 1) {
	goto L40;
    }
    io___2.ciunit = mpcom_1.lun;
    s_wsfe(&io___2);
    e_wsfe();
    mperr_();
    goto L10;
/* LOOK FOR FIRST NONZERO DIGIT */
L40:
    i__1 = i2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	is = i__ - 1;
	if (mpcom_1.r__[i__ - 1] > 0) {
	    goto L60;
	}
/* L50: */
    }
/* FRACTION ZERO */
    goto L10;
L60:
    if (is == 0) {
	goto L90;
    }
/* NORMALIZE */
    *re -= is;
    i2m = i2 - is;
    i__1 = i2m;
    for (j = 1; j <= i__1; ++j) {
	k = j + is;
/* L70: */
	mpcom_1.r__[j - 1] = mpcom_1.r__[k - 1];
    }
    i2p = i2m + 1;
    i__1 = i2;
    for (j = i2p; j <= i__1; ++j) {
/* L80: */
	mpcom_1.r__[j - 1] = 0;
    }
/* CHECK TO SEE IF TRUNCATION IS DESIRED */
L90:
    if (*trunc != 0) {
	goto L150;
    }
/* SEE IF ROUNDING NECESSARY */
/* TREAT EVEN AND ODD BASES DIFFERENTLY */
    b2 = mpcom_1.b / 2;
    if (b2 << 1 != mpcom_1.b) {
	goto L130;
    }
/* B EVEN.  ROUND IF R(T+1).GE.B2 UNLESS R(T) ODD AND ALL ZEROS */
/* AFTER R(T+2). */
    if ((i__1 = mpcom_1.r__[mpcom_1.t] - b2) < 0) {
	goto L150;
    } else if (i__1 == 0) {
	goto L100;
    } else {
	goto L110;
    }
L100:
    if (mpcom_1.r__[mpcom_1.t - 1] % 2 == 0) {
	goto L110;
    }
    if (mpcom_1.r__[mpcom_1.t + 1] + mpcom_1.r__[mpcom_1.t + 2] + mpcom_1.r__[
	    mpcom_1.t + 3] == 0) {
	goto L150;
    }
/* ROUND */
L110:
    i__1 = mpcom_1.t;
    for (j = 1; j <= i__1; ++j) {
	i__ = mpcom_1.t + 1 - j;
	++mpcom_1.r__[i__ - 1];
	if (mpcom_1.r__[i__ - 1] < mpcom_1.b) {
	    goto L150;
	}
/* L120: */
	mpcom_1.r__[i__ - 1] = 0;
    }
/* EXCEPTIONAL CASE, ROUNDED UP TO .10000... */
    ++(*re);
    mpcom_1.r__[0] = 1;
    goto L150;
/* ODD BASE, ROUND IF R(T+1)... .GT. 1/2 */
L130:
    for (i__ = 1; i__ <= 4; ++i__) {
	it = mpcom_1.t + i__;
	if ((i__1 = mpcom_1.r__[it - 1] - b2) < 0) {
	    goto L150;
	} else if (i__1 == 0) {
	    goto L140;
	} else {
	    goto L110;
	}
L140:
	;
    }
/* CHECK FOR OVERFLOW */
L150:
    if (*re <= mpcom_1.m) {
	goto L170;
    }
    io___11.ciunit = mpcom_1.lun;
    s_wsfe(&io___11);
    e_wsfe();
    mpovfl_(&z__[1]);
    return 0;
/* CHECK FOR UNDERFLOW */
L170:
    if (*re < -mpcom_1.m) {
	goto L190;
    }
/* STORE RESULT IN Z */
    z__[1] = *rs;
    z__[2] = *re;
    i__1 = mpcom_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L180: */
	z__[i__ + 2] = mpcom_1.r__[i__ - 1];
    }
    return 0;
/* UNDERFLOW HERE */
L190:
    mpunfl_(&z__[1]);
    return 0;
} /* mpnzr_ */

