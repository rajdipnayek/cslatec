/* mpmul.f -- translated by f2c (version 12.02.01).
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

/* DECK MPMUL */
/* Subroutine */ int mpmul_(integer *x, integer *y, integer *z__)
{
    /* Format strings */
    static char fmt_80[] = "(\002 *** INTEGER OVERFLOW IN MPMUL, B TOO LARGE"
	    " ***\002)";
    static char fmt_100[] = "(\002 *** ILLEGAL BASE B DIGIT IN CALL TO MPM"
	    "UL,\002,\002 POSSIBLE OVERWRITING PROBLEM ***\002)";

    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer c__, i__, j, i2, j1, re, ri, xi, rs, i2p;
    extern /* Subroutine */ int mpchk_(integer *, integer *), mpmlp_(integer *
	    , integer *, integer *, integer *), mperr_(void), mpnzr_(integer *
	    , integer *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 0, 0, fmt_80, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_100, 0 };


/* ***BEGIN PROLOGUE  MPMUL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPMUL-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Multiplies X and Y, returning result in Z, for 'mp' X, Y and Z. */
/*  The simple o(t**2) algorithm is used, with four guard digits and */
/*  R*-rounding. Advantage is taken of zero digits in X, but not in Y. */
/*  Asymptotically faster algorithms are known (see Knuth, VOL. 2), */
/*  but are difficult to implement in FORTRAN in an efficient and */
/*  machine-independent manner. In comments to other 'mp' routines, */
/*  M(t) is the time to perform t-digit 'mp' multiplication. Thus */
/*  M(t) = o(t**2) with the present version of MPMUL, but */
/*  M(t) = o(t.log(t).log(log(t))) is theoretically possible. */

/*  The arguments X(*), Y(*), and Z(*), and the variable R in COMMON are */
/*  all INTEGER arrays of size 30.  See the comments in the routine */
/*  MPBLAS for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPCHK, MPERR, MPMLP, MPNZR */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPMUL */
/* ***FIRST EXECUTABLE STATEMENT  MPMUL */
    /* Parameter adjustments */
    --z__;
    --y;
    --x;

    /* Function Body */
    mpchk_(&c__1, &c__4);
    i2 = mpcom_1.t + 4;
    i2p = i2 + 1;
/* FORM SIGN OF PRODUCT */
    rs = x[1] * y[1];
    if (rs != 0) {
	goto L10;
    }
/* SET RESULT TO ZERO */
    z__[1] = 0;
    return 0;
/* FORM EXPONENT OF PRODUCT */
L10:
    re = x[2] + y[2];
/* CLEAR ACCUMULATOR */
    i__1 = i2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	mpcom_1.r__[i__ - 1] = 0;
    }
/* PERFORM MULTIPLICATION */
    c__ = 8;
    i__1 = mpcom_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = x[i__ + 2];
/* FOR SPEED, PUT THE NUMBER WITH MANY ZEROS FIRST */
	if (xi == 0) {
	    goto L40;
	}
/* Computing MIN */
	i__3 = mpcom_1.t, i__4 = i2 - i__;
	i__2 = min(i__3,i__4);
	mpmlp_(&mpcom_1.r__[i__], &y[3], &xi, &i__2);
	--c__;
	if (c__ > 0) {
	    goto L40;
	}
/* CHECK FOR LEGAL BASE B DIGIT */
	if (xi < 0 || xi >= mpcom_1.b) {
	    goto L90;
	}
/* PROPAGATE CARRIES AT END AND EVERY EIGHTH TIME, */
/* FASTER THAN DOING IT EVERY TIME. */
	i__2 = i2;
	for (j = 1; j <= i__2; ++j) {
	    j1 = i2p - j;
	    ri = mpcom_1.r__[j1 - 1] + c__;
	    if (ri < 0) {
		goto L70;
	    }
	    c__ = ri / mpcom_1.b;
/* L30: */
	    mpcom_1.r__[j1 - 1] = ri - mpcom_1.b * c__;
	}
	if (c__ != 0) {
	    goto L90;
	}
	c__ = 8;
L40:
	;
    }
    if (c__ == 8) {
	goto L60;
    }
    if (xi < 0 || xi >= mpcom_1.b) {
	goto L90;
    }
    c__ = 0;
    i__1 = i2;
    for (j = 1; j <= i__1; ++j) {
	j1 = i2p - j;
	ri = mpcom_1.r__[j1 - 1] + c__;
	if (ri < 0) {
	    goto L70;
	}
	c__ = ri / mpcom_1.b;
/* L50: */
	mpcom_1.r__[j1 - 1] = ri - mpcom_1.b * c__;
    }
    if (c__ != 0) {
	goto L90;
    }
/* NORMALIZE AND ROUND RESULT */
L60:
    mpnzr_(&rs, &re, &z__[1], &c__0);
    return 0;
L70:
    io___11.ciunit = mpcom_1.lun;
    s_wsfe(&io___11);
    e_wsfe();
    goto L110;
L90:
    io___12.ciunit = mpcom_1.lun;
    s_wsfe(&io___12);
    e_wsfe();
L110:
    mperr_();
    z__[1] = 0;
    return 0;
} /* mpmul_ */

