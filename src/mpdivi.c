/* mpdivi.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;
static integer c__1 = 1;
static integer c__4 = 4;

/* DECK MPDIVI */
/* Subroutine */ int mpdivi_(integer *x, integer *iy, integer *z__)
{
    /* Format strings */
    static char fmt_20[] = "(\002 *** ATTEMPTED DIVISION BY ZERO IN CALL TO "
	    "MPDIVI ***\002)";
    static char fmt_220[] = "(\002 *** INTEGER OVERFLOW IN MPDIVI, B TOO LAR"
	    "GE ***\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer c__, i__, j, k, b2, c2, i2, j1, j2, r1, j11, kh, re, iq, 
	    ir, rs, iqj;
    extern /* Subroutine */ int mpchk_(integer *, integer *), mperr_(void), 
	    mpstr_(integer *, integer *), mpnzr_(integer *, integer *, 
	    integer *, integer *), mpunfl_(integer *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_220, 0 };


/* ***BEGIN PROLOGUE  MPDIVI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPDIVI-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Divides 'mp' X by the single-precision integer IY giving 'mp' Z. */
/*  This is much faster than division by an 'mp' number. */

/*  The arguments X(*) and Z(*), and the variable R in COMMON are all */
/*  INTEGER arrays of size 30.  See the comments in the routine MPBLAS */
/*  for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPCHK, MPERR, MPNZR, MPSTR, MPUNFL */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPDIVI */
/* ***FIRST EXECUTABLE STATEMENT  MPDIVI */
    /* Parameter adjustments */
    --z__;
    --x;

    /* Function Body */
    rs = x[1];
    j = *iy;
    if (j < 0) {
	goto L30;
    } else if (j == 0) {
	goto L10;
    } else {
	goto L40;
    }
L10:
    io___3.ciunit = mpcom_1.lun;
    s_wsfe(&io___3);
    e_wsfe();
    goto L230;
L30:
    j = -j;
    rs = -rs;
L40:
    re = x[2];
/* CHECK FOR ZERO DIVIDEND */
    if (rs == 0) {
	goto L120;
    }
/* CHECK FOR DIVISION BY B */
    if (j != mpcom_1.b) {
	goto L50;
    }
    mpstr_(&x[1], &z__[1]);
    if (re <= -mpcom_1.m) {
	goto L240;
    }
    z__[1] = rs;
    z__[2] = re - 1;
    return 0;
/* CHECK FOR DIVISION BY 1 OR -1 */
L50:
    if (j != 1) {
	goto L60;
    }
    mpstr_(&x[1], &z__[1]);
    z__[1] = rs;
    return 0;
L60:
    c__ = 0;
    i2 = mpcom_1.t + 4;
    i__ = 0;
/* IF J*B NOT REPRESENTABLE AS AN INTEGER HAVE TO SIMULATE */
/* LONG DIVISION.   ASSUME AT LEAST 16-BIT WORD. */
/* Computing MAX */
    i__1 = mpcom_1.b << 3, i__2 = 32767 / mpcom_1.b;
    b2 = max(i__1,i__2);
    if (j >= b2) {
	goto L130;
    }
/* LOOK FOR FIRST NONZERO DIGIT IN QUOTIENT */
L70:
    ++i__;
    c__ = mpcom_1.b * c__;
    if (i__ <= mpcom_1.t) {
	c__ += x[i__ + 2];
    }
    r1 = c__ / j;
    if (r1 < 0) {
	goto L210;
    } else if (r1 == 0) {
	goto L70;
    } else {
	goto L80;
    }
/* ADJUST EXPONENT AND GET T+4 DIGITS IN QUOTIENT */
L80:
    re = re + 1 - i__;
    mpcom_1.r__[0] = r1;
    c__ = mpcom_1.b * (c__ - j * r1);
    kh = 2;
    if (i__ >= mpcom_1.t) {
	goto L100;
    }
    kh = mpcom_1.t + 1 - i__;
    i__1 = kh;
    for (k = 2; k <= i__1; ++k) {
	++i__;
	c__ += x[i__ + 2];
	mpcom_1.r__[k - 1] = c__ / j;
/* L90: */
	c__ = mpcom_1.b * (c__ - j * mpcom_1.r__[k - 1]);
    }
    if (c__ < 0) {
	goto L210;
    }
    ++kh;
L100:
    i__1 = i2;
    for (k = kh; k <= i__1; ++k) {
	mpcom_1.r__[k - 1] = c__ / j;
/* L110: */
	c__ = mpcom_1.b * (c__ - j * mpcom_1.r__[k - 1]);
    }
    if (c__ < 0) {
	goto L210;
    }
/* NORMALIZE AND ROUND RESULT */
L120:
    mpnzr_(&rs, &re, &z__[1], &c__0);
    return 0;
/* HERE NEED SIMULATED DOUBLE-PRECISION DIVISION */
L130:
    c2 = 0;
    j1 = j / mpcom_1.b;
    j2 = j - j1 * mpcom_1.b;
    j11 = j1 + 1;
/* LOOK FOR FIRST NONZERO DIGIT */
L140:
    ++i__;
    c__ = mpcom_1.b * c__ + c2;
    c2 = 0;
    if (i__ <= mpcom_1.t) {
	c2 = x[i__ + 2];
    }
    if ((i__1 = c__ - j1) < 0) {
	goto L140;
    } else if (i__1 == 0) {
	goto L150;
    } else {
	goto L160;
    }
L150:
    if (c2 < j2) {
	goto L140;
    }
/* COMPUTE T+4 QUOTIENT DIGITS */
L160:
    re = re + 1 - i__;
    k = 1;
    goto L180;
/* MAIN LOOP FOR LARGE ABS(IY) CASE */
L170:
    ++k;
    if (k > i2) {
	goto L120;
    }
    ++i__;
/* GET APPROXIMATE QUOTIENT FIRST */
L180:
    ir = c__ / j11;
/* NOW REDUCE SO OVERFLOW DOES NOT OCCUR */
    iq = c__ - ir * j1;
    if (iq < b2) {
	goto L190;
    }
/* HERE IQ*B WOULD POSSIBLY OVERFLOW SO INCREASE IR */
    ++ir;
    iq -= j1;
L190:
    iq = iq * mpcom_1.b - ir * j2;
    if (iq >= 0) {
	goto L200;
    }
/* HERE IQ NEGATIVE SO IR WAS TOO LARGE */
    --ir;
    iq += j;
L200:
    if (i__ <= mpcom_1.t) {
	iq += x[i__ + 2];
    }
    iqj = iq / j;
/* R(K) = QUOTIENT, C = REMAINDER */
    mpcom_1.r__[k - 1] = iqj + ir;
    c__ = iq - j * iqj;
    if (c__ >= 0) {
	goto L170;
    }
/* CARRY NEGATIVE SO OVERFLOW MUST HAVE OCCURRED */
L210:
    mpchk_(&c__1, &c__4);
    io___19.ciunit = mpcom_1.lun;
    s_wsfe(&io___19);
    e_wsfe();
L230:
    mperr_();
    z__[1] = 0;
    return 0;
/* UNDERFLOW HERE */
L240:
    mpunfl_(&z__[1]);
    return 0;
} /* mpdivi_ */

