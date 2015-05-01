/* mpmul2.f -- translated by f2c (version 12.02.01).
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

/* DECK MPMUL2 */
/* Subroutine */ int mpmul2_(integer *x, integer *iy, integer *z__, integer *
	trunc)
{
    /* Format strings */
    static char fmt_30[] = "(\002 *** OVERFLOW OCCURRED IN MPMUL2 ***\002)";
    static char fmt_140[] = "(\002 *** INTEGER OVERFLOW IN MPMUL2, B TOO LAR"
	    "GE ***\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer c__, i__, j, c1, c2, j1, j2, t1, t3, t4, ij, re, ri, is, 
	    ix, rs;
    extern /* Subroutine */ int mpchk_(integer *, integer *), mperr_(void), 
	    mpstr_(integer *, integer *), mpnzr_(integer *, integer *, 
	    integer *, integer *), mpovfl_(integer *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_140, 0 };


/* ***BEGIN PROLOGUE  MPMUL2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPMUL2-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Multiplies 'mp' X by single-precision integer IY giving 'mp' Z. */
/*  Multiplication by 1 may be used to normalize a number even if some */
/*  digits are greater than B-1. Result is rounded if TRUNC.EQ.0, */
/*  otherwise truncated. */

/*  The arguments X(*) and Z(*), and the variable R in COMMON are all */
/*  INTEGER arrays of size 30.  See the comments in the routine MPBLAS */
/*  for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPCHK, MPERR, MPNZR, MPOVFL, MPSTR */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPMUL2 */
/* ***FIRST EXECUTABLE STATEMENT  MPMUL2 */
    /* Parameter adjustments */
    --z__;
    --x;

    /* Function Body */
    rs = x[1];
    if (rs == 0) {
	goto L10;
    }
    j = *iy;
    if (j < 0) {
	goto L20;
    } else if (j == 0) {
	goto L10;
    } else {
	goto L50;
    }
/* RESULT ZERO */
L10:
    z__[1] = 0;
    return 0;
L20:
    j = -j;
    rs = -rs;
/* CHECK FOR MULTIPLICATION BY B */
    if (j != mpcom_1.b) {
	goto L50;
    }
    if (x[2] < mpcom_1.m) {
	goto L40;
    }
    mpchk_(&c__1, &c__4);
    io___3.ciunit = mpcom_1.lun;
    s_wsfe(&io___3);
    e_wsfe();
    mpovfl_(&z__[1]);
    return 0;
L40:
    mpstr_(&x[1], &z__[1]);
    z__[1] = rs;
    z__[2] = x[2] + 1;
    return 0;
/* SET EXPONENT TO EXPONENT(X) + 4 */
L50:
    re = x[2] + 4;
/* FORM PRODUCT IN ACCUMULATOR */
    c__ = 0;
    t1 = mpcom_1.t + 1;
    t3 = mpcom_1.t + 3;
    t4 = mpcom_1.t + 4;
/* IF J*B NOT REPRESENTABLE AS AN INTEGER WE HAVE TO SIMULATE */
/* DOUBLE-PRECISION MULTIPLICATION. */
/* Computing MAX */
    i__1 = mpcom_1.b << 3, i__2 = 32767 / mpcom_1.b;
    if (j >= max(i__1,i__2)) {
	goto L110;
    }
    i__1 = mpcom_1.t;
    for (ij = 1; ij <= i__1; ++ij) {
	i__ = t1 - ij;
	ri = j * x[i__ + 2] + c__;
	c__ = ri / mpcom_1.b;
/* L60: */
	mpcom_1.r__[i__ + 3] = ri - mpcom_1.b * c__;
    }
/* CHECK FOR INTEGER OVERFLOW */
    if (ri < 0) {
	goto L130;
    }
/* HAVE TO TREAT FIRST FOUR WORDS OF R SEPARATELY */
    for (ij = 1; ij <= 4; ++ij) {
	i__ = 5 - ij;
	ri = c__;
	c__ = ri / mpcom_1.b;
/* L70: */
	mpcom_1.r__[i__ - 1] = ri - mpcom_1.b * c__;
    }
    if (c__ == 0) {
	goto L100;
    }
/* HAVE TO SHIFT RIGHT HERE AS CARRY OFF END */
L80:
    i__1 = t3;
    for (ij = 1; ij <= i__1; ++ij) {
	i__ = t4 - ij;
/* L90: */
	mpcom_1.r__[i__] = mpcom_1.r__[i__ - 1];
    }
    ri = c__;
    c__ = ri / mpcom_1.b;
    mpcom_1.r__[0] = ri - mpcom_1.b * c__;
    ++re;
    if (c__ < 0) {
	goto L130;
    } else if (c__ == 0) {
	goto L100;
    } else {
	goto L80;
    }
/* NORMALIZE AND ROUND OR TRUNCATE RESULT */
L100:
    mpnzr_(&rs, &re, &z__[1], trunc);
    return 0;
/* HERE J IS TOO LARGE FOR SINGLE-PRECISION MULTIPLICATION */
L110:
    j1 = j / mpcom_1.b;
    j2 = j - j1 * mpcom_1.b;
/* FORM PRODUCT */
    i__1 = t4;
    for (ij = 1; ij <= i__1; ++ij) {
	c1 = c__ / mpcom_1.b;
	c2 = c__ - mpcom_1.b * c1;
	i__ = t1 - ij;
	ix = 0;
	if (i__ > 0) {
	    ix = x[i__ + 2];
	}
	ri = j2 * ix + c2;
	is = ri / mpcom_1.b;
	c__ = j1 * ix + c1 + is;
/* L120: */
	mpcom_1.r__[i__ + 3] = ri - mpcom_1.b * is;
    }
    if (c__ < 0) {
	goto L130;
    } else if (c__ == 0) {
	goto L100;
    } else {
	goto L80;
    }
/* CAN ONLY GET HERE IF INTEGER OVERFLOW OCCURRED */
L130:
    mpchk_(&c__1, &c__4);
    io___18.ciunit = mpcom_1.lun;
    s_wsfe(&io___18);
    e_wsfe();
    mperr_();
    goto L10;
} /* mpmul2_ */

