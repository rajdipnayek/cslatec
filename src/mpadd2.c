/* mpadd2.f -- translated by f2c (version 12.02.01).
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

/* DECK MPADD2 */
/* Subroutine */ int mpadd2_(integer *x, integer *y, integer *z__, integer *
	y1, integer *trunc)
{
    /* Format strings */
    static char fmt_50[] = "(\002 *** SIGN NOT 0, +1 OR -1 IN CALL TO MPAD"
	    "D2,\002,\002 POSSIBLE OVERWRITING PROBLEM ***\002)";

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, s, ed, re, rs, med;
    extern /* Subroutine */ int mpchk_(integer *, integer *), mperr_(void), 
	    mpstr_(integer *, integer *), mpnzr_(integer *, integer *, 
	    integer *, integer *), mpadd3_(integer *, integer *, integer *, 
	    integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_50, 0 };


/* ***BEGIN PROLOGUE  MPADD2 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPADD2-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Called by MPADD, MPSUB etc. */
/*  X, Y and Z are MP numbers, Y1 and TRUNC are integers. */
/*  To force call by reference rather than value/result, Y1 is */
/*  declared as an array, but only Y1(1) is ever used. */
/*  Sets Z = X + Y1(1)*ABS(Y), where Y1(1) = +- Y(1). */
/*  If TRUNC .EQ. 0, R*-rounding is used;  otherwise, truncation. */
/*  R*-rounding is defined in the Kuki and Cody reference. */

/*  The arguments X(*), Y(*), and Z(*) are all INTEGER arrays of size */
/*  30.  See the comments in the routine MPBLAS for the reason for this */
/*  choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***REFERENCES  H. Kuki and W. J. Cody, A statistical study of floating */
/*                 point number systems, Communications of the ACM 16, 4 */
/*                 (April 1973), pp. 223-230. */
/*               R. P. Brent, On the precision attainable with various */
/*                 floating-point number systems, IEEE Transactions on */
/*                 Computers C-22, 6 (June 1973), pp. 601-607. */
/*               R. P. Brent, A Fortran multiple-precision arithmetic */
/*                 package, ACM Transactions on Mathematical Software 4, */
/*                 1 (March 1978), pp. 57-70. */
/*               R. P. Brent, MP, a Fortran multiple-precision arithmetic */
/*                 package, Algorithm 524, ACM Transactions on Mathema- */
/*                 tical Software 4, 1 (March 1978), pp. 71-81. */
/* ***ROUTINES CALLED  MPADD3, MPCHK, MPERR, MPNZR, MPSTR */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   920528  Added a REFERENCES section revised.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPADD2 */
/* ***FIRST EXECUTABLE STATEMENT  MPADD2 */
    /* Parameter adjustments */
    --y1;
    --z__;
    --y;
    --x;

    /* Function Body */
    if (x[1] != 0) {
	goto L20;
    }
L10:
    mpstr_(&y[1], &z__[1]);
    z__[1] = y1[1];
    return 0;
L20:
    if (y1[1] != 0) {
	goto L40;
    }
L30:
    mpstr_(&x[1], &z__[1]);
    return 0;
/* COMPARE SIGNS */
L40:
    s = x[1] * y1[1];
    if (abs(s) <= 1) {
	goto L60;
    }
    mpchk_(&c__1, &c__4);
    io___2.ciunit = mpcom_1.lun;
    s_wsfe(&io___2);
    e_wsfe();
    mperr_();
    z__[1] = 0;
    return 0;
/* COMPARE EXPONENTS */
L60:
    ed = x[2] - y[2];
    med = abs(ed);
    if (ed < 0) {
	goto L90;
    } else if (ed == 0) {
	goto L70;
    } else {
	goto L120;
    }
/* EXPONENTS EQUAL SO COMPARE SIGNS, THEN FRACTIONS IF NEC. */
L70:
    if (s > 0) {
	goto L100;
    }
    i__1 = mpcom_1.t;
    for (j = 1; j <= i__1; ++j) {
	if ((i__2 = x[j + 2] - y[j + 2]) < 0) {
	    goto L100;
	} else if (i__2 == 0) {
	    goto L80;
	} else {
	    goto L130;
	}
L80:
	;
    }
/* RESULT IS ZERO */
    z__[1] = 0;
    return 0;
/* HERE EXPONENT(Y) .GE. EXPONENT(X) */
L90:
    if (med > mpcom_1.t) {
	goto L10;
    }
L100:
    rs = y1[1];
    re = y[2];
    mpadd3_(&x[1], &y[1], &s, &med, &re);
/* NORMALIZE, ROUND OR TRUNCATE, AND RETURN */
L110:
    mpnzr_(&rs, &re, &z__[1], trunc);
    return 0;
/* ABS(X) .GT. ABS(Y) */
L120:
    if (med > mpcom_1.t) {
	goto L30;
    }
L130:
    rs = x[1];
    re = x[2];
    mpadd3_(&y[1], &x[1], &s, &med, &re);
    goto L110;
} /* mpadd2_ */

