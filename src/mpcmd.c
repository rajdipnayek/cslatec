/* mpcmd.f -- translated by f2c (version 12.02.01).
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

/* DECK MPCMD */
/* Subroutine */ int mpcmd_(integer *x, doublereal *dz)
{
    /* Format strings */
    static char fmt_40[] = "(\002 *** FLOATING-POINT OVER/UNDER-FLOW IN MPCM"
	    "D ***\002)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal db;
    static integer tm;
    static doublereal dz2;
    extern /* Subroutine */ int mpchk_(integer *, integer *), mperr_(void);

    /* Fortran I/O blocks */
    static cilist io___5 = { 0, 0, 0, fmt_40, 0 };


/* ***BEGIN PROLOGUE  MPCMD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DQDOTA and DQDOTI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (MPCMD-A) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*  Converts multiple-precision X to double-precision DZ. Assumes */
/*  X is in allowable range for double-precision numbers. There is */
/*  some loss of accuracy if the exponent is large. */

/*  The argument X(*) is INTEGER array of size 30.  See the comments in */
/*  the routine MPBLAS for the reason for this choice. */

/* ***SEE ALSO  DQDOTA, DQDOTI, MPBLAS */
/* ***ROUTINES CALLED  MPCHK, MPERR */
/* ***COMMON BLOCKS    MPCOM */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791001  DATE WRITTEN */
/*   ??????  Modified for use with BLAS.  Blank COMMON changed to named */
/*           COMMON.  R given dimension 12. */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/*   930124  Increased Array size in MPCON for SUN -r8.  (RWC) */
/* ***END PROLOGUE  MPCMD */
/* ***FIRST EXECUTABLE STATEMENT  MPCMD */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    mpchk_(&c__1, &c__4);
    *dz = 0.;
    if (x[1] == 0) {
	return 0;
    }
    db = (doublereal) mpcom_1.b;
    i__1 = mpcom_1.t;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*dz = db * *dz + (doublereal) x[i__ + 2];
	tm = i__;
/* CHECK IF FULL DOUBLE-PRECISION ACCURACY ATTAINED */
	dz2 = *dz + 1.;
/* TEST BELOW NOT ALWAYS EQUIVALENT TO - IF (DZ2.LE.DZ) GO TO 20, */
/* FOR EXAMPLE ON CYBER 76. */
	if (dz2 - *dz <= 0.) {
	    goto L20;
	}
/* L10: */
    }
/* NOW ALLOW FOR EXPONENT */
L20:
    i__1 = x[2] - tm;
    *dz *= pow_di(&db, &i__1);
/* CHECK REASONABLENESS OF RESULT. */
    if (*dz <= 0.) {
	goto L30;
    }
/* LHS SHOULD BE .LE. 0.5 BUT ALLOW FOR SOME ERROR IN LOG */
    if ((d__1 = (doublereal) x[2] - (log(*dz) / log((doublereal) mpcom_1.b) + 
	    .5), abs(d__1)) > .6) {
	goto L30;
    }
    if (x[1] < 0) {
	*dz = -(*dz);
    }
    return 0;
/* FOLLOWING MESSAGE INDICATES THAT X IS TOO LARGE OR SMALL - */
/* TRY USING MPCMDE INSTEAD. */
L30:
    io___5.ciunit = mpcom_1.lun;
    s_wsfe(&io___5);
    e_wsfe();
    mperr_();
    return 0;
} /* mpcmd_ */

