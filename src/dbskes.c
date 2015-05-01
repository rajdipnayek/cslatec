/* dbskes.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static real c_b18 = 1.f;
static doublereal c_b19 = 1.;
static integer c__5 = 5;

/* DECK DBSKES */
/* Subroutine */ int dbskes_(doublereal *xnu, doublereal *x, integer *nin, 
	doublereal *bke)
{
    /* Initialized data */

    static doublereal alnbig = 0.;

    /* System generated locals */
    integer i__1;
    real r__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, n;
    static doublereal v, vend, bknu1, vincr;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int d9knus_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal direct;
    static integer iswtch;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBSKES */
/* ***PURPOSE  Compute a sequence of exponentially scaled modified Bessel */
/*            functions of the third kind of fractional order. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B3 */
/* ***TYPE      DOUBLE PRECISION (BESKES-S, DBSKES-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, FRACTIONAL ORDER, */
/*             MODIFIED BESSEL FUNCTION, SEQUENCE OF BESSEL FUNCTIONS, */
/*             SPECIAL FUNCTIONS, THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBSKES(XNU,X,NIN,BKE) computes a double precision sequence */
/* of exponentially scaled modified Bessel functions */
/* of the third kind of order XNU + I at X, where X .GT. 0, */
/* XNU lies in (-1,1), and I = 0, 1, ... , NIN - 1, if NIN is positive */
/* and I = 0, -1, ... , NIN + 1, if NIN is negative.  On return, the */
/* vector BKE(.) contains the results at X for order starting at XNU. */
/* XNU, X, and BKE are double precision.  NIN is integer. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9KNUS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBSKES */
    /* Parameter adjustments */
    --bke;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DBSKES */
    if (alnbig == 0.) {
	alnbig = log(d1mach_(&c__2));
    }

    v = abs(*xnu);
    n = abs(*nin);

    if (v >= 1.) {
	xermsg_("SLATEC", "DBSKES", "ABS(XNU) MUST BE LT 1", &c__2, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x <= 0.) {
	xermsg_("SLATEC", "DBSKES", "X IS LE 0", &c__3, &c__2, (ftnlen)6, (
		ftnlen)6, (ftnlen)9);
    }
    if (n == 0) {
	xermsg_("SLATEC", "DBSKES", "N THE NUMBER IN THE SEQUENCE IS 0", &
		c__4, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)33);
    }

    d9knus_(&v, x, &bke[1], &bknu1, &iswtch);
    if (n == 1) {
	return 0;
    }

    r__1 = (real) (*nin);
    vincr = r_sign(&c_b18, &r__1);
    direct = vincr;
    if (*xnu != 0.) {
	direct = vincr * d_sign(&c_b19, xnu);
    }
    if (iswtch == 1 && direct > 0.f) {
	xermsg_("SLATEC", "DBSKES", "X SO SMALL BESSEL K-SUB-XNU+1 OVERFLOWS",
		 &c__5, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)39);
    }
    bke[2] = bknu1;

    if (direct < 0.f) {
	d__2 = (d__1 = *xnu + vincr, abs(d__1));
	d9knus_(&d__2, x, &bke[2], &bknu1, &iswtch);
    }
    if (n == 2) {
	return 0;
    }

    vend = (d__1 = *xnu + *nin, abs(d__1)) - 1.;
    if ((vend - .5) * log(vend) + .27 - vend * (log(*x) - .694) > alnbig) {
	xermsg_("SLATEC", "DBSKES", "X SO SMALL OR ABS(NU) SO BIG THAT BESSE"
		"L K-SUB-NU OVERFLOWS", &c__5, &c__2, (ftnlen)6, (ftnlen)6, (
		ftnlen)59);
    }

    v = *xnu;
    i__1 = n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	v += vincr;
	bke[i__] = v * 2. * bke[i__ - 1] / *x + bke[i__ - 2];
/* L10: */
    }

    return 0;
} /* dbskes_ */

