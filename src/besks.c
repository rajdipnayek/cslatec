/* besks.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESKS */
/* Subroutine */ int besks_(real *xnu, real *x, integer *nin, real *bk)
{
    /* Initialized data */

    static real xmax = 0.f;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n;
    static real expxi;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int beskes_(real *, real *, integer *, real *), 
	    xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESKS */
/* ***PURPOSE  Compute a sequence of modified Bessel functions of the */
/*            third kind of fractional order. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B3 */
/* ***TYPE      SINGLE PRECISION (BESKS-S, DBESKS-D) */
/* ***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION, */
/*             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESKS computes a sequence of modified Bessel functions of the third */
/* kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1), */
/* and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... , */
/* NIN + 1, if NIN is negative.  On return, the vector BK(.) Contains */
/* the results at X for order starting at XNU. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESKES, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESKS */
    /* Parameter adjustments */
    --bk;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  BESKS */
    if (xmax == 0.f) {
	xmax = -log(r1mach_(&c__1));
    }

    if (*x > xmax) {
	xermsg_("SLATEC", "BESKS", "X SO BIG BESSEL K UNDERFLOWS", &c__1, &
		c__2, (ftnlen)6, (ftnlen)5, (ftnlen)28);
    }

    beskes_(xnu, x, nin, &bk[1]);

    expxi = exp(-(*x));
    n = abs(*nin);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bk[i__] = expxi * bk[i__];
/* L20: */
    }

    return 0;
} /* besks_ */

