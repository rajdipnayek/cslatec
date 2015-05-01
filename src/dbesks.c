/* dbesks.f -- translated by f2c (version 12.02.01).
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

/* DECK DBESKS */
/* Subroutine */ int dbesks_(doublereal *xnu, doublereal *x, integer *nin, 
	doublereal *bk)
{
    /* Initialized data */

    static doublereal xmax = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, n;
    static doublereal expxi;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int dbskes_(doublereal *, doublereal *, integer *,
	     doublereal *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBESKS */
/* ***PURPOSE  Compute a sequence of modified Bessel functions of the */
/*            third kind of fractional order. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B3 */
/* ***TYPE      DOUBLE PRECISION (BESKS-S, DBESKS-D) */
/* ***KEYWORDS  FNLIB, FRACTIONAL ORDER, MODIFIED BESSEL FUNCTION, */
/*             SEQUENCE OF BESSEL FUNCTIONS, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESKS computes a sequence of modified Bessel functions of the third */
/* kind of order XNU + I at X, where X .GT. 0, XNU lies in (-1,1), */
/* and I = 0, 1, ... , NIN - 1, if NIN is positive and I = 0, 1, ... , */
/* NIN + 1, if NIN is negative.  On return, the vector BK(.) contains */
/* the results at X for order starting at XNU.  XNU, X, and BK are */
/* double precision.  NIN is an integer. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DBSKES, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESKS */
    /* Parameter adjustments */
    --bk;

    /* Function Body */
/* ***FIRST EXECUTABLE STATEMENT  DBESKS */
    if (xmax == 0.) {
	xmax = -log(d1mach_(&c__1));
    }

    if (*x > xmax) {
	xermsg_("SLATEC", "DBESKS", "X SO BIG BESSEL K UNDERFLOWS", &c__1, &
		c__2, (ftnlen)6, (ftnlen)6, (ftnlen)28);
    }

    dbskes_(xnu, x, nin, &bk[1]);

    expxi = exp(-(*x));
    n = abs(*nin);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bk[i__] = expxi * bk[i__];
/* L20: */
    }

    return 0;
} /* dbesks_ */

