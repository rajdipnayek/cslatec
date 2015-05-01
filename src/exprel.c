/* exprel.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;

/* DECK EXPREL */
doublereal exprel_(real *x)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__;
    static real xn, xln, xbnd, absx;
    extern doublereal r1mach_(integer *);
    static real alneps;
    static integer nterms;

/* ***BEGIN PROLOGUE  EXPREL */
/* ***PURPOSE  Calculate the relative error exponential (EXP(X)-1)/X. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      SINGLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the */
/* Taylor series is used.  If X is negative, the reflection formula */
/*         EXPREL(X) = EXP(X) * EXPREL(ABS(X)) */
/* may be used.  This reflection formula will be of use when the */
/* evaluation for small ABS(X) is done by Chebyshev series rather than */
/* Taylor series.  EXPREL and X are single precision. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  EXPREL */
/* ***FIRST EXECUTABLE STATEMENT  EXPREL */
    if (first) {
	alneps = log(r1mach_(&c__3));
	xn = 3.72f - alneps * .3f;
	xln = log((xn + 1.f) / 1.36f);
	nterms = xn - (xn * xln + alneps) / (xln + 1.36f) + 1.5f;
	xbnd = r1mach_(&c__3);
    }
    first = FALSE_;

    absx = dabs(*x);
    if (absx > .5f) {
	ret_val = (exp(*x) - 1.f) / *x;
    }
    if (absx > .5f) {
	return ret_val;
    }

    ret_val = 1.f;
    if (absx < xbnd) {
	return ret_val;
    }

    ret_val = 0.f;
    i__1 = nterms;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val = ret_val * *x / (nterms + 2 - i__) + 1.f;
/* L20: */
    }

    return ret_val;
} /* exprel_ */

