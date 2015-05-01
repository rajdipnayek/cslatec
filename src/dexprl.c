/* dexprl.f -- translated by f2c (version 12.02.01).
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

/* DECK DEXPRL */
doublereal dexprl_(doublereal *x)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__;
    static doublereal xn, xln, xbnd, absx;
    extern doublereal d1mach_(integer *);
    static doublereal alneps;
    static integer nterms;

/* ***BEGIN PROLOGUE  DEXPRL */
/* ***PURPOSE  Calculate the relative error exponential (EXP(X)-1)/X. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      DOUBLE PRECISION (EXPREL-S, DEXPRL-D, CEXPRL-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, EXPONENTIAL, FIRST ORDER, FNLIB */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  EXPREL(X) = (EXP(X) - 1.0) / X.   For small ABS(X) the */
/* Taylor series is used.  If X is negative the reflection formula */
/*         EXPREL(X) = EXP(X) * EXPREL(ABS(X)) */
/* may be used.  This reflection formula will be of use when the */
/* evaluation for small ABS(X) is done by Chebyshev series rather than */
/* Taylor series. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DEXPRL */
/* ***FIRST EXECUTABLE STATEMENT  DEXPRL */
    if (first) {
	alneps = log(d1mach_(&c__3));
	xn = 3.72 - alneps * .3;
	xln = log((xn + 1.) / 1.36);
	nterms = (integer) (xn - (xn * xln + alneps) / (xln + 1.36) + 1.5);
	xbnd = d1mach_(&c__3);
    }
    first = FALSE_;

    absx = abs(*x);
    if (absx > .5) {
	ret_val = (exp(*x) - 1.) / *x;
    }
    if (absx > .5) {
	return ret_val;
    }

    ret_val = 1.;
    if (absx < xbnd) {
	return ret_val;
    }

    ret_val = 0.;
    i__1 = nterms;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val = ret_val * *x / (nterms + 2 - i__) + 1.;
/* L20: */
    }

    return ret_val;
} /* dexprl_ */

