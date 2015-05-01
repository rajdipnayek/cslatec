/* alngam.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__1 = 1;

/* DECK ALNGAM */
doublereal alngam_(real *x)
{
    /* Initialized data */

    static real sq2pil = .91893853320467274f;
    static real sqpi2l = .22579135264472743f;
    static real pi = 3.14159265358979324f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real y, xmax;
    extern doublereal gamma_(real *);
    static real dxrel;
    extern doublereal r1mach_(integer *), r9lgmc_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real sinpiy;

/* ***BEGIN PROLOGUE  ALNGAM */
/* ***PURPOSE  Compute the logarithm of the absolute value of the Gamma */
/*            function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      SINGLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C) */
/* ***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ALNGAM(X) computes the logarithm of the absolute value of the */
/* gamma function at X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  GAMMA, R1MACH, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  ALNGAM */
/* ***FIRST EXECUTABLE STATEMENT  ALNGAM */
    if (first) {
	xmax = r1mach_(&c__2) / log(r1mach_(&c__2));
	dxrel = sqrt(r1mach_(&c__4));
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 10.f) {
	goto L20;
    }

/* LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0 */

    ret_val = log((r__1 = gamma_(x), dabs(r__1)));
    return ret_val;

/* LOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0 */

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "ALNGAM", "ABS(X) SO BIG ALNGAM OVERFLOWS", &c__2, &
		c__2, (ftnlen)6, (ftnlen)6, (ftnlen)30);
    }

    if (*x > 0.f) {
	ret_val = sq2pil + (*x - .5f) * log(*x) - *x + r9lgmc_(&y);
    }
    if (*x > 0.f) {
	return ret_val;
    }

    sinpiy = (r__1 = sin(pi * y), dabs(r__1));
    if (sinpiy == 0.f) {
	xermsg_("SLATEC", "ALNGAM", "X IS A NEGATIVE INTEGER", &c__3, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)23);
    }

    r__2 = *x - .5f;
    if ((r__1 = (*x - r_int(&r__2)) / *x, dabs(r__1)) < dxrel) {
	xermsg_("SLATEC", "ALNGAM", "ANSWER LT HALF PRECISION BECAUSE X TOO "
		"NEAR NEGATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (
		ftnlen)60);
    }

    ret_val = sqpi2l + (*x - .5f) * log(y) - *x - log(sinpiy) - r9lgmc_(&y);
    return ret_val;

} /* alngam_ */

