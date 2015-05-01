/* atanh.f -- translated by f2c (version 12.02.01).
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
static integer c__15 = 15;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK ATANH */
doublereal atanh_(real *x)
{
    /* Initialized data */

    static real atnhcs[15] = { .094395102393195492f,.049198437055786159f,
	    .002102593522455432f,1.07355444977611e-4f,5.978267249293e-6f,
	    3.50506203088e-7f,2.1263743437e-8f,1.321694535e-9f,8.3658755e-11f,
	    5.370503e-12f,3.48665e-13f,2.2845e-14f,1.508e-15f,1e-16f,6e-18f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y;
    extern doublereal csevl_(real *, real *, integer *);
    static real dxrel;
    extern integer inits_(real *, integer *, real *);
    static real sqeps;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;

/* ***BEGIN PROLOGUE  ATANH */
/* ***PURPOSE  Compute the arc hyperbolic tangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      SINGLE PRECISION (ATANH-S, DATANH-D, CATANH-C) */
/* ***KEYWORDS  ARC HYPERBOLIC TANGENT, ATANH, ELEMENTARY FUNCTIONS, */
/*             FNLIB, INVERSE HYPERBOLIC TANGENT */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ATANH(X) computes the arc hyperbolic tangent of X. */

/* Series for ATNH       on the interval  0.          to  2.50000D-01 */
/*                                        with weighted error   6.70E-18 */
/*                                         log weighted error  17.17 */
/*                               significant figures required  16.01 */
/*                                    decimal places required  17.76 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  ATANH */
/* ***FIRST EXECUTABLE STATEMENT  ATANH */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nterms = inits_(atnhcs, &c__15, &r__1);
	dxrel = sqrt(r1mach_(&c__4));
	sqeps = sqrt(r1mach_(&c__3) * 3.f);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y >= 1.f) {
	xermsg_("SLATEC", "ATANH", "ABS(X) GE 1", &c__2, &c__2, (ftnlen)6, (
		ftnlen)5, (ftnlen)11);
    }

    if (1.f - y < dxrel) {
	xermsg_("SLATEC", "ATANH", "ANSWER LT HALF PRECISION BECAUSE ABS(X) "
		"TOO NEAR 1", &c__1, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)50);
    }

    ret_val = *x;
    if (y > sqeps && y <= .5f) {
	r__1 = *x * 8.f * *x - 1.f;
	ret_val = *x * (csevl_(&r__1, atnhcs, &nterms) + 1.f);
    }
    if (y > .5f) {
	ret_val = log((*x + 1.f) / (1.f - *x)) * .5f;
    }

    return ret_val;
} /* atanh_ */

