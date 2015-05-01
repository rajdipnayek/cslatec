/* gamma.f -- translated by f2c (version 12.02.01).
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
static integer c__23 = 23;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK GAMMA */
doublereal gamma_(real *x)
{
    /* Initialized data */

    static real gcs[23] = { .008571195590989331f,.004415381324841007f,
	    .05685043681599363f,-.004219835396418561f,.00132680818121246f,
	    -1.89302452979888e-4f,3.60692532744124e-5f,-6.0567619044608e-6f,
	    1.0558295463022e-6f,-1.811967365542e-7f,3.11772496471e-8f,
	    -5.354219639e-9f,9.193275519e-10f,-1.57794128e-10f,
	    2.70798062e-11f,-4.6468186e-12f,7.97335e-13f,-1.368078e-13f,
	    2.34731e-14f,-4.0274e-15f,6.91e-16f,-1.185e-16f,2.03e-17f };
    static real pi = 3.14159265358979324f;
    static real sq2pil = .91893853320467274f;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2;

    /* Local variables */
    static integer i__, n;
    static real y;
    static integer ngcs;
    static real xmin, xmax;
    extern doublereal csevl_(real *, real *, integer *);
    static real dxrel;
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *), r9lgmc_(real *);
    extern /* Subroutine */ int gamlim_(real *, real *), xermsg_(char *, char 
	    *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    static real sinpiy;

/* ***BEGIN PROLOGUE  GAMMA */
/* ***PURPOSE  Compute the complete Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      SINGLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* GAMMA computes the gamma function at X, where X is not 0, -1, -2, .... */
/* GAMMA and X are single precision. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, GAMLIM, INITS, R1MACH, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  GAMMA */
/* SQ2PIL IS LOG (SQRT (2.*PI) ) */

/* LANL DEPENDENT CODE REMOVED 81.02.04 */

/* ***FIRST EXECUTABLE STATEMENT  GAMMA */
    if (first) {

/* --------------------------------------------------------------------- */
/* INITIALIZE.  FIND LEGAL BOUNDS FOR X, AND DETERMINE THE NUMBER OF */
/* TERMS IN THE SERIES REQUIRED TO ATTAIN AN ACCURACY TEN TIMES BETTER */
/* THAN MACHINE PRECISION. */

	r__1 = r1mach_(&c__3) * .1f;
	ngcs = inits_(gcs, &c__23, &r__1);

	gamlim_(&xmin, &xmax);
	dxrel = sqrt(r1mach_(&c__4));

/* --------------------------------------------------------------------- */
/* FINISH INITIALIZATION.  START EVALUATING GAMMA(X). */

    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 10.f) {
	goto L50;
    }

/* COMPUTE GAMMA(X) FOR ABS(X) .LE. 10.0.  REDUCE INTERVAL AND */
/* FIND GAMMA(1+Y) FOR 0. .LE. Y .LT. 1. FIRST OF ALL. */

    n = *x;
    if (*x < 0.f) {
	--n;
    }
    y = *x - n;
    --n;
    r__1 = y * 2.f - 1.f;
    ret_val = csevl_(&r__1, gcs, &ngcs) + .9375f;
    if (n == 0) {
	return ret_val;
    }

    if (n > 0) {
	goto L30;
    }

/* COMPUTE GAMMA(X) FOR X .LT. 1. */

    n = -n;
    if (*x == 0.f) {
	xermsg_("SLATEC", "GAMMA", "X IS 0", &c__4, &c__2, (ftnlen)6, (ftnlen)
		5, (ftnlen)6);
    }
    if (*x < 0.f && *x + n - 2 == 0.f) {
	xermsg_("SLATEC", "GAMMA", "X IS A NEGATIVE INTEGER", &c__4, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)23);
    }
    r__2 = *x - .5f;
    if (*x < -.5f && (r__1 = (*x - r_int(&r__2)) / *x, dabs(r__1)) < dxrel) {
	xermsg_("SLATEC", "GAMMA", "ANSWER LT HALF PRECISION BECAUSE X TOO N"
		"EAR NEGATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)5, (
		ftnlen)60);
    }

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val /= *x + i__ - 1;
/* L20: */
    }
    return ret_val;

/* GAMMA(X) FOR X .GE. 2. */

L30:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val = (y + i__) * ret_val;
/* L40: */
    }
    return ret_val;

/* COMPUTE GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X). */

L50:
    if (*x > xmax) {
	xermsg_("SLATEC", "GAMMA", "X SO BIG GAMMA OVERFLOWS", &c__3, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)24);
    }

    ret_val = 0.f;
    if (*x < xmin) {
	xermsg_("SLATEC", "GAMMA", "X SO SMALL GAMMA UNDERFLOWS", &c__2, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)27);
    }
    if (*x < xmin) {
	return ret_val;
    }

    ret_val = exp((y - .5f) * log(y) - y + sq2pil + r9lgmc_(&y));
    if (*x > 0.f) {
	return ret_val;
    }

    r__2 = *x - .5f;
    if ((r__1 = (*x - r_int(&r__2)) / *x, dabs(r__1)) < dxrel) {
	xermsg_("SLATEC", "GAMMA", "ANSWER LT HALF PRECISION, X TOO NEAR NEG"
		"ATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)
		53);
    }

    sinpiy = sin(pi * y);
    if (sinpiy == 0.f) {
	xermsg_("SLATEC", "GAMMA", "X IS A NEGATIVE INTEGER", &c__4, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)23);
    }

    ret_val = -pi / (y * sinpiy * ret_val);

    return ret_val;
} /* gamma_ */

