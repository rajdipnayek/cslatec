/* binom.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;

/* DECK BINOM */
doublereal binom_(integer *n, integer *m)
{
    /* Initialized data */

    static real sq2pil = .91893853320467274f;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1;

    /* Local variables */
    static integer i__, k;
    static real xk, xn, xnk, corr;
    extern doublereal r1mach_(integer *), r9lgmc_(real *), alnrel_(real *);
    static real bilnmx, fintmx;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BINOM */
/* ***PURPOSE  Compute the binomial coefficients. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1 */
/* ***TYPE      SINGLE PRECISION (BINOM-S, DBINOM-D) */
/* ***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BINOM(N,M) calculates the binomial coefficient (N!)/((M!)*(N-M)!). */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALNREL, R1MACH, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BINOM */
/* ***FIRST EXECUTABLE STATEMENT  BINOM */
    if (first) {
	bilnmx = log(r1mach_(&c__2));
	fintmx = .9f / r1mach_(&c__3);
    }
    first = FALSE_;

    if (*n < 0 || *m < 0) {
	xermsg_("SLATEC", "BINOM", "N OR M LT ZERO", &c__1, &c__2, (ftnlen)6, 
		(ftnlen)5, (ftnlen)14);
    }
    if (*n < *m) {
	xermsg_("SLATEC", "BINOM", "N LT M", &c__2, &c__2, (ftnlen)6, (ftnlen)
		5, (ftnlen)6);
    }

/* Computing MIN */
    i__1 = *m, i__2 = *n - *m;
    k = min(i__1,i__2);
    if (k > 20) {
	goto L30;
    }
    if (k * log((real) max(*n,1)) > bilnmx) {
	goto L30;
    }

    ret_val = 1.f;
    if (k == 0) {
	return ret_val;
    }

    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val = ret_val * (real) (*n - i__ + 1) / i__;
/* L20: */
    }

    if (ret_val < fintmx) {
	r__1 = ret_val + .5f;
	ret_val = r_int(&r__1);
    }
    return ret_val;

/* IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM */
L30:
    if (k < 9) {
	xermsg_("SLATEC", "BINOM", "RESULT OVERFLOWS BECAUSE N AND/OR M TOO "
		"BIG", &c__3, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)43);
    }

    xn = (real) (*n + 1);
    xk = (real) (k + 1);
    xnk = (real) (*n - k + 1);

    corr = r9lgmc_(&xn) - r9lgmc_(&xk) - r9lgmc_(&xnk);
    r__1 = -(xk - 1.f) / xn;
    ret_val = xk * log(xnk / xk) - xn * alnrel_(&r__1) - log(xn * xnk / xk) * 
	    .5f + 1.f - sq2pil + corr;

    if (ret_val > bilnmx) {
	xermsg_("SLATEC", "BINOM", "RESULT OVERFLOWS BECAUSE N AND/OR M TOO "
		"BIG", &c__3, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)43);
    }

    ret_val = exp(ret_val);
    if (ret_val < fintmx) {
	r__1 = ret_val + .5f;
	ret_val = r_int(&r__1);
    }

    return ret_val;
} /* binom_ */

