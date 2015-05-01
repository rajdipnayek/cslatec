/* dbinom.f -- translated by f2c (version 12.02.01).
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

/* DECK DBINOM */
doublereal dbinom_(integer *n, integer *m)
{
    /* Initialized data */

    static doublereal sq2pil = .91893853320467274178032973640562;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__, k;
    static doublereal xk, xn, xnk, corr;
    extern doublereal d1mach_(integer *), d9lgmc_(doublereal *), dlnrel_(
	    doublereal *);
    static doublereal bilnmx, fintmx;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBINOM */
/* ***PURPOSE  Compute the binomial coefficients. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1 */
/* ***TYPE      DOUBLE PRECISION (BINOM-S, DBINOM-D) */
/* ***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBINOM(N,M) calculates the double precision binomial coefficient */
/* for integer arguments N and M.  The result is (N!)/((M!)(N-M)!). */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9LGMC, DLNREL, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBINOM */
/* ***FIRST EXECUTABLE STATEMENT  DBINOM */
    if (first) {
	bilnmx = log(d1mach_(&c__2)) - 1e-4;
	fintmx = .9 / d1mach_(&c__3);
    }
    first = FALSE_;

    if (*n < 0 || *m < 0) {
	xermsg_("SLATEC", "DBINOM", "N OR M LT ZERO", &c__1, &c__2, (ftnlen)6,
		 (ftnlen)6, (ftnlen)14);
    }
    if (*n < *m) {
	xermsg_("SLATEC", "DBINOM", "N LT M", &c__2, &c__2, (ftnlen)6, (
		ftnlen)6, (ftnlen)6);
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

    ret_val = 1.;
    if (k == 0) {
	return ret_val;
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xn = (doublereal) (*n - i__ + 1);
	xk = (doublereal) i__;
	ret_val *= xn / xk;
/* L20: */
    }

    if (ret_val < fintmx) {
	d__1 = ret_val + .5;
	ret_val = d_int(&d__1);
    }
    return ret_val;

/* IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM */
L30:
    if (k < 9) {
	xermsg_("SLATEC", "DBINOM", "RESULT OVERFLOWS BECAUSE N AND/OR M TOO"
		" BIG", &c__3, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)43);
    }

    xn = (doublereal) (*n + 1);
    xk = (doublereal) (k + 1);
    xnk = (doublereal) (*n - k + 1);

    corr = d9lgmc_(&xn) - d9lgmc_(&xk) - d9lgmc_(&xnk);
    d__1 = -(xk - 1.) / xn;
    ret_val = xk * log(xnk / xk) - xn * dlnrel_(&d__1) - log(xn * xnk / xk) * 
	    .5 + 1. - sq2pil + corr;

    if (ret_val > bilnmx) {
	xermsg_("SLATEC", "DBINOM", "RESULT OVERFLOWS BECAUSE N AND/OR M TOO"
		" BIG", &c__3, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)43);
    }

    ret_val = exp(ret_val);
    if (ret_val < fintmx) {
	d__1 = ret_val + .5;
	ret_val = d_int(&d__1);
    }

    return ret_val;
} /* dbinom_ */

