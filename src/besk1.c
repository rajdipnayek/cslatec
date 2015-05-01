/* besk1.f -- translated by f2c (version 12.02.01).
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
static integer c__11 = 11;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESK1 */
doublereal besk1_(real *x)
{
    /* Initialized data */

    static real bk1cs[11] = { .0253002273389477705f,-.353155960776544876f,
	    -.122611180822657148f,-.0069757238596398643f,
	    -1.730288957513052e-4f,-2.4334061415659e-6f,-2.21338763073e-8f,
	    -1.411488392e-10f,-6.666901e-13f,-2.4274e-15f,-7e-18f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real y;
    static integer ntk1;
    static real xmin, xmax, xsml;
    extern doublereal besi1_(real *), csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    static real xmaxt;
    extern doublereal besk1e_(real *), r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESK1 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            third kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK1-S, DBESK1-D) */
/* ***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK1(X) computes the modified (hyperbolic) Bessel function of third */
/* kind of order one for real argument X, where X .GT. 0. */

/* Series for BK1        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   7.02E-18 */
/*                                         log weighted error  17.15 */
/*                               significant figures required  16.73 */
/*                                    decimal places required  17.67 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI1, BESK1E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK1 */
/* ***FIRST EXECUTABLE STATEMENT  BESK1 */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntk1 = inits_(bk1cs, &c__11, &r__1);
/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	xmin = exp(dmax(r__1,r__2) + .01f);
	xsml = sqrt(r1mach_(&c__3) * 4.f);
	xmaxt = -log(r1mach_(&c__1));
	xmax = xmaxt - xmaxt * .5f * log(xmaxt) / (xmaxt + .5f);
    }
    first = FALSE_;

    if (*x <= 0.f) {
	xermsg_("SLATEC", "BESK1", "X IS ZERO OR NEGATIVE", &c__2, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)21);
    }
    if (*x > 2.f) {
	goto L20;
    }

    if (*x < xmin) {
	xermsg_("SLATEC", "BESK1", "X SO SMALL K1 OVERFLOWS", &c__3, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)23);
    }
    y = 0.f;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * .5f - 1.f;
    ret_val = log(*x * .5f) * besi1_(x) + (csevl_(&r__1, bk1cs, &ntk1) + .75f)
	     / *x;
    return ret_val;

L20:
    ret_val = 0.f;
    if (*x > xmax) {
	xermsg_("SLATEC", "BESK1", "X SO BIG K1 UNDERFLOWS", &c__1, &c__1, (
		ftnlen)6, (ftnlen)5, (ftnlen)22);
    }
    if (*x > xmax) {
	return ret_val;
    }

    ret_val = exp(-(*x)) * besk1e_(x);

    return ret_val;
} /* besk1_ */

