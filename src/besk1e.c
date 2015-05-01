/* besk1e.f -- translated by f2c (version 12.02.01).
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
static integer c__17 = 17;
static integer c__14 = 14;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESK1E */
doublereal besk1e_(real *x)
{
    /* Initialized data */

    static real bk1cs[11] = { .0253002273389477705f,-.353155960776544876f,
	    -.122611180822657148f,-.0069757238596398643f,
	    -1.730288957513052e-4f,-2.4334061415659e-6f,-2.21338763073e-8f,
	    -1.411488392e-10f,-6.666901e-13f,-2.4274e-15f,-7e-18f };
    static real ak1cs[17] = { .2744313406973883f,.07571989953199368f,
	    -.0014410515564754f,6.650116955125e-5f,-4.36998470952e-6f,
	    3.5402774997e-7f,-3.311163779e-8f,3.44597758e-9f,-3.8989323e-10f,
	    4.720819e-11f,-6.04783e-12f,8.1284e-13f,-1.1386e-13f,1.654e-14f,
	    -2.48e-15f,3.8e-16f,-6e-17f };
    static real ak12cs[14] = { .06379308343739001f,.02832887813049721f,
	    -2.4753706739052e-4f,5.7719724516e-6f,-2.0689392195e-7f,
	    9.73998344e-9f,-5.5853361e-10f,3.732996e-11f,-2.82505e-12f,
	    2.372e-13f,-2.176e-14f,2.15e-15f,-2.2e-16f,2e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real y;
    static integer ntk1;
    static real xmin, xsml;
    extern doublereal besi1_(real *);
    static integer ntak1, ntak12;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESK1E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the third kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK1E-S, DBSK1E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK1E(X) computes the exponentially scaled modified (hyperbolic) */
/* Bessel function of third kind of order one for real argument */
/* X .GT. 0.0, i.e., EXP(X)*K1(X). */

/* Series for BK1        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   7.02E-18 */
/*                                         log weighted error  17.15 */
/*                               significant figures required  16.73 */
/*                                    decimal places required  17.67 */

/* Series for AK1        on the interval  1.25000D-01 to  5.00000D-01 */
/*                                        with weighted error   6.06E-17 */
/*                                         log weighted error  16.22 */
/*                               significant figures required  15.41 */
/*                                    decimal places required  16.83 */

/* Series for AK12       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   2.58E-17 */
/*                                         log weighted error  16.59 */
/*                               significant figures required  15.22 */
/*                                    decimal places required  17.16 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI1, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK1E */
/* ***FIRST EXECUTABLE STATEMENT  BESK1E */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntk1 = inits_(bk1cs, &c__11, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntak1 = inits_(ak1cs, &c__17, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntak12 = inits_(ak12cs, &c__14, &r__1);

/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	xmin = exp(dmax(r__1,r__2) + .01f);
	xsml = sqrt(r1mach_(&c__3) * 4.f);
    }
    first = FALSE_;

    if (*x <= 0.f) {
	xermsg_("SLATEC", "BESK1E", "X IS ZERO OR NEGATIVE", &c__2, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x > 2.f) {
	goto L20;
    }

    if (*x < xmin) {
	xermsg_("SLATEC", "BESK1E", "X SO SMALL K1 OVERFLOWS", &c__3, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)23);
    }
    y = 0.f;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * .5f - 1.f;
    ret_val = exp(*x) * (log(*x * .5f) * besi1_(x) + (csevl_(&r__1, bk1cs, &
	    ntk1) + .75f) / *x);
    return ret_val;

L20:
    if (*x <= 8.f) {
	r__1 = (16.f / *x - 5.f) / 3.f;
	ret_val = (csevl_(&r__1, ak1cs, &ntak1) + 1.25f) / sqrt(*x);
    }
    if (*x > 8.f) {
	r__1 = 16.f / *x - 1.f;
	ret_val = (csevl_(&r__1, ak12cs, &ntak12) + 1.25f) / sqrt(*x);
    }

    return ret_val;
} /* besk1e_ */

