/* besj0.f -- translated by f2c (version 12.02.01).
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
static integer c__13 = 13;
static integer c__21 = 21;
static integer c__24 = 24;
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESJ0 */
doublereal besj0_(real *x)
{
    /* Initialized data */

    static real bj0cs[13] = { .100254161968939137f,-.665223007764405132f,
	    .248983703498281314f,-.0332527231700357697f,.0023114179304694015f,
	    -9.9112774199508e-5f,2.8916708643998e-6f,-6.1210858663e-8f,
	    9.838650793e-10f,-1.24235515e-11f,1.265433e-13f,-1.0619e-15f,
	    7.4e-18f };
    static real bm0cs[21] = { .09284961637381644f,-.00142987707403484f,
	    2.830579271257e-5f,-1.43300611424e-6f,1.2028628046e-7f,
	    -1.397113013e-8f,2.04076188e-9f,-3.5399669e-10f,7.024759e-11f,
	    -1.554107e-11f,3.76226e-12f,-9.8282e-13f,2.7408e-13f,-8.091e-14f,
	    2.511e-14f,-8.14e-15f,2.75e-15f,-9.6e-16f,3.4e-16f,-1.2e-16f,
	    4e-17f };
    static real bth0cs[24] = { -.24639163774300119f,.001737098307508963f,
	    -6.2183633402968e-5f,4.368050165742e-6f,-4.56093019869e-7f,
	    6.2197400101e-8f,-1.0300442889e-8f,1.979526776e-9f,
	    -4.28198396e-10f,1.0203584e-10f,-2.6363898e-11f,7.297935e-12f,
	    -2.144188e-12f,6.63693e-13f,-2.15126e-13f,7.2659e-14f,
	    -2.5465e-14f,9.229e-15f,-3.448e-15f,1.325e-15f,-5.22e-16f,
	    2.1e-16f,-8.7e-17f,3.6e-17f };
    static real pi4 = .78539816339744831f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y, z__;
    static integer ntj0, ntm0;
    static real ampl, xmax, xsml;
    static integer ntth0;
    static real theta;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESJ0 */
/* ***PURPOSE  Compute the Bessel function of the first kind of order */
/*            zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10A1 */
/* ***TYPE      SINGLE PRECISION (BESJ0-S, DBESJ0-D) */
/* ***KEYWORDS  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESJ0(X) calculates the Bessel function of the first kind of */
/* order zero for real argument X. */

/* Series for BJ0        on the interval  0.          to  1.60000D+01 */
/*                                        with weighted error   7.47E-18 */
/*                                         log weighted error  17.13 */
/*                               significant figures required  16.98 */
/*                                    decimal places required  17.68 */

/* Series for BM0        on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   4.98E-17 */
/*                                         log weighted error  16.30 */
/*                               significant figures required  14.97 */
/*                                    decimal places required  16.96 */

/* Series for BTH0       on the interval  0.          to  6.25000D-02 */
/*                                        with weighted error   3.67E-17 */
/*                                         log weighted error  16.44 */
/*                               significant figures required  15.53 */
/*                                    decimal places required  17.13 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890210  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESJ0 */
/* ***FIRST EXECUTABLE STATEMENT  BESJ0 */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntj0 = inits_(bj0cs, &c__13, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntm0 = inits_(bm0cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntth0 = inits_(bth0cs, &c__24, &r__1);

	xsml = sqrt(r1mach_(&c__3) * 8.f);
	xmax = 1.f / r1mach_(&c__4);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 4.f) {
	goto L20;
    }

    ret_val = 1.f;
    if (y > xsml) {
	r__1 = y * .125f * y - 1.f;
	ret_val = csevl_(&r__1, bj0cs, &ntj0);
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESJ0", "NO PRECISION BECAUSE ABS(X) IS TOO BIG", &
		c__1, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)38);
    }

/* Computing 2nd power */
    r__1 = y;
    z__ = 32.f / (r__1 * r__1) - 1.f;
    ampl = (csevl_(&z__, bm0cs, &ntm0) + .75f) / sqrt(y);
    theta = y - pi4 + csevl_(&z__, bth0cs, &ntth0) / y;
    ret_val = ampl * cos(theta);

    return ret_val;
} /* besj0_ */

