/* besk0e.f -- translated by f2c (version 12.02.01).
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

/* DECK BESK0E */
doublereal besk0e_(real *x)
{
    /* Initialized data */

    static real bk0cs[11] = { -.03532739323390276872f,.3442898999246284869f,
	    .03597993651536150163f,.00126461541144692592f,
	    2.286212103119451e-5f,2.5347910790261e-7f,1.90451637722e-9f,
	    1.034969525e-11f,4.259816e-14f,1.3744e-16f,3.5e-19f };
    static real ak0cs[17] = { -.07643947903327941f,-.02235652605699819f,
	    7.7341811546938e-4f,-4.281006688886e-5f,3.08170017386e-6f,
	    -2.639367222e-7f,2.563713036e-8f,-2.74270554e-9f,3.1694296e-10f,
	    -3.902353e-11f,5.06804e-12f,-6.8895e-13f,9.744e-14f,-1.427e-14f,
	    2.15e-15f,-3.3e-16f,5e-17f };
    static real ak02cs[14] = { -.01201869826307592f,-.00917485269102569f,
	    1.444550931775e-4f,-4.01361417543e-6f,1.5678318108e-7f,
	    -7.77011043e-9f,4.6111825e-10f,-3.158592e-11f,2.43501e-12f,
	    -2.0743e-13f,1.925e-14f,-1.92e-15f,2e-16f,-2e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y;
    static integer ntk0;
    static real xsml;
    extern doublereal besi0_(real *);
    static integer ntak0, ntak02;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESK0E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the third kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK0E-S, DBSK0E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK0E(X) computes the exponentially scaled modified (hyperbolic) */
/* Bessel function of third kind of order zero for real argument */
/* X .GT. 0.0, i.e., EXP(X)*K0(X). */

/* Series for BK0        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   3.57E-19 */
/*                                         log weighted error  18.45 */
/*                               significant figures required  17.99 */
/*                                    decimal places required  18.97 */

/* Series for AK0        on the interval  1.25000D-01 to  5.00000D-01 */
/*                                        with weighted error   5.34E-17 */
/*                                         log weighted error  16.27 */
/*                               significant figures required  14.92 */
/*                                    decimal places required  16.89 */

/* Series for AK02       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   2.34E-17 */
/*                                         log weighted error  16.63 */
/*                               significant figures required  14.67 */
/*                                    decimal places required  17.20 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI0, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK0E */
/* ***FIRST EXECUTABLE STATEMENT  BESK0E */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntk0 = inits_(bk0cs, &c__11, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntak0 = inits_(ak0cs, &c__17, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntak02 = inits_(ak02cs, &c__14, &r__1);
	xsml = sqrt(r1mach_(&c__3) * 4.f);
    }
    first = FALSE_;

    if (*x <= 0.f) {
	xermsg_("SLATEC", "BESK0E", "X IS ZERO OR NEGATIVE", &c__2, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x > 2.f) {
	goto L20;
    }

    y = 0.f;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * .5f - 1.f;
    ret_val = exp(*x) * (-log(*x * .5f) * besi0_(x) - .25f + csevl_(&r__1, 
	    bk0cs, &ntk0));
    return ret_val;

L20:
    if (*x <= 8.f) {
	r__1 = (16.f / *x - 5.f) / 3.f;
	ret_val = (csevl_(&r__1, ak0cs, &ntak0) + 1.25f) / sqrt(*x);
    }
    if (*x > 8.f) {
	r__1 = 16.f / *x - 1.f;
	ret_val = (csevl_(&r__1, ak02cs, &ntak02) + 1.25f) / sqrt(*x);
    }

    return ret_val;
} /* besk0e_ */

