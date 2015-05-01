/* besi0e.f -- translated by f2c (version 12.02.01).
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
static integer c__12 = 12;
static integer c__21 = 21;
static integer c__22 = 22;

/* DECK BESI0E */
doublereal besi0e_(real *x)
{
    /* Initialized data */

    static real bi0cs[12] = { -.07660547252839144951f,1.92733795399380827f,
	    .2282644586920301339f,.01304891466707290428f,
	    4.3442709008164874e-4f,9.42265768600193e-6f,1.4340062895106e-7f,
	    1.61384906966e-9f,1.396650044e-11f,9.579451e-14f,5.3339e-16f,
	    2.45e-18f };
    static real ai0cs[21] = { .07575994494023796f,.00759138081082334f,
	    4.1531313389237e-4f,1.070076463439e-5f,-7.90117997921e-6f,
	    -7.8261435014e-7f,2.7838499429e-7f,8.2524726e-9f,-1.204463945e-8f,
	    1.55964859e-9f,2.2925563e-10f,-1.1916228e-10f,1.757854e-11f,
	    1.12822e-12f,-1.14684e-12f,2.7155e-13f,-2.415e-14f,-6.08e-15f,
	    3.14e-15f,-7.1e-16f,7e-17f };
    static real ai02cs[22] = { .05449041101410882f,.00336911647825569f,
	    6.889758346918e-5f,2.89137052082e-6f,2.0489185893e-7f,
	    2.266668991e-8f,3.39623203e-9f,4.9406022e-10f,1.188914e-11f,
	    -3.149915e-11f,-1.32158e-11f,-1.79419e-12f,7.1801e-13f,
	    3.8529e-13f,1.539e-14f,-4.151e-14f,-9.54e-15f,3.82e-15f,1.76e-15f,
	    -3.4e-16f,-2.7e-16f,3e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y;
    static integer nti0;
    static real xsml;
    static integer ntai0, ntai02;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  BESI0E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the first kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI0E-S, DBSI0E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB, */
/*             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION, */
/*             ORDER ZERO, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI0E(X) calculates the exponentially scaled modified (hyperbolic) */
/* Bessel function of the first kind of order zero for real argument X; */
/* i.e., EXP(-ABS(X))*I0(X). */


/* Series for BI0        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.46E-18 */
/*                                         log weighted error  17.61 */
/*                               significant figures required  17.90 */
/*                                    decimal places required  18.15 */


/* Series for AI0        on the interval  1.25000D-01 to  3.33333D-01 */
/*                                        with weighted error   7.87E-17 */
/*                                         log weighted error  16.10 */
/*                               significant figures required  14.69 */
/*                                    decimal places required  16.76 */


/* Series for AI02       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   3.79E-17 */
/*                                         log weighted error  16.42 */
/*                               significant figures required  14.86 */
/*                                    decimal places required  17.09 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890313  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  BESI0E */
/* ***FIRST EXECUTABLE STATEMENT  BESI0E */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nti0 = inits_(bi0cs, &c__12, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntai0 = inits_(ai0cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntai02 = inits_(ai02cs, &c__22, &r__1);
	xsml = sqrt(r1mach_(&c__3) * 4.5f);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 3.f) {
	goto L20;
    }

    ret_val = 1.f - *x;
    if (y > xsml) {
	r__1 = y * y / 4.5f - 1.f;
	ret_val = exp(-y) * (csevl_(&r__1, bi0cs, &nti0) + 2.75f);
    }
    return ret_val;

L20:
    if (y <= 8.f) {
	r__1 = (48.f / y - 11.f) / 5.f;
	ret_val = (csevl_(&r__1, ai0cs, &ntai0) + .375f) / sqrt(y);
    }
    if (y > 8.f) {
	r__1 = 16.f / y - 1.f;
	ret_val = (csevl_(&r__1, ai02cs, &ntai02) + .375f) / sqrt(y);
    }

    return ret_val;
} /* besi0e_ */

