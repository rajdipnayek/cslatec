/* asinh.f -- translated by f2c (version 12.02.01).
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
static integer c__20 = 20;

/* DECK ASINH */
doublereal asinh_(real *x)
{
    /* Initialized data */

    static real aln2 = .69314718055994530942f;
    static real asnhcs[20] = { -.12820039911738186f,-.058811761189951768f,
	    .004727465432212481f,-4.93836316265361e-4f,5.8506207058557e-5f,
	    -7.466998328931e-6f,1.001169358355e-6f,-1.39035438587e-7f,
	    1.9823169483e-8f,-2.884746841e-9f,4.26729654e-10f,-6.3976084e-11f,
	    9.699168e-12f,-1.484427e-12f,2.29037e-13f,-3.5588e-14f,5.563e-15f,
	    -8.74e-16f,1.38e-16f,-2.1e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y, xmax;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    static real sqeps;
    extern doublereal r1mach_(integer *);
    static integer nterms;

/* ***BEGIN PROLOGUE  ASINH */
/* ***PURPOSE  Compute the arc hyperbolic sine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      SINGLE PRECISION (ASINH-S, DASINH-D, CASINH-C) */
/* ***KEYWORDS  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB, */
/*             INVERSE HYPERBOLIC SINE */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ASINH(X) computes the arc hyperbolic sine of X. */

/* Series for ASNH       on the interval  0.          to  1.00000D+00 */
/*                                        with weighted error   2.19E-17 */
/*                                         log weighted error  16.66 */
/*                               significant figures required  15.60 */
/*                                    decimal places required  17.31 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  ASINH */
/* ***FIRST EXECUTABLE STATEMENT  ASINH */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nterms = inits_(asnhcs, &c__20, &r__1);
	sqeps = sqrt(r1mach_(&c__3));
	xmax = 1.f / sqeps;
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 1.f) {
	goto L20;
    }

    ret_val = *x;
    if (y > sqeps) {
	r__1 = *x * 2.f * *x - 1.f;
	ret_val = *x * (csevl_(&r__1, asnhcs, &nterms) + 1.f);
    }
    return ret_val;

L20:
    if (y < xmax) {
/* Computing 2nd power */
	r__1 = y;
	ret_val = log(y + sqrt(r__1 * r__1 + 1.f));
    }
    if (y >= xmax) {
	ret_val = aln2 + log(y);
    }
    ret_val = r_sign(&ret_val, x);

    return ret_val;
} /* asinh_ */

