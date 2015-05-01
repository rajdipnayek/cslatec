/* erf.f -- translated by f2c (version 12.02.01).
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
static real c_b7 = 1.f;

/* DECK ERF */
doublereal erf_(real *x)
{
    /* Initialized data */

    static real erfcs[13] = { -.049046121234691808f,-.14226120510371364f,
	    .010035582187599796f,-5.76876469976748e-4f,2.7419931252196e-5f,
	    -1.104317550734e-6f,3.848875542e-8f,-1.180858253e-9f,
	    3.2334215e-11f,-7.99101e-13f,1.799e-14f,-3.71e-16f,7e-18f };
    static real sqrtpi = 1.772453850905516f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real y;
    extern doublereal erfc_(real *);
    static real xbig;
    extern doublereal csevl_(real *, real *, integer *);
    static integer nterf;
    extern integer inits_(real *, integer *, real *);
    static real sqeps;
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  ERF */
/* ***PURPOSE  Compute the error function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C8A, L5A1E */
/* ***TYPE      SINGLE PRECISION (ERF-S, DERF-D) */
/* ***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ERF(X) calculates the single precision error function for */
/* single precision argument X. */

/* Series for ERF        on the interval  0.          to  1.00000D+00 */
/*                                        with weighted error   7.10E-18 */
/*                                         log weighted error  17.15 */
/*                               significant figures required  16.31 */
/*                                    decimal places required  17.71 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, ERFC, INITS, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/*   920618  Removed space from variable name.  (RWC, WRB) */
/* ***END PROLOGUE  ERF */
/* ***FIRST EXECUTABLE STATEMENT  ERF */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nterf = inits_(erfcs, &c__13, &r__1);
	xbig = sqrt(-log(sqrtpi * r1mach_(&c__3)));
	sqeps = sqrt(r1mach_(&c__3) * 2.f);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 1.f) {
	goto L20;
    }

/* ERF(X) = 1. - ERFC(X) FOR -1. .LE. X .LE. 1. */

    if (y <= sqeps) {
	ret_val = *x * 2.f / sqrtpi;
    }
    if (y > sqeps) {
/* Computing 2nd power */
	r__2 = *x;
	r__1 = r__2 * r__2 * 2.f - 1.f;
	ret_val = *x * (csevl_(&r__1, erfcs, &nterf) + 1.f);
    }
    return ret_val;

/* ERF(X) = 1. - ERFC(X) FOR  ABS(X) .GT. 1. */

L20:
    if (y <= xbig) {
	r__1 = 1.f - erfc_(&y);
	ret_val = r_sign(&r__1, x);
    }
    if (y > xbig) {
	ret_val = r_sign(&c_b7, x);
    }

    return ret_val;
} /* erf_ */

