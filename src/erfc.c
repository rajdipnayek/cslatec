/* erfc.f -- translated by f2c (version 12.02.01).
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
static integer c__24 = 24;
static integer c__23 = 23;
static integer c__1 = 1;

/* DECK ERFC */
doublereal erfc_(real *x)
{
    /* Initialized data */

    static real erfcs[13] = { -.049046121234691808f,-.14226120510371364f,
	    .010035582187599796f,-5.76876469976748e-4f,2.7419931252196e-5f,
	    -1.104317550734e-6f,3.848875542e-8f,-1.180858253e-9f,
	    3.2334215e-11f,-7.99101e-13f,1.799e-14f,-3.71e-16f,7e-18f };
    static real erc2cs[23] = { -.069601346602309501f,-.041101339362620893f,
	    .003914495866689626f,-4.90639565054897e-4f,7.157479001377e-5f,
	    -1.1530716341312e-5f,1.994670590201e-6f,-3.64266647159e-7f,
	    6.94437261e-8f,-1.3712209021e-8f,2.788389661e-9f,-5.81416472e-10f,
	    1.23892049e-10f,-2.6906391e-11f,5.942614e-12f,-1.332386e-12f,
	    3.02804e-13f,-6.9666e-14f,1.6208e-14f,-3.809e-15f,9.04e-16f,
	    -2.16e-16f,5.2e-17f };
    static real erfccs[24] = { .0715179310202925f,-.026532434337606719f,
	    .001711153977920853f,-1.63751663458512e-4f,1.9871293500549e-5f,
	    -2.843712412769e-6f,4.60616130901e-7f,-8.2277530261e-8f,
	    1.5921418724e-8f,-3.295071356e-9f,7.22343973e-10f,
	    -1.66485584e-10f,4.0103931e-11f,-1.0048164e-11f,2.608272e-12f,
	    -6.99105e-13f,1.92946e-13f,-5.4704e-14f,1.5901e-14f,-4.729e-15f,
	    1.432e-15f,-4.39e-16f,1.38e-16f,-4.8e-17f };
    static real sqrtpi = 1.772453850905516f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y, eta, xmax, xsml;
    extern doublereal csevl_(real *, real *, integer *);
    static integer nterf;
    extern integer inits_(real *, integer *, real *);
    static real sqeps, txmax;
    extern doublereal r1mach_(integer *);
    static integer nterc2, nterfc;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  ERFC */
/* ***PURPOSE  Compute the complementary error function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C8A, L5A1E */
/* ***TYPE      SINGLE PRECISION (ERFC-S, DERFC-D) */
/* ***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ERFC(X) calculates the single precision complementary error */
/* function for single precision argument X. */

/* Series for ERF        on the interval  0.          to  1.00000D+00 */
/*                                        with weighted error   7.10E-18 */
/*                                         log weighted error  17.15 */
/*                               significant figures required  16.31 */
/*                                    decimal places required  17.71 */

/* Series for ERFC       on the interval  0.          to  2.50000D-01 */
/*                                        with weighted error   4.81E-17 */
/*                                         log weighted error  16.32 */
/*                        approx significant figures required  15.0 */


/* Series for ERC2       on the interval  2.50000D-01 to  1.00000D+00 */
/*                                        with weighted error   5.22E-17 */
/*                                         log weighted error  16.28 */
/*                        approx significant figures required  15.0 */
/*                                    decimal places required  16.96 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  ERFC */
/* ***FIRST EXECUTABLE STATEMENT  ERFC */
    if (first) {
	eta = r1mach_(&c__3) * .1f;
	nterf = inits_(erfcs, &c__13, &eta);
	nterfc = inits_(erfccs, &c__24, &eta);
	nterc2 = inits_(erc2cs, &c__23, &eta);

	xsml = -sqrt(-log(sqrtpi * r1mach_(&c__3)));
	txmax = sqrt(-log(sqrtpi * r1mach_(&c__1)));
	xmax = txmax - log(txmax) * .5f / txmax - .01f;
	sqeps = sqrt(r1mach_(&c__3) * 2.f);
    }
    first = FALSE_;

    if (*x > xsml) {
	goto L20;
    }

/* ERFC(X) = 1.0 - ERF(X) FOR X .LT. XSML */

    ret_val = 2.f;
    return ret_val;

L20:
    if (*x > xmax) {
	goto L40;
    }
    y = dabs(*x);
    if (y > 1.f) {
	goto L30;
    }

/* ERFC(X) = 1.0 - ERF(X) FOR -1. .LE. X .LE. 1. */

    if (y < sqeps) {
	ret_val = 1.f - *x * 2.f / sqrtpi;
    }
    if (y >= sqeps) {
	r__1 = *x * 2.f * *x - 1.f;
	ret_val = 1.f - *x * (csevl_(&r__1, erfcs, &nterf) + 1.f);
    }
    return ret_val;

/* ERFC(X) = 1.0 - ERF(X) FOR 1. .LT. ABS(X) .LE. XMAX */

L30:
    y *= y;
    if (y <= 4.f) {
	r__1 = (8.f / y - 5.f) / 3.f;
	ret_val = exp(-y) / dabs(*x) * (csevl_(&r__1, erc2cs, &nterc2) + .5f);
    }
    if (y > 4.f) {
	r__1 = 8.f / y - 1.f;
	ret_val = exp(-y) / dabs(*x) * (csevl_(&r__1, erfccs, &nterfc) + .5f);
    }
    if (*x < 0.f) {
	ret_val = 2.f - ret_val;
    }
    return ret_val;

L40:
    xermsg_("SLATEC", "ERFC", "X SO BIG ERFC UNDERFLOWS", &c__1, &c__1, (
	    ftnlen)6, (ftnlen)4, (ftnlen)24);
    ret_val = 0.f;
    return ret_val;

} /* erfc_ */

