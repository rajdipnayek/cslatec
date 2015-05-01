/* aie.f -- translated by f2c (version 12.02.01).
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
static integer c__9 = 9;
static integer c__8 = 8;
static integer c__34 = 34;
static doublereal c_b6 = .3333;
static integer c__2 = 2;
static doublereal c_b8 = .6666;

/* DECK AIE */
doublereal aie_(real *x)
{
    /* Initialized data */

    static real aifcs[9] = { -.0379713584966699975f,.05919188853726363857f,
	    9.8629280577279975e-4f,6.84884381907656e-6f,2.594202596219e-8f,
	    6.176612774e-11f,1.0092454e-13f,1.2014e-16f,1e-19f };
    static real aigcs[8] = { .01815236558116127f,.02157256316601076f,
	    2.5678356987483e-4f,1.42652141197e-6f,4.57211492e-9f,9.52517e-12f,
	    1.392e-14f,1e-17f };
    static real aipcs[34] = { -.0187519297793868f,-.0091443848250055f,
	    9.010457337825e-4f,-1.394184127221e-4f,2.73815815785e-5f,
	    -6.2750421119e-6f,1.6064844184e-6f,-4.476392158e-7f,
	    1.334635874e-7f,-4.20735334e-8f,1.3902199e-8f,-4.7831848e-9f,
	    1.7047897e-9f,-6.268389e-10f,2.369824e-10f,-9.18641e-11f,
	    3.64278e-11f,-1.47475e-11f,6.0851e-12f,-2.5552e-12f,1.0906e-12f,
	    -4.725e-13f,2.076e-13f,-9.24e-14f,4.17e-14f,-1.9e-14f,8.7e-15f,
	    -4e-15f,1.9e-15f,-9e-16f,4e-16f,-2e-16f,1e-16f,-0.f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;
    doublereal d__1;

    /* Local variables */
    static real z__, xm, eta;
    static integer naif, naig, naip;
    static real xbig, x3sml, theta;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    static real x32sml;
    extern doublereal r1mach_(integer *);
    static real sqrtx;
    extern /* Subroutine */ int r9aimp_(real *, real *, real *);

/* ***BEGIN PROLOGUE  AIE */
/* ***PURPOSE  Calculate the Airy function for a negative argument and an */
/*            exponentially scaled Airy function for a non-negative */
/*            argument. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      SINGLE PRECISION (AIE-S, DAIE-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED AIRY FUNCTION, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* AIE(X) computes the exponentially scaled Airy function for */
/* non-negative X.  It evaluates AI(X) for X .LE. 0.0 and */
/* EXP(ZETA)*AI(X) for X .GE. 0.0 where ZETA = (2.0/3.0)*(X**1.5). */

/* Series for AIF        on the interval -1.00000D+00 to  1.00000D+00 */
/*                                        with weighted error   1.09E-19 */
/*                                         log weighted error  18.96 */
/*                               significant figures required  17.76 */
/*                                    decimal places required  19.44 */

/* Series for AIG        on the interval -1.00000D+00 to  1.00000D+00 */
/*                                        with weighted error   1.51E-17 */
/*                                         log weighted error  16.82 */
/*                               significant figures required  15.19 */
/*                                    decimal places required  17.27 */

/* Series for AIP        on the interval  0.          to  1.00000D+00 */
/*                                        with weighted error   5.10E-17 */
/*                                         log weighted error  16.29 */
/*                               significant figures required  14.41 */
/*                                    decimal places required  17.06 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  AIE */
/* ***FIRST EXECUTABLE STATEMENT  AIE */
    if (first) {
	eta = r1mach_(&c__3) * .1f;
	naif = inits_(aifcs, &c__9, &eta);
	naig = inits_(aigcs, &c__8, &eta);
	naip = inits_(aipcs, &c__34, &eta);

	d__1 = (doublereal) eta;
	x3sml = pow_dd(&d__1, &c_b6);
/* Computing 2nd power */
	r__1 = x3sml;
	x32sml = r__1 * r__1 * 1.3104f;
	d__1 = (doublereal) r1mach_(&c__2);
	xbig = pow_dd(&d__1, &c_b8);
    }
    first = FALSE_;

    if (*x >= -1.f) {
	goto L20;
    }
    r9aimp_(x, &xm, &theta);
    ret_val = xm * cos(theta);
    return ret_val;

L20:
    if (*x > 1.f) {
	goto L30;
    }
    z__ = 0.f;
    if (dabs(*x) > x3sml) {
/* Computing 3rd power */
	r__1 = *x;
	z__ = r__1 * (r__1 * r__1);
    }
    ret_val = csevl_(&z__, aifcs, &naif) - *x * (csevl_(&z__, aigcs, &naig) + 
	    .25f) + .375f;
    if (*x > x32sml) {
	ret_val *= exp(*x * 2.f * sqrt(*x) / 3.f);
    }
    return ret_val;

L30:
    sqrtx = sqrt(*x);
    z__ = -1.f;
    if (*x < xbig) {
	z__ = 2.f / (*x * sqrtx) - 1.f;
    }
    ret_val = (csevl_(&z__, aipcs, &naip) + .28125f) / sqrt(sqrtx);
    return ret_val;

} /* aie_ */

