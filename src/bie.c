/* bie.f -- translated by f2c (version 12.02.01).
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
static integer c__10 = 10;
static integer c__24 = 24;
static integer c__29 = 29;
static doublereal c_b9 = .3333;
static integer c__2 = 2;
static doublereal c_b11 = .6666;

/* DECK BIE */
doublereal bie_(real *x)
{
    /* Initialized data */

    static real bifcs[9] = { -.01673021647198664948f,.1025233583424944561f,
	    .00170830925073815165f,1.186254546774468e-5f,4.493290701779e-8f,
	    1.0698207143e-10f,1.7480643e-13f,2.081e-16f,1.8e-19f };
    static real bigcs[8] = { .02246622324857452f,.03736477545301955f,
	    4.4476218957212e-4f,2.47080756363e-6f,7.91913533e-9f,
	    1.649807e-11f,2.411e-14f,2e-17f };
    static real bif2cs[10] = { .09984572693816041f,.478624977863005538f,
	    .0251552119604330118f,5.820693885232645e-4f,7.4997659644377e-6f,
	    6.13460287034e-8f,3.462753885e-10f,1.428891e-12f,4.4962e-15f,
	    1.11e-17f };
    static real big2cs[10] = { .03330566214551434f,.161309215123197068f,
	    .0063190073096134286f,1.187904568162517e-4f,1.30453458862e-6f,
	    9.3741259955e-9f,4.74580188e-11f,1.783107e-13f,5.167e-16f,
	    1.1e-18f };
    static real bipcs[24] = { -.08322047477943447f,.01146118927371174f,
	    4.2896440718911e-4f,-1.490663937995e-4f,-1.307659726787e-5f,
	    6.3275983961e-6f,-4.2226696982e-7f,-1.9147186298e-7f,
	    6.453106284e-8f,-7.84485467e-9f,-9.6077216e-10f,7.0004713e-10f,
	    -1.7731789e-10f,2.272089e-11f,1.65404e-12f,-1.85171e-12f,
	    5.9576e-13f,-1.2194e-13f,1.334e-14f,1.72e-15f,-1.45e-15f,4.9e-16f,
	    -1.1e-16f,1e-17f };
    static real bip2cs[29] = { -.113596737585988679f,.0041381473947881595f,
	    1.353470622119332e-4f,1.04273166530153e-5f,1.3474954767849e-6f,
	    1.696537405438e-7f,-1.00965008656e-8f,-1.67291194937e-8f,
	    -4.5815364485e-9f,3.736681366e-10f,5.76693032e-10f,6.2181265e-11f,
	    -6.32941202e-11f,-1.49150479e-11f,7.8896213e-12f,2.4960513e-12f,
	    -1.2130075e-12f,-3.740493e-13f,2.237727e-13f,4.74902e-14f,
	    -4.52616e-14f,-3.0172e-15f,9.1058e-15f,-9.814e-16f,-1.6429e-15f,
	    5.533e-16f,2.175e-16f,-1.737e-16f,-1e-18f };
    static real atr = 8.7506905708484345f;
    static real btr = -2.093836321356054f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;
    doublereal d__1;

    /* Local variables */
    static real z__, xm, eta;
    static integer nbif, nbig, nbip;
    static real xbig;
    static integer nbif2, nbig2, nbip2;
    static real x3sml, theta;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    static real x32sml;
    extern doublereal r1mach_(integer *);
    static real sqrtx;
    extern /* Subroutine */ int r9aimp_(real *, real *, real *);

/* ***BEGIN PROLOGUE  BIE */
/* ***PURPOSE  Calculate the Bairy function for a negative argument and an */
/*            exponentially scaled Bairy function for a non-negative */
/*            argument. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      SINGLE PRECISION (BIE-S, DBIE-D) */
/* ***KEYWORDS  BAIRY FUNCTION, EXPONENTIALLY SCALED, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate BI(X) for X .LE. 0  and  BI(X)*EXP(ZETA)  where */
/* ZETA = 2/3 * X**(3/2)  for X .GE. 0.0 */

/* Series for BIF        on the interval -1.00000D+00 to  1.00000D+00 */
/*                                        with weighted error   1.88E-19 */
/*                                         log weighted error  18.72 */
/*                               significant figures required  17.74 */
/*                                    decimal places required  19.20 */

/* Series for BIG        on the interval -1.00000D+00 to  1.00000D+00 */
/*                                        with weighted error   2.61E-17 */
/*                                         log weighted error  16.58 */
/*                               significant figures required  15.17 */
/*                                    decimal places required  17.03 */

/* Series for BIF2       on the interval  1.00000D+00 to  8.00000D+00 */
/*                                        with weighted error   1.11E-17 */
/*                                         log weighted error  16.95 */
/*                        approx significant figures required  16.5 */
/*                                    decimal places required  17.45 */

/* Series for BIG2       on the interval  1.00000D+00 to  8.00000D+00 */
/*                                        with weighted error   1.19E-18 */
/*                                         log weighted error  17.92 */
/*                        approx significant figures required  17.2 */
/*                                    decimal places required  18.42 */

/* Series for BIP        on the interval  1.25000D-01 to  3.53553D-01 */
/*                                        with weighted error   1.91E-17 */
/*                                         log weighted error  16.72 */
/*                               significant figures required  15.35 */
/*                                    decimal places required  17.41 */

/* Series for BIP2       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   1.05E-18 */
/*                                         log weighted error  17.98 */
/*                               significant figures required  16.74 */
/*                                    decimal places required  18.71 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  BIE */
/* ***FIRST EXECUTABLE STATEMENT  BIE */
    if (first) {
	eta = r1mach_(&c__3) * .1f;
	nbif = inits_(bifcs, &c__9, &eta);
	nbig = inits_(bigcs, &c__8, &eta);
	nbif2 = inits_(bif2cs, &c__10, &eta);
	nbig2 = inits_(big2cs, &c__10, &eta);
	nbip = inits_(bipcs, &c__24, &eta);
	nbip2 = inits_(bip2cs, &c__29, &eta);

	d__1 = (doublereal) eta;
	x3sml = pow_dd(&d__1, &c_b9);
/* Computing 2nd power */
	r__1 = x3sml;
	x32sml = r__1 * r__1 * 1.3104f;
	d__1 = (doublereal) r1mach_(&c__2);
	xbig = pow_dd(&d__1, &c_b11);
    }
    first = FALSE_;

    if (*x >= -1.f) {
	goto L20;
    }
    r9aimp_(x, &xm, &theta);
    ret_val = xm * sin(theta);
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
    ret_val = csevl_(&z__, bifcs, &nbif) + .625f + *x * (csevl_(&z__, bigcs, &
	    nbig) + .4375f);
    if (*x > x32sml) {
	ret_val *= exp(*x * -2.f * sqrt(*x) / 3.f);
    }
    return ret_val;

L30:
    if (*x > 2.f) {
	goto L40;
    }
/* Computing 3rd power */
    r__1 = *x;
    z__ = (r__1 * (r__1 * r__1) * 2.f - 9.f) / 7.f;
    ret_val = exp(*x * -2.f * sqrt(*x) / 3.f) * (csevl_(&z__, bif2cs, &nbif2) 
	    + 1.125f + *x * (csevl_(&z__, big2cs, &nbig2) + .625f));
    return ret_val;

L40:
    if (*x > 4.f) {
	goto L50;
    }
    sqrtx = sqrt(*x);
    z__ = atr / (*x * sqrtx) + btr;
    ret_val = (csevl_(&z__, bipcs, &nbip) + .625f) / sqrt(sqrtx);
    return ret_val;

L50:
    sqrtx = sqrt(*x);
    z__ = -1.f;
    if (*x < xbig) {
	z__ = 16.f / (*x * sqrtx) - 1.f;
    }
    ret_val = (csevl_(&z__, bip2cs, &nbip2) + .625f) / sqrt(sqrtx);
    return ret_val;

} /* bie_ */

