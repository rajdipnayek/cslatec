/* r9aimp.f -- translated by f2c (version 12.02.01).
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
static integer c__40 = 40;
static integer c__36 = 36;
static integer c__33 = 33;
static integer c__32 = 32;
static doublereal c_b8 = .3333;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK R9AIMP */
/* Subroutine */ int r9aimp_(real *x, real *ampl, real *theta)
{
    /* Initialized data */

    static real am21cs[40] = { .0065809191761485f,.0023675984685722f,
	    1.324741670371e-4f,1.57600904043e-5f,2.7529702663e-6f,
	    6.102679017e-7f,1.595088468e-7f,4.71033947e-8f,1.52933871e-8f,
	    5.3590722e-9f,2.000091e-9f,7.872292e-10f,3.243103e-10f,
	    1.390106e-10f,6.17011e-11f,2.82491e-11f,1.32979e-11f,6.4188e-12f,
	    3.1697e-12f,1.5981e-12f,8.213e-13f,4.296e-13f,2.284e-13f,
	    1.232e-13f,6.75e-14f,3.74e-14f,2.1e-14f,1.19e-14f,6.8e-15f,
	    3.9e-15f,2.3e-15f,1.3e-15f,8e-16f,5e-16f,3e-16f,1e-16f,1e-16f,0.f,
	    0.f,0.f };
    static real ath1cs[36] = { -.07125837815669365f,-.00590471979831451f,
	    -1.2114544069499e-4f,-9.8860854227e-6f,-1.38084097352e-6f,
	    -2.6142640172e-7f,-6.050432589e-8f,-1.618436223e-8f,
	    -4.83464911e-9f,-1.57655272e-9f,-5.5231518e-10f,-2.0545441e-10f,
	    -8.043412e-11f,-3.291252e-11f,-1.399875e-11f,-6.16151e-12f,
	    -2.79614e-12f,-1.30428e-12f,-6.2373e-13f,-3.0512e-13f,
	    -1.5239e-13f,-7.758e-14f,-4.02e-14f,-2.117e-14f,-1.132e-14f,
	    -6.14e-15f,-3.37e-15f,-1.88e-15f,-1.05e-15f,-6e-16f,-3.4e-16f,
	    -2e-16f,-1.1e-16f,-7e-17f,-4e-17f,-2e-17f };
    static real am22cs[33] = { -.01562844480625341f,.00778336445239681f,
	    8.6705777047718e-4f,1.5696627315611e-4f,3.563962571432e-5f,
	    9.24598335425e-6f,2.6211016185e-6f,7.9188221651e-7f,
	    2.5104152792e-7f,8.265223206e-8f,2.805711662e-8f,9.7682109e-9f,
	    3.47407923e-9f,1.25828132e-9f,4.6298826e-10f,1.7272825e-10f,
	    6.523192e-11f,2.490471e-11f,9.60156e-12f,3.73448e-12f,
	    1.46417e-12f,5.7826e-13f,2.2991e-13f,9.197e-14f,3.7e-14f,
	    1.496e-14f,6.08e-15f,2.48e-15f,1.01e-15f,4.1e-16f,1.7e-16f,7e-17f,
	    2e-17f };
    static real ath2cs[32] = { .00440527345871877f,-.03042919452318455f,
	    -.00138565328377179f,-1.8044439089549e-4f,-3.380847108327e-5f,
	    -7.67818353522e-6f,-1.96783944371e-6f,-5.4837271158e-7f,
	    -1.6254615505e-7f,-5.053049981e-8f,-1.631580701e-8f,
	    -5.43420411e-9f,-1.85739855e-9f,-6.489512e-10f,-2.3105948e-10f,
	    -8.363282e-11f,-3.071196e-11f,-1.142367e-11f,-4.29811e-12f,
	    -1.63389e-12f,-6.2693e-13f,-2.426e-13f,-9.461e-14f,-3.716e-14f,
	    -1.469e-14f,-5.84e-15f,-2.33e-15f,-9.3e-16f,-3.7e-16f,-1.5e-16f,
	    -6e-17f,-2e-17f };
    static real pi4 = .78539816339744831f;
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal d__1;

    /* Local variables */
    static real z__, eta;
    static integer nam21, nam22;
    static real xsml;
    static integer nath1, nath2;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    static real sqrtx;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9AIMP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Evaluate the Airy modulus and phase. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      SINGLE PRECISION (R9AIMP-S, D9AIMP-D) */
/* ***KEYWORDS  AIRY FUNCTION, FNLIB, MODULUS, PHASE, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate the Airy modulus and phase for X .LE. -1.0 */

/* Series for AM21       on the interval -1.25000D-01 to  0. */
/*                                        with weighted error   2.89E-17 */
/*                                         log weighted error  16.54 */
/*                               significant figures required  14.15 */
/*                                    decimal places required  17.34 */

/* Series for ATH1       on the interval -1.25000D-01 to  0. */
/*                                        with weighted error   2.53E-17 */
/*                                         log weighted error  16.60 */
/*                               significant figures required  15.15 */
/*                                    decimal places required  17.38 */

/* Series for AM22       on the interval -1.00000D+00 to -1.25000D-01 */
/*                                        with weighted error   2.99E-17 */
/*                                         log weighted error  16.52 */
/*                               significant figures required  14.57 */
/*                                    decimal places required  17.28 */

/* Series for ATH2       on the interval -1.00000D+00 to -1.25000D-01 */
/*                                        with weighted error   2.57E-17 */
/*                                         log weighted error  16.59 */
/*                               significant figures required  15.07 */
/*                                    decimal places required  17.34 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9AIMP */
/* ***FIRST EXECUTABLE STATEMENT  R9AIMP */
    if (first) {
	eta = r1mach_(&c__3) * .1f;
	nam21 = inits_(am21cs, &c__40, &eta);
	nath1 = inits_(ath1cs, &c__36, &eta);
	nam22 = inits_(am22cs, &c__33, &eta);
	nath2 = inits_(ath2cs, &c__32, &eta);

	d__1 = (doublereal) r1mach_(&c__3);
	xsml = -1.f / pow_dd(&d__1, &c_b8);
    }
    first = FALSE_;

    if (*x >= -2.f) {
	goto L20;
    }
    z__ = 1.f;
    if (*x > xsml) {
/* Computing 3rd power */
	r__1 = *x;
	z__ = 16.f / (r__1 * (r__1 * r__1)) + 1.f;
    }
    *ampl = csevl_(&z__, am21cs, &nam21) + .3125f;
    *theta = csevl_(&z__, ath1cs, &nath1) - .625f;
    goto L30;

L20:
    if (*x > -1.f) {
	xermsg_("SLATEC", "R9AIMP", "X MUST BE LE -1.0", &c__1, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)17);
    }

/* Computing 3rd power */
    r__1 = *x;
    z__ = (16.f / (r__1 * (r__1 * r__1)) + 9.f) / 7.f;
    *ampl = csevl_(&z__, am22cs, &nam22) + .3125f;
    *theta = csevl_(&z__, ath2cs, &nath2) - .625f;

L30:
    sqrtx = sqrt(-(*x));
    *ampl = sqrt(*ampl / sqrtx);
    *theta = pi4 - *x * sqrtx * *theta;

    return 0;
} /* r9aimp_ */

