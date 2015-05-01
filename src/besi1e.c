/* besi1e.f -- translated by f2c (version 12.02.01).
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
static integer c__21 = 21;
static integer c__22 = 22;
static integer c__1 = 1;

/* DECK BESI1E */
doublereal besi1e_(real *x)
{
    /* Initialized data */

    static real bi1cs[11] = { -.001971713261099859f,.40734887667546481f,
	    .034838994299959456f,.001545394556300123f,4.1888521098377e-5f,
	    7.64902676483e-7f,1.0042493924e-8f,9.9322077e-11f,7.6638e-13f,
	    4.741e-15f,2.4e-17f };
    static real ai1cs[21] = { -.02846744181881479f,-.01922953231443221f,
	    -6.1151858579437e-4f,-2.06997125335e-5f,8.58561914581e-6f,
	    1.04949824671e-6f,-2.9183389184e-7f,-1.559378146e-8f,
	    1.318012367e-8f,-1.44842341e-9f,-2.9085122e-10f,1.2663889e-10f,
	    -1.664947e-11f,-1.66665e-12f,1.2426e-12f,-2.7315e-13f,2.023e-14f,
	    7.3e-15f,-3.33e-15f,7.1e-16f,-6e-17f };
    static real ai12cs[22] = { .02857623501828014f,-.00976109749136147f,
	    -1.1058893876263e-4f,-3.88256480887e-6f,-2.5122362377e-7f,
	    -2.631468847e-8f,-3.83538039e-9f,-5.5897433e-10f,-1.897495e-11f,
	    3.252602e-11f,1.41258e-11f,2.03564e-12f,-7.1985e-13f,-4.0836e-13f,
	    -2.101e-14f,4.273e-14f,1.041e-14f,-3.82e-15f,-1.86e-15f,3.3e-16f,
	    2.8e-16f,-3e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y;
    static integer nti1;
    static real xmin, xsml;
    static integer ntai1, ntai12;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESI1E */
/* ***PURPOSE  Compute the exponentially scaled modified (hyperbolic) */
/*            Bessel function of the first kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI1E-S, DBSI1E-D) */
/* ***KEYWORDS  EXPONENTIALLY SCALED, FIRST KIND, FNLIB, */
/*             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION, */
/*             ORDER ONE, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI1E(X) calculates the exponentially scaled modified (hyperbolic) */
/* Bessel function of the first kind of order one for real argument X; */
/* i.e., EXP(-ABS(X))*I1(X). */

/* Series for BI1        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.40E-17 */
/*                                         log weighted error  16.62 */
/*                               significant figures required  16.23 */
/*                                    decimal places required  17.14 */

/* Series for AI1        on the interval  1.25000D-01 to  3.33333D-01 */
/*                                        with weighted error   6.98E-17 */
/*                                         log weighted error  16.16 */
/*                               significant figures required  14.53 */
/*                                    decimal places required  16.82 */

/* Series for AI12       on the interval  0.          to  1.25000D-01 */
/*                                        with weighted error   3.55E-17 */
/*                                         log weighted error  16.45 */
/*                               significant figures required  14.69 */
/*                                    decimal places required  17.12 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890210  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  BESI1E */
/* ***FIRST EXECUTABLE STATEMENT  BESI1E */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nti1 = inits_(bi1cs, &c__11, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntai1 = inits_(ai1cs, &c__21, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	ntai12 = inits_(ai12cs, &c__22, &r__1);

	xmin = r1mach_(&c__1) * 2.f;
	xsml = sqrt(r1mach_(&c__3) * 4.5f);
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 3.f) {
	goto L20;
    }

    ret_val = 0.f;
    if (y == 0.f) {
	return ret_val;
    }

    if (y <= xmin) {
	xermsg_("SLATEC", "BESI1E", "ABS(X) SO SMALL I1 UNDERFLOWS", &c__1, &
		c__1, (ftnlen)6, (ftnlen)6, (ftnlen)29);
    }
    if (y > xmin) {
	ret_val = *x * .5f;
    }
    if (y > xsml) {
	r__1 = y * y / 4.5f - 1.f;
	ret_val = *x * (csevl_(&r__1, bi1cs, &nti1) + .875f);
    }
    ret_val = exp(-y) * ret_val;
    return ret_val;

L20:
    if (y <= 8.f) {
	r__1 = (48.f / y - 11.f) / 5.f;
	ret_val = (csevl_(&r__1, ai1cs, &ntai1) + .375f) / sqrt(y);
    }
    if (y > 8.f) {
	r__1 = 16.f / y - 1.f;
	ret_val = (csevl_(&r__1, ai12cs, &ntai12) + .375f) / sqrt(y);
    }
    ret_val = r_sign(&ret_val, x);

    return ret_val;
} /* besi1e_ */

