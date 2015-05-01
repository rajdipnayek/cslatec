/* besi1.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BESI1 */
doublereal besi1_(real *x)
{
    /* Initialized data */

    static real bi1cs[11] = { -.001971713261099859f,.40734887667546481f,
	    .034838994299959456f,.001545394556300123f,4.1888521098377e-5f,
	    7.64902676483e-7f,1.0042493924e-8f,9.9322077e-11f,7.6638e-13f,
	    4.741e-15f,2.4e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y;
    static integer nti1;
    static real xmin, xmax, xsml;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal besi1e_(real *), r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESI1 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            first kind of order one. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI1-S, DBESI1-D) */
/* ***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI1(X) calculates the modified (hyperbolic) Bessel function */
/* of the first kind of order one for real argument X. */

/* Series for BI1        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.40E-17 */
/*                                         log weighted error  16.62 */
/*                               significant figures required  16.23 */
/*                                    decimal places required  17.14 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI1E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESI1 */
/* ***FIRST EXECUTABLE STATEMENT  BESI1 */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nti1 = inits_(bi1cs, &c__11, &r__1);
	xmin = r1mach_(&c__1) * 2.f;
	xsml = sqrt(r1mach_(&c__3) * 4.5f);
	xmax = log(r1mach_(&c__2));
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
	xermsg_("SLATEC", "BESI1", "ABS(X) SO SMALL I1 UNDERFLOWS", &c__1, &
		c__1, (ftnlen)6, (ftnlen)5, (ftnlen)29);
    }
    if (y > xmin) {
	ret_val = *x * .5f;
    }
    if (y > xsml) {
	r__1 = y * y / 4.5f - 1.f;
	ret_val = *x * (csevl_(&r__1, bi1cs, &nti1) + .875f);
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESI1", "ABS(X) SO BIG I1 OVERFLOWS", &c__2, &c__2,
		 (ftnlen)6, (ftnlen)5, (ftnlen)26);
    }

    ret_val = exp(y) * besi1e_(x);

    return ret_val;
} /* besi1_ */

