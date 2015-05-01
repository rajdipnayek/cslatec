/* besk0.f -- translated by f2c (version 12.02.01).
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

/* DECK BESK0 */
doublereal besk0_(real *x)
{
    /* Initialized data */

    static real bk0cs[11] = { -.03532739323390276872f,.3442898999246284869f,
	    .03597993651536150163f,.00126461541144692592f,
	    2.286212103119451e-5f,2.5347910790261e-7f,1.90451637722e-9f,
	    1.034969525e-11f,4.259816e-14f,1.3744e-16f,3.5e-19f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y;
    static integer ntk0;
    static real xmax, xsml;
    extern doublereal besi0_(real *), csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    static real xmaxt;
    extern doublereal besk0e_(real *), r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESK0 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            third kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESK0-S, DBESK0-D) */
/* ***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESK0(X) calculates the modified (hyperbolic) Bessel function */
/* of the third kind of order zero for real argument X .GT. 0.0. */

/* Series for BK0        on the interval  0.          to  4.00000D+00 */
/*                                        with weighted error   3.57E-19 */
/*                                         log weighted error  18.45 */
/*                               significant figures required  17.99 */
/*                                    decimal places required  18.97 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI0, BESK0E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESK0 */
/* ***FIRST EXECUTABLE STATEMENT  BESK0 */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	ntk0 = inits_(bk0cs, &c__11, &r__1);
	xsml = sqrt(r1mach_(&c__3) * 4.f);
	xmaxt = -log(r1mach_(&c__1));
	xmax = xmaxt - xmaxt * .5f * log(xmaxt) / (xmaxt + .5f) - .01f;
    }
    first = FALSE_;

    if (*x <= 0.f) {
	xermsg_("SLATEC", "BESK0", "X IS ZERO OR NEGATIVE", &c__2, &c__2, (
		ftnlen)6, (ftnlen)5, (ftnlen)21);
    }
    if (*x > 2.f) {
	goto L20;
    }

    y = 0.f;
    if (*x > xsml) {
	y = *x * *x;
    }
    r__1 = y * .5f - 1.f;
    ret_val = -log(*x * .5f) * besi0_(x) - .25f + csevl_(&r__1, bk0cs, &ntk0);
    return ret_val;

L20:
    ret_val = 0.f;
    if (*x > xmax) {
	xermsg_("SLATEC", "BESK0", "X SO BIG K0 UNDERFLOWS", &c__1, &c__1, (
		ftnlen)6, (ftnlen)5, (ftnlen)22);
    }
    if (*x > xmax) {
	return ret_val;
    }

    ret_val = exp(-(*x)) * besk0e_(x);

    return ret_val;
} /* besk0_ */

