/* ai.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b7 = .3334;
static integer c__1 = 1;
static doublereal c_b9 = .6667;

/* DECK AI */
doublereal ai_(real *x)
{
    /* Initialized data */

    static real aifcs[9] = { -.0379713584966699975f,.05919188853726363857f,
	    9.8629280577279975e-4f,6.84884381907656e-6f,2.594202596219e-8f,
	    6.176612774e-11f,1.0092454e-13f,1.2014e-16f,1e-19f };
    static real aigcs[8] = { .01815236558116127f,.02157256316601076f,
	    2.5678356987483e-4f,1.42652141197e-6f,4.57211492e-9f,9.52517e-12f,
	    1.392e-14f,1e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;
    doublereal d__1;

    /* Local variables */
    static real z__, xm;
    extern doublereal aie_(real *);
    static integer naif, naig;
    static real xmax, x3sml, theta;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    static real xmaxt;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int r9aimp_(real *, real *, real *), xermsg_(char 
	    *, char *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  AI */
/* ***PURPOSE  Evaluate the Airy function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10D */
/* ***TYPE      SINGLE PRECISION (AI-S, DAI-D) */
/* ***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* AI(X) computes the Airy function Ai(X) */
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

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  AIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  AI */
/* ***FIRST EXECUTABLE STATEMENT  AI */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	naif = inits_(aifcs, &c__9, &r__1);
	r__1 = r1mach_(&c__3) * .1f;
	naig = inits_(aigcs, &c__8, &r__1);

	d__1 = (doublereal) r1mach_(&c__3);
	x3sml = pow_dd(&d__1, &c_b7);
	d__1 = (doublereal) (log(r1mach_(&c__1)) * -1.5f);
	xmaxt = pow_dd(&d__1, &c_b9);
	xmax = xmaxt - xmaxt * log(xmaxt) / (sqrt(xmaxt) * 4.f + 1.f) - .01f;
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
    return ret_val;

L30:
    if (*x > xmax) {
	goto L40;
    }
    ret_val = aie_(x) * exp(*x * -2.f * sqrt(*x) / 3.f);
    return ret_val;

L40:
    ret_val = 0.f;
    xermsg_("SLATEC", "AI", "X SO BIG AI UNDERFLOWS", &c__1, &c__1, (ftnlen)6,
	     (ftnlen)2, (ftnlen)22);
    return ret_val;

} /* ai_ */

