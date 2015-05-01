/* dcot.f -- translated by f2c (version 12.02.01).
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
static integer c__15 = 15;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b19 = 2.;

/* DECK DCOT */
doublereal dcot_(doublereal *x)
{
    /* Initialized data */

    static doublereal cotcs[15] = { .24025916098295630250955361774497,
	    -.0165330316015002278454746025255758,
	    -4.29983919317240189356476228239895e-5,
	    -1.59283223327541046023490851122445e-7,
	    -6.19109313512934872588620579343187e-10,
	    -2.43019741507264604331702590579575e-12,
	    -9.560936758800080984270620831e-15,
	    -3.76353798194580580416291539706666e-17,
	    -1.48166574646746578852176794666666e-19,
	    -5.83335658903666579477984e-22,
	    -2.29662646964645773928533333333333e-24,
	    -9.04197057307483326719999999999999e-27,-3.55988551920600064e-29,
	    -1.40155139824298666666666666666666e-31,
	    -5.51800436872533333333333333333333e-34 };
    static doublereal pi2rec = .011619772367581343075535053490057;
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static doublereal y;
    static integer ifn;
    static doublereal xmin, yrem, xmax, xsml, ainty, sqeps;
    extern doublereal d1mach_(integer *);
    static doublereal ainty2, prodbg;
    extern doublereal dcsevl_(doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer nterms;

/* ***BEGIN PROLOGUE  DCOT */
/* ***PURPOSE  Compute the cotangent. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      DOUBLE PRECISION (COT-S, DCOT-D, CCOT-C) */
/* ***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DCOT(X) calculates the double precision trigonometric cotangent */
/* for double precision argument X.  X is in units of radians. */

/* Series for COT        on the interval  0.          to  6.25000E-02 */
/*                                        with weighted error   5.52E-34 */
/*                                         log weighted error  33.26 */
/*                               significant figures required  32.34 */
/*                                    decimal places required  33.85 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920618  Removed space from variable names.  (RWC, WRB) */
/* ***END PROLOGUE  DCOT */
/* ***FIRST EXECUTABLE STATEMENT  DCOT */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	nterms = initds_(cotcs, &c__15, &r__1);
	xmax = 1. / d1mach_(&c__4);
	xsml = sqrt(d1mach_(&c__3) * 3.);
/* Computing MAX */
	d__1 = log(d1mach_(&c__1)), d__2 = -log(d1mach_(&c__2));
	xmin = exp(max(d__1,d__2) + .01);
	sqeps = sqrt(d1mach_(&c__4));
    }
    first = FALSE_;

    y = abs(*x);
    if (y < xmin) {
	xermsg_("SLATEC", "DCOT", "ABS(X) IS ZERO OR SO SMALL DCOT OVERFLOWS",
		 &c__2, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)41);
    }
    if (y > xmax) {
	xermsg_("SLATEC", "DCOT", "NO PRECISION BECAUSE ABS(X) IS TOO BIG", &
		c__3, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)38);
    }

/* CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC) */
/* = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z */
/* = AINT(.625*Y) + AINT(Z) + REM(Z) */

    ainty = d_int(&y);
    yrem = y - ainty;
    prodbg = ainty * .625;
    ainty = d_int(&prodbg);
    y = prodbg - ainty + yrem * .625 + pi2rec * y;
    ainty2 = d_int(&y);
    ainty += ainty2;
    y -= ainty2;

    ifn = (integer) d_mod(&ainty, &c_b19);
    if (ifn == 1) {
	y = 1. - y;
    }

    if (abs(*x) > .5 && y < abs(*x) * sqeps) {
	xermsg_("SLATEC", "DCOT", "ANSWER LT HALF PRECISION, ABS(X) TOO BIG "
		"OR X NEAR N*PI (N.NE.0)", &c__1, &c__1, (ftnlen)6, (ftnlen)4, 
		(ftnlen)64);
    }

    if (y > .25) {
	goto L20;
    }
    ret_val = 1. / *x;
    if (y > xsml) {
	d__1 = y * 32. * y - 1.;
	ret_val = (dcsevl_(&d__1, cotcs, &nterms) + .5) / y;
    }
    goto L40;

L20:
    if (y > .5) {
	goto L30;
    }
    d__1 = y * 8. * y - 1.;
    ret_val = (dcsevl_(&d__1, cotcs, &nterms) + .5) / (y * .5);
    ret_val = (ret_val * ret_val - 1.) * .5 / ret_val;
    goto L40;

L30:
    d__1 = y * 2. * y - 1.;
    ret_val = (dcsevl_(&d__1, cotcs, &nterms) + .5) / (y * .25);
    ret_val = (ret_val * ret_val - 1.) * .5 / ret_val;
    ret_val = (ret_val * ret_val - 1.) * .5 / ret_val;

L40:
    if (*x != 0.) {
	ret_val = d_sign(&ret_val, x);
    }
    if (ifn == 1) {
	ret_val = -ret_val;
    }

    return ret_val;
} /* dcot_ */

