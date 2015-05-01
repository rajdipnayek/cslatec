/* dbesk0.f -- translated by f2c (version 12.02.01).
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
static integer c__16 = 16;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK DBESK0 */
doublereal dbesk0_(doublereal *x)
{
    /* Initialized data */

    static doublereal bk0cs[16] = { -.0353273932339027687201140060063153,
	    .344289899924628486886344927529213,
	    .0359799365153615016265721303687231,
	    .00126461541144692592338479508673447,
	    2.28621210311945178608269830297585e-5,
	    2.53479107902614945730790013428354e-7,
	    1.90451637722020885897214059381366e-9,
	    1.03496952576336245851008317853089e-11,
	    4.25981614279108257652445327170133e-14,
	    1.3744654358807508969423832544e-16,
	    3.57089652850837359099688597333333e-19,
	    7.63164366011643737667498666666666e-22,
	    1.36542498844078185908053333333333e-24,
	    2.07527526690666808319999999999999e-27,2.7128142180729856e-30,
	    3.08259388791466666666666666666666e-33 };
    static logical first = TRUE_;

    /* System generated locals */
    real r__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal y;
    static integer ntk0;
    static doublereal xmax, xsml;
    extern doublereal d1mach_(integer *);
    static doublereal xmaxt;
    extern doublereal dbesi0_(doublereal *), dbsk0e_(doublereal *), dcsevl_(
	    doublereal *, doublereal *, integer *);
    extern integer initds_(doublereal *, integer *, real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBESK0 */
/* ***PURPOSE  Compute the modified (hyperbolic) Bessel function of the */
/*            third kind of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      DOUBLE PRECISION (BESK0-S, DBESK0-D) */
/* ***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS, */
/*             THIRD KIND */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DBESK0(X) calculates the double precision modified (hyperbolic) */
/* Bessel function of the third kind of order zero for double */
/* precision argument X.  The argument must be greater than zero */
/* but not so large that the result underflows. */

/* Series for BK0        on the interval  0.          to  4.00000E+00 */
/*                                        with weighted error   3.08E-33 */
/*                                         log weighted error  32.51 */
/*                               significant figures required  32.05 */
/*                                    decimal places required  33.11 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DBESK0 */
/* ***FIRST EXECUTABLE STATEMENT  DBESK0 */
    if (first) {
	r__1 = (real) d1mach_(&c__3) * .1f;
	ntk0 = initds_(bk0cs, &c__16, &r__1);
	xsml = sqrt(d1mach_(&c__3) * 4.);
	xmaxt = -log(d1mach_(&c__1));
	xmax = xmaxt - xmaxt * .5 * log(xmaxt) / (xmaxt + .5);
    }
    first = FALSE_;

    if (*x <= 0.) {
	xermsg_("SLATEC", "DBESK0", "X IS ZERO OR NEGATIVE", &c__2, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)21);
    }
    if (*x > 2.) {
	goto L20;
    }

    y = 0.;
    if (*x > xsml) {
	y = *x * *x;
    }
    d__1 = y * .5 - 1.;
    ret_val = -log(*x * .5) * dbesi0_(x) - .25 + dcsevl_(&d__1, bk0cs, &ntk0);
    return ret_val;

L20:
    ret_val = 0.;
    if (*x > xmax) {
	xermsg_("SLATEC", "DBESK0", "X SO BIG K0 UNDERFLOWS", &c__1, &c__1, (
		ftnlen)6, (ftnlen)6, (ftnlen)22);
    }
    if (*x > xmax) {
	return ret_val;
    }

    ret_val = exp(-(*x)) * dbsk0e_(x);

    return ret_val;
} /* dbesk0_ */

