/* gamit.f -- translated by f2c (version 12.02.01).
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
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__2 = 2;
static real c_b10 = 1.f;

/* DECK GAMIT */
doublereal gamit_(real *a, real *x)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real h__, t, sga, bot, alx, alng;
    extern doublereal gamr_(real *);
    static real aeps, ainta, sqeps, algap1;
    extern doublereal r1mach_(integer *), r9lgic_(real *, real *, real *), 
	    r9lgit_(real *, real *, real *), r9gmit_(real *, real *, real *, 
	    real *, real *), alngam_(real *);
    extern /* Subroutine */ int algams_(real *, real *, real *);
    static real sgngam, alneps;
    extern /* Subroutine */ int xerclr_(void), xermsg_(char *, char *, char *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  GAMIT */
/* ***PURPOSE  Calculate Tricomi's form of the incomplete Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (GAMIT-S, DGAMIT-D) */
/* ***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, */
/*             SPECIAL FUNCTIONS, TRICOMI */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*   Evaluate Tricomi's incomplete gamma function defined by */

/*   GAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) * */
/*             T**(A-1.) */

/*   for A .GT. 0.0 and by analytic continuation for A .LE. 0.0. */
/*   GAMMA(X) is the complete gamma function of X. */

/*   GAMIT is evaluated for arbitrary real values of A and for non- */
/*   negative values of X (even though GAMIT is defined for X .LT. */
/*   0.0), except that for X = 0 and A .LE. 0.0, GAMIT is infinite, */
/*   which is a fatal error. */

/*   The function and both arguments are REAL. */

/*   A slight deterioration of 2 or 3 digits accuracy will occur when */
/*   GAMIT is very large or very small in absolute value, because log- */
/*   arithmic variables are used.  Also, if the parameter  A  is very */
/*   close to a negative integer (but not a negative integer), there is */
/*   a loss of accuracy, which is reported if the result is less than */
/*   half machine precision. */

/* ***REFERENCES  W. Gautschi, A computational procedure for incomplete */
/*                 gamma functions, ACM Transactions on Mathematical */
/*                 Software 5, 4 (December 1979), pp. 466-481. */
/*               W. Gautschi, Incomplete gamma functions, Algorithm 542, */
/*                 ACM Transactions on Mathematical Software 5, 4 */
/*                 (December 1979), pp. 482-489. */
/* ***ROUTINES CALLED  ALGAMS, ALNGAM, GAMR, R1MACH, R9GMIT, R9LGIC, */
/*                    R9LGIT, XERCLR, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920528  DESCRIPTION and REFERENCES sections revised.  (WRB) */
/* ***END PROLOGUE  GAMIT */
/* ***FIRST EXECUTABLE STATEMENT  GAMIT */
    if (first) {
	alneps = -log(r1mach_(&c__3));
	sqeps = sqrt(r1mach_(&c__4));
	bot = log(r1mach_(&c__1));
    }
    first = FALSE_;

    if (*x < 0.f) {
	xermsg_("SLATEC", "GAMIT", "X IS NEGATIVE", &c__2, &c__2, (ftnlen)6, (
		ftnlen)5, (ftnlen)13);
    }

    if (*x != 0.f) {
	alx = log(*x);
    }
    sga = 1.f;
    if (*a != 0.f) {
	sga = r_sign(&c_b10, a);
    }
    r__1 = *a + sga * .5f;
    ainta = r_int(&r__1);
    aeps = *a - ainta;

    if (*x > 0.f) {
	goto L20;
    }
    ret_val = 0.f;
    if (ainta > 0.f || aeps != 0.f) {
	r__1 = *a + 1.f;
	ret_val = gamr_(&r__1);
    }
    return ret_val;

L20:
    if (*x > 1.f) {
	goto L40;
    }
    if (*a >= -.5f || aeps != 0.f) {
	r__1 = *a + 1.f;
	algams_(&r__1, &algap1, &sgngam);
    }
    ret_val = r9gmit_(a, x, &algap1, &sgngam, &alx);
    return ret_val;

L40:
    if (*a < *x) {
	goto L50;
    }
    r__2 = *a + 1.f;
    r__1 = alngam_(&r__2);
    t = r9lgit_(a, x, &r__1);
    if (t < bot) {
	xerclr_();
    }
    ret_val = exp(t);
    return ret_val;

L50:
    alng = r9lgic_(a, x, &alx);

/* EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X)) */

    h__ = 1.f;
    if (aeps == 0.f && ainta <= 0.f) {
	goto L60;
    }
    r__1 = *a + 1.f;
    algams_(&r__1, &algap1, &sgngam);
    t = log((dabs(*a))) + alng - algap1;
    if (t > alneps) {
	goto L70;
    }
    if (t > -alneps) {
	h__ = 1.f - sga * sgngam * exp(t);
    }
    if (dabs(h__) > sqeps) {
	goto L60;
    }
    xerclr_();
    xermsg_("SLATEC", "GAMIT", "RESULT LT HALF PRECISION", &c__1, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)24);

L60:
    t = -(*a) * alx + log((dabs(h__)));
    if (t < bot) {
	xerclr_();
    }
    r__1 = exp(t);
    ret_val = r_sign(&r__1, &h__);
    return ret_val;

L70:
    t -= *a * alx;
    if (t < bot) {
	xerclr_();
    }
    ret_val = -sga * sgngam * exp(t);
    return ret_val;

} /* gamit_ */

