/* dgamic.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b17 = 1.;
static doublereal c_b20 = -.001;

/* DECK DGAMIC */
doublereal dgamic_(doublereal *a, doublereal *x)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal e, h__, t, sga, alx, bot, eps, aeps, sgng, ainta, alngs,
	     gstar, sgngs;
    static integer izero;
    static doublereal sqeps;
    extern doublereal d1mach_(integer *);
    static doublereal algap1;
    extern doublereal d9lgic_(doublereal *, doublereal *, doublereal *), 
	    d9gmic_(doublereal *, doublereal *, doublereal *), d9lgit_(
	    doublereal *, doublereal *, doublereal *), d9gmit_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlngam_(
	    doublereal *);
    extern /* Subroutine */ int dlgams_(doublereal *, doublereal *, 
	    doublereal *);
    static doublereal sgngam, alneps;
    extern /* Subroutine */ int xerclr_(void), xermsg_(char *, char *, char *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DGAMIC */
/* ***PURPOSE  Calculate the complementary incomplete Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      DOUBLE PRECISION (GAMIC-S, DGAMIC-D) */
/* ***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*   Evaluate the complementary incomplete Gamma function */

/*   DGAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  . */

/*   DGAMIC is evaluated for arbitrary real values of A and for non- */
/*   negative values of X (even though DGAMIC is defined for X .LT. */
/*   0.0), except that for X = 0 and A .LE. 0.0, DGAMIC is undefined. */

/*   DGAMIC, A, and X are DOUBLE PRECISION. */

/*   A slight deterioration of 2 or 3 digits accuracy will occur when */
/*   DGAMIC is very large or very small in absolute value, because log- */
/*   arithmic variables are used.  Also, if the parameter A is very close */
/*   to a negative INTEGER (but not a negative integer), there is a loss */
/*   of accuracy, which is reported if the result is less than half */
/*   machine precision. */

/* ***REFERENCES  W. Gautschi, A computational procedure for incomplete */
/*                 gamma functions, ACM Transactions on Mathematical */
/*                 Software 5, 4 (December 1979), pp. 466-481. */
/*               W. Gautschi, Incomplete gamma functions, Algorithm 542, */
/*                 ACM Transactions on Mathematical Software 5, 4 */
/*                 (December 1979), pp. 482-489. */
/* ***ROUTINES CALLED  D1MACH, D9GMIC, D9GMIT, D9LGIC, D9LGIT, DLGAMS, */
/*                    DLNGAM, XERCLR, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920528  DESCRIPTION and REFERENCES sections revised.  (WRB) */
/* ***END PROLOGUE  DGAMIC */
/* ***FIRST EXECUTABLE STATEMENT  DGAMIC */
    if (first) {
	eps = d1mach_(&c__3) * .5;
	sqeps = sqrt(d1mach_(&c__4));
	alneps = -log(d1mach_(&c__3));
	bot = log(d1mach_(&c__1));
    }
    first = FALSE_;

    if (*x < 0.) {
	xermsg_("SLATEC", "DGAMIC", "X IS NEGATIVE", &c__2, &c__2, (ftnlen)6, 
		(ftnlen)6, (ftnlen)13);
    }

    if (*x > 0.) {
	goto L20;
    }
    if (*a <= 0.) {
	xermsg_("SLATEC", "DGAMIC", "X = 0 AND A LE 0 SO DGAMIC IS UNDEFINED",
		 &c__3, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)39);
    }

    d__1 = *a + 1.;
    ret_val = exp(dlngam_(&d__1) - log(*a));
    return ret_val;

L20:
    alx = log(*x);
    sga = 1.;
    if (*a != 0.) {
	sga = d_sign(&c_b17, a);
    }
    d__1 = *a + sga * .5;
    ainta = d_int(&d__1);
    aeps = *a - ainta;

    izero = 0;
    if (*x >= 1.) {
	goto L40;
    }

    if (*a > .5 || abs(aeps) > .001) {
	goto L30;
    }
    e = 2.;
    if (-ainta > 1.) {
	e = (-ainta + 2.) * 2. / (ainta * ainta - 1.);
    }
    e -= alx * pow_dd(x, &c_b20);
    if (e * abs(aeps) > eps) {
	goto L30;
    }

    ret_val = d9gmic_(a, x, &alx);
    return ret_val;

L30:
    d__1 = *a + 1.;
    dlgams_(&d__1, &algap1, &sgngam);
    gstar = d9gmit_(a, x, &algap1, &sgngam, &alx);
    if (gstar == 0.) {
	izero = 1;
    }
    if (gstar != 0.) {
	alngs = log((abs(gstar)));
    }
    if (gstar != 0.) {
	sgngs = d_sign(&c_b17, &gstar);
    }
    goto L50;

L40:
    if (*a < *x) {
	ret_val = exp(d9lgic_(a, x, &alx));
    }
    if (*a < *x) {
	return ret_val;
    }

    sgngam = 1.;
    d__1 = *a + 1.;
    algap1 = dlngam_(&d__1);
    sgngs = 1.;
    alngs = d9lgit_(a, x, &algap1);

/* EVALUATION OF DGAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN. */

L50:
    h__ = 1.;
    if (izero == 1) {
	goto L60;
    }

    t = *a * alx + alngs;
    if (t > alneps) {
	goto L70;
    }
    if (t > -alneps) {
	h__ = 1. - sgngs * exp(t);
    }

    if (abs(h__) < sqeps) {
	xerclr_();
    }
    if (abs(h__) < sqeps) {
	xermsg_("SLATEC", "DGAMIC", "RESULT LT HALF PRECISION", &c__1, &c__1, 
		(ftnlen)6, (ftnlen)6, (ftnlen)24);
    }

L60:
    sgng = d_sign(&c_b17, &h__) * sga * sgngam;
    t = log((abs(h__))) + algap1 - log((abs(*a)));
    if (t < bot) {
	xerclr_();
    }
    ret_val = sgng * exp(t);
    return ret_val;

L70:
    sgng = -sgngs * sga * sgngam;
    t = t + algap1 - log((abs(*a)));
    if (t < bot) {
	xerclr_();
    }
    ret_val = sgng * exp(t);
    return ret_val;

} /* dgamic_ */

