/* gamic.f -- translated by f2c (version 12.02.01).
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
static real c_b17 = 1.f;
static doublereal c_b20 = -.001;

/* DECK GAMIC */
doublereal gamic_(real *a, real *x)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;
    doublereal d__1;

    /* Local variables */
    static real e, h__, t;
    static integer ma;
    static real fm, sga, bot, alx, eps, aeps, sgng, alngs, gstar, sgngs;
    static integer izero;
    static real sqeps, algap1;
    extern doublereal r1mach_(integer *), r9lgic_(real *, real *, real *), 
	    r9gmic_(real *, real *, real *), r9lgit_(real *, real *, real *), 
	    r9gmit_(real *, real *, real *, real *, real *), alngam_(real *);
    extern /* Subroutine */ int algams_(real *, real *, real *);
    static real sgngam, alneps;
    extern /* Subroutine */ int xerclr_(void), xermsg_(char *, char *, char *,
	     integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  GAMIC */
/* ***PURPOSE  Calculate the complementary incomplete Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (GAMIC-S, DGAMIC-D) */
/* ***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*   Evaluate the complementary incomplete gamma function */

/*   GAMIC = integral from X to infinity of EXP(-T) * T**(A-1.)  . */

/*   GAMIC is evaluated for arbitrary real values of A and for non- */
/*   negative values of X (even though GAMIC is defined for X .LT. */
/*   0.0), except that for X = 0 and A .LE. 0.0, GAMIC is undefined. */

/*   GAMIC, A, and X are REAL. */

/*   A slight deterioration of 2 or 3 digits accuracy will occur when */
/*   GAMIC is very large or very small in absolute value, because log- */
/*   arithmic variables are used.  Also, if the parameter A is very close */
/*   to a negative integer (but not a negative integer), there is a loss */
/*   of accuracy, which is reported if the result is less than half */
/*   machine precision. */

/* ***REFERENCES  W. Gautschi, A computational procedure for incomplete */
/*                 gamma functions, ACM Transactions on Mathematical */
/*                 Software 5, 4 (December 1979), pp. 466-481. */
/*               W. Gautschi, Incomplete gamma functions, Algorithm 542, */
/*                 ACM Transactions on Mathematical Software 5, 4 */
/*                 (December 1979), pp. 482-489. */
/* ***ROUTINES CALLED  ALGAMS, ALNGAM, R1MACH, R9GMIC, R9GMIT, R9LGIC, */
/*                    R9LGIT, XERCLR, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920528  DESCRIPTION and REFERENCES sections revised.  (WRB) */
/* ***END PROLOGUE  GAMIC */
/* ***FIRST EXECUTABLE STATEMENT  GAMIC */
    if (first) {
	eps = r1mach_(&c__3) * .5f;
	sqeps = sqrt(r1mach_(&c__4));
	alneps = -log(r1mach_(&c__3));
	bot = log(r1mach_(&c__1));
    }
    first = FALSE_;

    if (*x < 0.f) {
	xermsg_("SLATEC", "GAMIC", "X IS NEGATIVE", &c__2, &c__2, (ftnlen)6, (
		ftnlen)5, (ftnlen)13);
    }

    if (*x > 0.f) {
	goto L20;
    }
    if (*a <= 0.f) {
	xermsg_("SLATEC", "GAMIC", "X = 0 AND A LE 0 SO GAMIC IS UNDEFINED", &
		c__3, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)38);
    }

    r__1 = *a + 1.f;
    ret_val = exp(alngam_(&r__1) - log(*a));
    return ret_val;

L20:
    alx = log(*x);
    sga = 1.f;
    if (*a != 0.f) {
	sga = r_sign(&c_b17, a);
    }
    ma = *a + sga * .5f;
    aeps = *a - ma;

    izero = 0;
    if (*x >= 1.f) {
	goto L60;
    }

    if (*a > .5f || dabs(aeps) > .001f) {
	goto L50;
    }
    fm = (real) (-ma);
    e = 2.f;
    if (fm > 1.f) {
	e = (fm + 2.f) * 2.f / (fm * fm - 1.f);
    }
    d__1 = (doublereal) (*x);
    e -= alx * pow_dd(&d__1, &c_b20);
    if (e * dabs(aeps) > eps) {
	goto L50;
    }

    ret_val = r9gmic_(a, x, &alx);
    return ret_val;

L50:
    r__1 = *a + 1.f;
    algams_(&r__1, &algap1, &sgngam);
    gstar = r9gmit_(a, x, &algap1, &sgngam, &alx);
    if (gstar == 0.f) {
	izero = 1;
    }
    if (gstar != 0.f) {
	alngs = log((dabs(gstar)));
    }
    if (gstar != 0.f) {
	sgngs = r_sign(&c_b17, &gstar);
    }
    goto L70;

L60:
    if (*a < *x) {
	ret_val = exp(r9lgic_(a, x, &alx));
    }
    if (*a < *x) {
	return ret_val;
    }

    sgngam = 1.f;
    r__1 = *a + 1.f;
    algap1 = alngam_(&r__1);
    sgngs = 1.f;
    alngs = r9lgit_(a, x, &algap1);

/* EVALUATION OF GAMIC(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN. */

L70:
    h__ = 1.f;
    if (izero == 1) {
	goto L80;
    }

    t = *a * alx + alngs;
    if (t > alneps) {
	goto L90;
    }
    if (t > -alneps) {
	h__ = 1.f - sgngs * exp(t);
    }

    if (dabs(h__) < sqeps) {
	xerclr_();
    }
    if (dabs(h__) < sqeps) {
	xermsg_("SLATEC", "GAMIC", "RESULT LT HALF PRECISION", &c__1, &c__1, (
		ftnlen)6, (ftnlen)5, (ftnlen)24);
    }

L80:
    sgng = r_sign(&c_b17, &h__) * sga * sgngam;
    t = log((dabs(h__))) + algap1 - log((dabs(*a)));
    if (t < bot) {
	xerclr_();
    }
    ret_val = sgng * exp(t);
    return ret_val;

L90:
    sgng = -sgngs * sga * sgngam;
    t = t + algap1 - log((dabs(*a)));
    if (t < bot) {
	xerclr_();
    }
    ret_val = sgng * exp(t);
    return ret_val;

} /* gamic_ */

