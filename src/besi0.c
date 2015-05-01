/* besi0.f -- translated by f2c (version 12.02.01).
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
static integer c__12 = 12;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BESI0 */
doublereal besi0_(real *x)
{
    /* Initialized data */

    static real bi0cs[12] = { -.07660547252839144951f,1.92733795399380827f,
	    .2282644586920301339f,.01304891466707290428f,
	    4.3442709008164874e-4f,9.42265768600193e-6f,1.4340062895106e-7f,
	    1.61384906966e-9f,1.396650044e-11f,9.579451e-14f,5.3339e-16f,
	    2.45e-18f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y;
    static integer nti0;
    static real xmax, xsml;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal besi0e_(real *), r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BESI0 */
/* ***PURPOSE  Compute the hyperbolic Bessel function of the first kind */
/*            of order zero. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C10B1 */
/* ***TYPE      SINGLE PRECISION (BESI0-S, DBESI0-D) */
/* ***KEYWORDS  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION, */
/*             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BESI0(X) computes the modified (hyperbolic) Bessel function */
/* of the first kind of order zero and real argument X. */

/* Series for BI0        on the interval  0.          to  9.00000D+00 */
/*                                        with weighted error   2.46E-18 */
/*                                         log weighted error  17.61 */
/*                               significant figures required  17.90 */
/*                                    decimal places required  18.15 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BESI0E, CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  BESI0 */
/* ***FIRST EXECUTABLE STATEMENT  BESI0 */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nti0 = inits_(bi0cs, &c__12, &r__1);
	xsml = sqrt(r1mach_(&c__3) * 4.5f);
	xmax = log(r1mach_(&c__2));
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 3.f) {
	goto L20;
    }

    ret_val = 1.f;
    if (y > xsml) {
	r__1 = y * y / 4.5f - 1.f;
	ret_val = csevl_(&r__1, bi0cs, &nti0) + 2.75f;
    }
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "BESI0", "ABS(X) SO BIG I0 OVERFLOWS", &c__1, &c__2,
		 (ftnlen)6, (ftnlen)5, (ftnlen)26);
    }

    ret_val = exp(y) * besi0e_(x);

    return ret_val;
} /* besi0_ */

