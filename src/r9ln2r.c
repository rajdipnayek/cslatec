/* r9ln2r.f -- translated by f2c (version 12.02.01).
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
static integer c__26 = 26;
static integer c__20 = 20;
static integer c__4 = 4;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK R9LN2R */
doublereal r9ln2r_(real *x)
{
    /* Initialized data */

    static real ln21cs[26] = { .1811196251347881f,-.15627123192872463f,
	    .028676305361557275f,-.005558699655948139f,.001117897665229983f,
	    -2.30805089823279e-4f,4.85988533411e-5f,-1.0390127388903e-5f,
	    2.248456370739e-6f,-4.91405927392e-7f,1.0828256507e-7f,
	    -2.4025872763e-8f,5.362460047e-9f,-1.202995136e-9f,
	    2.71078892e-10f,-6.1323562e-11f,1.3920858e-11f,-3.16993e-12f,
	    7.23837e-13f,-1.657e-13f,3.8018e-14f,-8.741e-15f,2.013e-15f,
	    -4.64e-16f,1.07e-16f,-2.4e-17f };
    static real ln22cs[20] = { -.22242532535020461f,-.061047100108078624f,
	    .007427235009750394f,-9.33501826163697e-4f,1.2004990768726e-4f,
	    -1.570472295282e-5f,2.081874781051e-6f,-2.78919557764e-7f,
	    3.7693558237e-8f,-5.130902896e-9f,7.02714117e-10f,-9.6748595e-11f,
	    1.3381046e-11f,-1.858102e-12f,2.58929e-13f,-3.6195e-14f,
	    5.074e-15f,-7.13e-16f,1e-16f,-1.4e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real eps, xbig, xmin, xmax;
    extern doublereal csevl_(real *, real *, integer *);
    static real txbig;
    static integer ntln21, ntln22;
    extern integer inits_(real *, integer *, real *);
    static real sqeps, txmax;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9LN2R */
/* ***SUBSIDIARY */
/* ***PURPOSE  Evaluate LOG(1+X) from second order relative accuracy so */
/*            that LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      SINGLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  LOG(1+X)  from 2-nd order with relative error accuracy so */
/* that    LOG(1+X) = X - X**2/2 + X**3*R9LN2R(X) */

/* Series for LN21       on the interval -6.25000D-01 to  0. */
/*                                        with weighted error   2.49E-17 */
/*                                         log weighted error  16.60 */
/*                               significant figures required  15.87 */
/*                                    decimal places required  17.31 */

/* Series for LN22       on the interval  0.          to  8.12500D-01 */
/*                                        with weighted error   1.42E-17 */
/*                                         log weighted error  16.85 */
/*                               significant figures required  15.95 */
/*                                    decimal places required  17.50 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9LN2R */
/* ***FIRST EXECUTABLE STATEMENT  R9LN2R */
    if (first) {
	eps = r1mach_(&c__3);
	r__1 = eps * .1f;
	ntln21 = inits_(ln21cs, &c__26, &r__1);
	r__1 = eps * .1f;
	ntln22 = inits_(ln22cs, &c__20, &r__1);

	xmin = sqrt(r1mach_(&c__4)) - 1.f;
	sqeps = sqrt(eps);
	txmax = 6.f / sqeps;
/* Computing 2nd power */
	r__1 = txmax;
	xmax = txmax - (eps * (r__1 * r__1) - log(txmax) * 2.f) / (eps * 2.f *
		 txmax);
	txbig = 4.f / sqrt(sqeps);
/* Computing 2nd power */
	r__1 = txbig;
	xbig = txbig - (sqeps * (r__1 * r__1) - log(txbig) * 2.f) / (sqeps * 
		2.f * txbig);
    }
    first = FALSE_;

    if (*x < -.625f || *x > .8125f) {
	goto L20;
    }

    if (*x < 0.f) {
	r__1 = *x * 16.f / 5.f + 1.f;
	ret_val = csevl_(&r__1, ln21cs, &ntln21) + .375f;
    }
    if (*x >= 0.f) {
	r__1 = *x * 32.f / 13.f - 1.f;
	ret_val = csevl_(&r__1, ln22cs, &ntln22) + .375f;
    }
    return ret_val;

L20:
    if (*x < xmin) {
	xermsg_("SLATEC", "R9LN2R", "ANSWER LT HALF PRECISION BECAUSE X IS T"
		"OO NEAR -1", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)49);
    }
    if (*x > xmax) {
	xermsg_("SLATEC", "R9LN2R", "NO PRECISION IN ANSWER BECAUSE X IS TOO"
		" BIG", &c__3, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)43);
    }
    if (*x > xbig) {
	xermsg_("SLATEC", "R9LN2R", "ANSWER LT HALF PRECISION BECAUSE X IS T"
		"OO BIG", &c__2, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)45);
    }

/* Computing 3rd power */
    r__1 = *x;
    ret_val = (log(*x + 1.f) - *x * (1.f - *x * .5f)) / (r__1 * (r__1 * r__1))
	    ;
    return ret_val;

} /* r9ln2r_ */

