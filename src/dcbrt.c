/* dcbrt.f -- translated by f2c (version 12.02.01).
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

/* DECK DCBRT */
doublereal dcbrt_(doublereal *x)
{
    /* Initialized data */

    static doublereal cbrt2[5] = { .62996052494743658238360530363911,
	    .79370052598409973737585281963615,1.,
	    1.25992104989487316476721060727823,
	    1.58740105196819947475170563927231 };
    static integer niter = 0;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer n;
    static doublereal y;
    static real z__;
    static integer irem, iter;
    extern doublereal d9pak_(doublereal *, integer *), d1mach_(integer *);
    static integer ixpnt;
    extern /* Subroutine */ int d9upak_(doublereal *, doublereal *, integer *)
	    ;
    static doublereal cbrtsq;

/* ***BEGIN PROLOGUE  DCBRT */
/* ***PURPOSE  Compute the cube root. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C2 */
/* ***TYPE      DOUBLE PRECISION (CBRT-S, DCBRT-D, CCBRT-C) */
/* ***KEYWORDS  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DCBRT(X) calculates the double precision cube root for */
/* double precision argument X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9PAK, D9UPAK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  DCBRT */
/* ***FIRST EXECUTABLE STATEMENT  DCBRT */
    if (niter == 0) {
	niter = 1.443f * log(-.106f * log(.1f * (real) d1mach_(&c__3))) + 1.f;
    }

    ret_val = 0.;
    if (*x == 0.) {
	return ret_val;
    }

    d__1 = abs(*x);
    d9upak_(&d__1, &y, &n);
    ixpnt = n / 3;
    irem = n - ixpnt * 3 + 3;

/* THE APPROXIMATION BELOW IS A GENERALIZED CHEBYSHEV SERIES CONVERTED */
/* TO POLYNOMIAL FORM.  THE APPROX IS NEARLY BEST IN THE SENSE OF */
/* RELATIVE ERROR WITH 4.085 DIGITS ACCURACY. */

    z__ = y;
    ret_val = z__ * (z__ * (z__ * .144586f - .512653f) + .928549f) + .439581f;

    i__1 = niter;
    for (iter = 1; iter <= i__1; ++iter) {
	cbrtsq = ret_val * ret_val;
	ret_val += (y - ret_val * cbrtsq) / (cbrtsq * 3.);
/* L10: */
    }

    d__1 = cbrt2[irem - 1] * d_sign(&ret_val, x);
    ret_val = d9pak_(&d__1, &ixpnt);
    return ret_val;

} /* dcbrt_ */

