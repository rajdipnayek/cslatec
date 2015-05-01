/* cbrt.f -- translated by f2c (version 12.02.01).
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

/* DECK CBRT */
doublereal cbrt_(real *x)
{
    /* Initialized data */

    static real cbrt2[5] = { .62996052494743658f,.79370052598409974f,1.f,
	    1.25992104989487316f,1.58740105196819947f };
    static integer niter = 0;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Local variables */
    static integer n;
    static real y;
    static integer irem, iter;
    extern doublereal r9pak_(real *, integer *);
    static integer ixpnt;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int r9upak_(real *, real *, integer *);
    static real cbrtsq;

/* ***BEGIN PROLOGUE  CBRT */
/* ***PURPOSE  Compute the cube root. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C2 */
/* ***TYPE      SINGLE PRECISION (CBRT-S, DCBRT-D, CCBRT-C) */
/* ***KEYWORDS  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CBRT(X) calculates the cube root of X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, R9PAK, R9UPAK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CBRT */
/* ***FIRST EXECUTABLE STATEMENT  CBRT */
    if (niter == 0) {
	niter = 1.443f * log(-.106f * log(.1f * r1mach_(&c__3))) + 1.f;
    }

    ret_val = 0.f;
    if (*x == 0.f) {
	return ret_val;
    }

    r__1 = dabs(*x);
    r9upak_(&r__1, &y, &n);
    ixpnt = n / 3;
    irem = n - ixpnt * 3 + 3;

/* THE APPROXIMATION BELOW IS A GENERALIZED CHEBYSHEV SERIES CONVERTED */
/* TO POLYNOMIAL FORM.  THE APPROX IS NEARLY BEST IN THE SENSE OF */
/* RELATIVE ERROR WITH 4.085 DIGITS ACCURACY. */

    ret_val = y * (y * (y * .144586f - .512653f) + .928549f) + .439581f;

    i__1 = niter;
    for (iter = 1; iter <= i__1; ++iter) {
	cbrtsq = ret_val * ret_val;
	ret_val += (y - ret_val * cbrtsq) / (cbrtsq * 3.f);
/* L10: */
    }

    r__1 = cbrt2[irem - 1] * r_sign(&ret_val, x);
    ret_val = r9pak_(&r__1, &ixpnt);
    return ret_val;

} /* cbrt_ */

