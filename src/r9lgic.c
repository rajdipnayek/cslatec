/* r9lgic.f -- translated by f2c (version 12.02.01).
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
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK R9LGIC */
doublereal r9lgic_(real *a, real *x, real *alx)
{
    /* Initialized data */

    static real eps = 0.f;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer k;
    static real p, r__, s, t, fk, xma, xpa;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9LGIC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the log complementary incomplete Gamma function */
/*            for large X and for A .LE. X. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (R9LGIC-S, D9LGIC-D) */
/* ***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X, */
/*             LOGARITHM, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the log complementary incomplete gamma function for large X */
/* and for A .LE. X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9LGIC */
/* ***FIRST EXECUTABLE STATEMENT  R9LGIC */
    if (eps == 0.f) {
	eps = .5f * r1mach_(&c__3);
    }

    xpa = *x + 1.f - *a;
    xma = *x - 1.f - *a;

    r__ = 0.f;
    p = 1.f;
    s = p;
    for (k = 1; k <= 200; ++k) {
	fk = (real) k;
	t = fk * (*a - fk) * (r__ + 1.f);
	r__ = -t / ((xma + fk * 2.f) * (xpa + fk * 2.f) + t);
	p = r__ * p;
	s += p;
	if (dabs(p) < eps * s) {
	    goto L20;
	}
/* L10: */
    }
    xermsg_("SLATEC", "R9LGIC", "NO CONVERGENCE IN 200 TERMS OF CONTINUED FR"
	    "ACTION", &c__1, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)49);

L20:
    ret_val = *a * *alx - *x + log(s / xpa);

    return ret_val;
} /* r9lgic_ */

