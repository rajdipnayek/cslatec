/* r9gmit.f -- translated by f2c (version 12.02.01).
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
static real c_b19 = 1.f;

/* DECK R9GMIT */
doublereal r9gmit_(real *a, real *x, real *algap1, real *sgngam, real *alx)
{
    /* Initialized data */

    static real eps = 0.f;
    static real bot = 0.f;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1;

    /* Local variables */
    static integer k, m;
    static real s, t, ae;
    static integer ma;
    static real fk, te, alg2, algs, aeps, sgng2;
    extern doublereal r1mach_(integer *), alngam_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9GMIT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute Tricomi's incomplete Gamma function for small */
/*            arguments. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (R9GMIT-S, D9GMIT-D) */
/* ***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X, */
/*             SPECIAL FUNCTIONS, TRICOMI */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute Tricomi's incomplete gamma function for small X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALNGAM, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9GMIT */
/* ***FIRST EXECUTABLE STATEMENT  R9GMIT */
    if (eps == 0.f) {
	eps = .5f * r1mach_(&c__3);
    }
    if (bot == 0.f) {
	bot = log(r1mach_(&c__1));
    }

    if (*x <= 0.f) {
	xermsg_("SLATEC", "R9GMIT", "X SHOULD BE GT 0", &c__1, &c__2, (ftnlen)
		6, (ftnlen)6, (ftnlen)16);
    }

    ma = *a + .5f;
    if (*a < 0.f) {
	ma = *a - .5f;
    }
    aeps = *a - ma;

    ae = *a;
    if (*a < -.5f) {
	ae = aeps;
    }

    t = 1.f;
    te = ae;
    s = t;
    for (k = 1; k <= 200; ++k) {
	fk = (real) k;
	te = -(*x) * te / fk;
	t = te / (ae + fk);
	s += t;
	if (dabs(t) < eps * dabs(s)) {
	    goto L30;
	}
/* L20: */
    }
    xermsg_("SLATEC", "R9GMIT", "NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SER"
	    "IES", &c__2, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)46);

L30:
    if (*a >= -.5f) {
	algs = -(*algap1) + log(s);
    }
    if (*a >= -.5f) {
	goto L60;
    }

    r__1 = aeps + 1.f;
    algs = -alngam_(&r__1) + log(s);
    s = 1.f;
    m = -ma - 1;
    if (m == 0) {
	goto L50;
    }
    t = 1.f;
    i__1 = m;
    for (k = 1; k <= i__1; ++k) {
	t = *x * t / (aeps - m - 1 + k);
	s += t;
	if (dabs(t) < eps * dabs(s)) {
	    goto L50;
	}
/* L40: */
    }

L50:
    ret_val = 0.f;
    algs = -ma * log(*x) + algs;
    if (s == 0.f || aeps == 0.f) {
	goto L60;
    }

    sgng2 = *sgngam * r_sign(&c_b19, &s);
    alg2 = -(*x) - *algap1 + log((dabs(s)));

    if (alg2 > bot) {
	ret_val = sgng2 * exp(alg2);
    }
    if (algs > bot) {
	ret_val += exp(algs);
    }
    return ret_val;

L60:
    ret_val = exp(algs);
    return ret_val;

} /* r9gmit_ */

