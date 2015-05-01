/* r9lgit.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK R9LGIT */
doublereal r9lgit_(real *a, real *x, real *algap1)
{
    /* Initialized data */

    static real eps = 0.f;
    static real sqeps = 0.f;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    static integer k;
    static real p, r__, s, t, fk, ax, a1x, hstar;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9LGIT */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the logarithm of Tricomi's incomplete Gamma */
/*            function with Perron's continued fraction for large X and */
/*            A .GE. X. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (R9LGIT-S, D9LGIT-D) */
/* ***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM, */
/*             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the log of Tricomi's incomplete gamma function with Perron's */
/* continued fraction for large X and for A .GE. X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9LGIT */
/* ***FIRST EXECUTABLE STATEMENT  R9LGIT */
    if (eps == 0.f) {
	eps = .5f * r1mach_(&c__3);
    }
    if (sqeps == 0.f) {
	sqeps = sqrt(r1mach_(&c__4));
    }

    if (*x <= 0.f || *a < *x) {
	xermsg_("SLATEC", "R9LGIT", "X SHOULD BE GT 0.0 AND LE A", &c__2, &
		c__2, (ftnlen)6, (ftnlen)6, (ftnlen)27);
    }

    ax = *a + *x;
    a1x = ax + 1.f;
    r__ = 0.f;
    p = 1.f;
    s = p;
    for (k = 1; k <= 200; ++k) {
	fk = (real) k;
	t = (*a + fk) * *x * (r__ + 1.f);
	r__ = t / ((ax + fk) * (a1x + fk) - t);
	p = r__ * p;
	s += p;
	if (dabs(p) < eps * s) {
	    goto L30;
	}
/* L20: */
    }
    xermsg_("SLATEC", "R9LGIT", "NO CONVERGENCE IN 200 TERMS OF CONTINUED FR"
	    "ACTION", &c__3, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)49);

L30:
    hstar = 1.f - *x * s / a1x;
    if (hstar < sqeps) {
	xermsg_("SLATEC", "R9LGIT", "RESULT LESS THAN HALF PRECISION", &c__1, 
		&c__1, (ftnlen)6, (ftnlen)6, (ftnlen)31);
    }

    ret_val = -(*x) - *algap1 - log(hstar);

    return ret_val;
} /* r9lgit_ */

