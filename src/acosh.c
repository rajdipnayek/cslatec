/* acosh.f -- translated by f2c (version 12.02.01).
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

/* DECK ACOSH */
doublereal acosh_(real *x)
{
    /* Initialized data */

    static real aln2 = .69314718055994530942f;
    static real xmax = 0.f;

    /* System generated locals */
    real ret_val;

    /* Local variables */
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  ACOSH */
/* ***PURPOSE  Compute the arc hyperbolic cosine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      SINGLE PRECISION (ACOSH-S, DACOSH-D, CACOSH-C) */
/* ***KEYWORDS  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB, */
/*             INVERSE HYPERBOLIC COSINE */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ACOSH(X) computes the arc hyperbolic cosine of X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  ACOSH */
/* ***FIRST EXECUTABLE STATEMENT  ACOSH */
    if (xmax == 0.f) {
	xmax = 1.f / sqrt(r1mach_(&c__3));
    }

    if (*x < 1.f) {
	xermsg_("SLATEC", "ACOSH", "X LESS THAN 1", &c__1, &c__2, (ftnlen)6, (
		ftnlen)5, (ftnlen)13);
    }

    if (*x < xmax) {
	ret_val = log(*x + sqrt(*x * *x - 1.f));
    }
    if (*x >= xmax) {
	ret_val = aln2 + log(*x);
    }

    return ret_val;
} /* acosh_ */

