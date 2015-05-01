/* dacosh.f -- translated by f2c (version 12.02.01).
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

/* DECK DACOSH */
doublereal dacosh_(doublereal *x)
{
    /* Initialized data */

    static doublereal dln2 = .69314718055994530941723212145818;
    static doublereal xmax = 0.;

    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DACOSH */
/* ***PURPOSE  Compute the arc hyperbolic cosine. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4C */
/* ***TYPE      DOUBLE PRECISION (ACOSH-S, DACOSH-D, CACOSH-C) */
/* ***KEYWORDS  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB, */
/*             INVERSE HYPERBOLIC COSINE */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DACOSH(X) calculates the double precision arc hyperbolic cosine for */
/* double precision argument X.  The result is returned on the */
/* positive branch. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  DACOSH */
/* ***FIRST EXECUTABLE STATEMENT  DACOSH */
    if (xmax == 0.) {
	xmax = 1. / sqrt(d1mach_(&c__3));
    }

    if (*x < 1.) {
	xermsg_("SLATEC", "DACOSH", "X LESS THAN 1", &c__1, &c__2, (ftnlen)6, 
		(ftnlen)6, (ftnlen)13);
    }

    if (*x < xmax) {
	ret_val = log(*x + sqrt(*x * *x - 1.));
    }
    if (*x >= xmax) {
	ret_val = dln2 + log(*x);
    }

    return ret_val;
} /* dacosh_ */

