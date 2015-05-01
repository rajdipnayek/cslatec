/* gamlim.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK GAMLIM */
/* Subroutine */ int gamlim_(real *xmin, real *xmax)
{
    /* System generated locals */
    real r__1, r__2;

    /* Local variables */
    static integer i__;
    static real xln, xold;
    extern doublereal r1mach_(integer *);
    static real alnbig, alnsml;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  GAMLIM */
/* ***PURPOSE  Compute the minimum and maximum bounds for the argument in */
/*            the Gamma function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A, R2 */
/* ***TYPE      SINGLE PRECISION (GAMLIM-S, DGAMLM-D) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Calculate the minimum and maximum legal bounds for X in GAMMA(X). */
/* XMIN and XMAX are not the only bounds, but they are the only non- */
/* trivial ones to calculate. */

/*             Output Arguments -- */
/* XMIN   minimum legal value of X in GAMMA(X).  Any smaller value of */
/*        X might result in underflow. */
/* XMAX   maximum legal value of X in GAMMA(X).  Any larger value will */
/*        cause overflow. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  GAMLIM */
/* ***FIRST EXECUTABLE STATEMENT  GAMLIM */
    alnsml = log(r1mach_(&c__1));
    *xmin = -alnsml;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmin;
	xln = log(*xmin);
	*xmin -= *xmin * ((*xmin + .5f) * xln - *xmin - .2258f + alnsml) / (*
		xmin * xln + .5f);
	if ((r__1 = *xmin - xold, dabs(r__1)) < .005f) {
	    goto L20;
	}
/* L10: */
    }
    xermsg_("SLATEC", "GAMLIM", "UNABLE TO FIND XMIN", &c__1, &c__2, (ftnlen)
	    6, (ftnlen)6, (ftnlen)19);

L20:
    *xmin = -(*xmin) + .01f;

    alnbig = log(r1mach_(&c__2));
    *xmax = alnbig;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmax;
	xln = log(*xmax);
	*xmax -= *xmax * ((*xmax - .5f) * xln - *xmax + .9189f - alnbig) / (*
		xmax * xln - .5f);
	if ((r__1 = *xmax - xold, dabs(r__1)) < .005f) {
	    goto L40;
	}
/* L30: */
    }
    xermsg_("SLATEC", "GAMLIM", "UNABLE TO FIND XMAX", &c__2, &c__2, (ftnlen)
	    6, (ftnlen)6, (ftnlen)19);

L40:
    *xmax += -.01f;
/* Computing MAX */
    r__1 = *xmin, r__2 = -(*xmax) + 1.f;
    *xmin = dmax(r__1,r__2);

    return 0;
} /* gamlim_ */

