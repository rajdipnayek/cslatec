/* r9lgmc.f -- translated by f2c (version 12.02.01).
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

static integer c__6 = 6;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__2 = 2;

/* DECK R9LGMC */
doublereal r9lgmc_(real *x)
{
    /* Initialized data */

    static real algmcs[6] = { .166638948045186f,-1.38494817606e-5f,
	    9.8108256e-9f,-1.80912e-11f,6.22e-14f,-3e-16f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real xbig, xmax;
    static integer nalgm;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9LGMC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Compute the log Gamma correction factor so that */
/*            LOG(GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X */
/*            + R9LGMC(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7E */
/* ***TYPE      SINGLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C) */
/* ***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB, */
/*             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Compute the log gamma correction factor for X .GE. 10.0 so that */
/*  LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X) */

/* Series for ALGM       on the interval  0.          to  1.00000D-02 */
/*                                        with weighted error   3.40E-16 */
/*                                         log weighted error  15.47 */
/*                               significant figures required  14.39 */
/*                                    decimal places required  15.86 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9LGMC */
/* ***FIRST EXECUTABLE STATEMENT  R9LGMC */
    if (first) {
	r__1 = r1mach_(&c__3);
	nalgm = inits_(algmcs, &c__6, &r__1);
	xbig = 1.f / sqrt(r1mach_(&c__3));
/* Computing MIN */
	r__1 = log(r1mach_(&c__2) / 12.f), r__2 = -log(r1mach_(&c__1) * 12.f);
	xmax = exp((dmin(r__1,r__2)));
    }
    first = FALSE_;

    if (*x < 10.f) {
	xermsg_("SLATEC", "R9LGMC", "X MUST BE GE 10", &c__1, &c__2, (ftnlen)
		6, (ftnlen)6, (ftnlen)15);
    }
    if (*x >= xmax) {
	goto L20;
    }

    ret_val = 1.f / (*x * 12.f);
    if (*x < xbig) {
/* Computing 2nd power */
	r__2 = 10.f / *x;
	r__1 = r__2 * r__2 * 2.f - 1.f;
	ret_val = csevl_(&r__1, algmcs, &nalgm) / *x;
    }
    return ret_val;

L20:
    ret_val = 0.f;
    xermsg_("SLATEC", "R9LGMC", "X SO BIG R9LGMC UNDERFLOWS", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)26);
    return ret_val;

} /* r9lgmc_ */

