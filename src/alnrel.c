/* alnrel.f -- translated by f2c (version 12.02.01).
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
static integer c__23 = 23;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK ALNREL */
doublereal alnrel_(real *x)
{
    /* Initialized data */

    static real alnrcs[23] = { 1.037869356274377f,-.13364301504908918f,
	    .019408249135520563f,-.003010755112753577f,4.86946147971548e-4f,
	    -8.1054881893175e-5f,1.3778847799559e-5f,-2.380221089435e-6f,
	    4.16404162138e-7f,-7.3595828378e-8f,1.3117611876e-8f,
	    -2.354670931e-9f,4.25227732e-10f,-7.7190894e-11f,1.4075746e-11f,
	    -2.576907e-12f,4.73424e-13f,-8.7249e-14f,1.6124e-14f,-2.987e-15f,
	    5.54e-16f,-1.03e-16f,1.9e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real xmin;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    static integer nlnrel;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  ALNREL */
/* ***PURPOSE  Evaluate ln(1+X) accurate in the sense of relative error. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4B */
/* ***TYPE      SINGLE PRECISION (ALNREL-S, DLNREL-D, CLNREL-C) */
/* ***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ALNREL(X) evaluates ln(1+X) accurately in the sense of relative */
/* error when X is very small.  This routine must be used to */
/* maintain relative error accuracy whenever X is small and */
/* accurately known. */

/* Series for ALNR       on the interval -3.75000D-01 to  3.75000D-01 */
/*                                        with weighted error   1.93E-17 */
/*                                         log weighted error  16.72 */
/*                               significant figures required  16.44 */
/*                                    decimal places required  17.40 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  ALNREL */
/* ***FIRST EXECUTABLE STATEMENT  ALNREL */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nlnrel = inits_(alnrcs, &c__23, &r__1);
	xmin = sqrt(r1mach_(&c__4)) - 1.f;
    }
    first = FALSE_;

    if (*x <= -1.f) {
	xermsg_("SLATEC", "ALNREL", "X IS LE -1", &c__2, &c__2, (ftnlen)6, (
		ftnlen)6, (ftnlen)10);
    }
    if (*x < xmin) {
	xermsg_("SLATEC", "ALNREL", "ANSWER LT HALF PRECISION BECAUSE X TOO "
		"NEAR -1", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)46);
    }

    if (dabs(*x) <= .375f) {
	r__1 = *x / .375f;
	ret_val = *x * (1.f - *x * csevl_(&r__1, alnrcs, &nlnrel));
    }
    if (dabs(*x) > .375f) {
	ret_val = log(*x + 1.f);
    }

    return ret_val;
} /* alnrel_ */

