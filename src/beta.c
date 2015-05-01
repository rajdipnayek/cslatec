/* beta.f -- translated by f2c (version 12.02.01).
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

/* DECK BETA */
doublereal beta_(real *a, real *b)
{
    /* Initialized data */

    static real xmax = 0.f;
    static real alnsml = 0.f;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real xmin;
    extern doublereal gamma_(real *), r1mach_(integer *), albeta_(real *, 
	    real *);
    extern /* Subroutine */ int gamlim_(real *, real *), xermsg_(char *, char 
	    *, char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BETA */
/* ***PURPOSE  Compute the complete Beta function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7B */
/* ***TYPE      SINGLE PRECISION (BETA-S, DBETA-D, CBETA-C) */
/* ***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* BETA computes the complete beta function. */

/* Input Parameters: */
/*       A   real and positive */
/*       B   real and positive */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALBETA, GAMLIM, GAMMA, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  BETA */
/* ***FIRST EXECUTABLE STATEMENT  BETA */
    if (alnsml != 0.f) {
	goto L10;
    }
    gamlim_(&xmin, &xmax);
    alnsml = log(r1mach_(&c__1));

L10:
    if (*a <= 0.f || *b <= 0.f) {
	xermsg_("SLATEC", "BETA", "BOTH ARGUMENTS MUST BE GT 0", &c__2, &c__2,
		 (ftnlen)6, (ftnlen)4, (ftnlen)27);
    }

    if (*a + *b < xmax) {
	r__1 = *a + *b;
	ret_val = gamma_(a) * gamma_(b) / gamma_(&r__1);
    }
    if (*a + *b < xmax) {
	return ret_val;
    }

    ret_val = albeta_(a, b);
    if (ret_val < alnsml) {
	xermsg_("SLATEC", "BETA", "A AND/OR B SO BIG BETA UNDERFLOWS", &c__1, 
		&c__2, (ftnlen)6, (ftnlen)4, (ftnlen)33);
    }

    ret_val = exp(ret_val);

    return ret_val;
} /* beta_ */

