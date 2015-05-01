/* albeta.f -- translated by f2c (version 12.02.01).
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

/* DECK ALBETA */
doublereal albeta_(real *a, real *b)
{
    /* Initialized data */

    static real sq2pil = .91893853320467274f;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real p, q, corr;
    extern doublereal gamma_(real *), r9lgmc_(real *), alngam_(real *), 
	    alnrel_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  ALBETA */
/* ***PURPOSE  Compute the natural logarithm of the complete Beta */
/*            function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7B */
/* ***TYPE      SINGLE PRECISION (ALBETA-S, DLBETA-D, CLBETA-C) */
/* ***KEYWORDS  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* ALBETA computes the natural log of the complete beta function. */

/* Input Parameters: */
/*       A   real and positive */
/*       B   real and positive */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALNGAM, ALNREL, GAMMA, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  ALBETA */
/* ***FIRST EXECUTABLE STATEMENT  ALBETA */
    p = dmin(*a,*b);
    q = dmax(*a,*b);

    if (p <= 0.f) {
	xermsg_("SLATEC", "ALBETA", "BOTH ARGUMENTS MUST BE GT ZERO", &c__1, &
		c__2, (ftnlen)6, (ftnlen)6, (ftnlen)30);
    }
    if (p >= 10.f) {
	goto L30;
    }
    if (q >= 10.f) {
	goto L20;
    }

/* P AND Q ARE SMALL. */

    r__1 = p + q;
    ret_val = log(gamma_(&p) * (gamma_(&q) / gamma_(&r__1)));
    return ret_val;

/* P IS SMALL, BUT Q IS BIG. */

L20:
    r__1 = p + q;
    corr = r9lgmc_(&q) - r9lgmc_(&r__1);
    r__1 = -p / (p + q);
    ret_val = alngam_(&p) + corr + p - p * log(p + q) + (q - .5f) * alnrel_(&
	    r__1);
    return ret_val;

/* P AND Q ARE BIG. */

L30:
    r__1 = p + q;
    corr = r9lgmc_(&p) + r9lgmc_(&q) - r9lgmc_(&r__1);
    r__1 = -p / (p + q);
    ret_val = log(q) * -.5f + sq2pil + corr + (p - .5f) * log(p / (p + q)) + 
	    q * alnrel_(&r__1);
    return ret_val;

} /* albeta_ */

