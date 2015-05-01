/* r9atn1.f -- translated by f2c (version 12.02.01).
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
static integer c__21 = 21;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK R9ATN1 */
doublereal r9atn1_(real *x)
{
    /* Initialized data */

    static real atn1cs[21] = { -.03283997535355202f,.05833432343172412f,
	    -.00740036969671964f,.00100978419933728f,-1.4397871635652e-4f,
	    2.114512648992e-5f,-3.17232107425e-6f,4.8366203654e-7f,
	    -7.467746546e-8f,1.164800896e-8f,-1.83208837e-9f,2.9019082e-10f,
	    -4.623885e-11f,7.40552e-12f,-1.19135e-12f,1.924e-13f,-3.118e-14f,
	    5.06e-15f,-8.2e-16f,1.3e-16f,-2e-17f };
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real y, eps, xbig, xmax, xsml;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    static integer ntatn1;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  R9ATN1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Evaluate ATAN(X) from first order relative accuracy so that */
/*            ATAN(X) = X + X**3*R9ATN1(X). */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C4A */
/* ***TYPE      SINGLE PRECISION (R9ATN1-S, D9ATN1-D) */
/* ***KEYWORDS  ARC TANGENT, ELEMENTARY FUNCTIONS, FIRST ORDER, FNLIB, */
/*             TRIGONOMETRIC */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate  ATAN(X)  from first order, that is, evaluate */
/* (ATAN(X)-X)/X**3  with relative error accuracy so that */
/*        ATAN(X) = X + X**3*R9ATN1(X). */

/* Series for ATN1       on the interval  0.          to  1.00000D+00 */
/*                                        with weighted error   2.21E-17 */
/*                                         log weighted error  16.66 */
/*                               significant figures required  15.44 */
/*                                    decimal places required  17.32 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  R9ATN1 */
/* ***FIRST EXECUTABLE STATEMENT  R9ATN1 */
    if (first) {
	eps = r1mach_(&c__3);
	r__1 = eps * .1f;
	ntatn1 = inits_(atn1cs, &c__21, &r__1);

	xsml = sqrt(eps * .1f);
	xbig = 1.571f / sqrt(eps);
	xmax = 1.571f / eps;
    }
    first = FALSE_;

    y = dabs(*x);
    if (y > 1.f) {
	goto L20;
    }

    if (y <= xsml) {
	ret_val = -.33333333333333331f;
    }
    if (y <= xsml) {
	return ret_val;
    }

    r__1 = y * 2.f * y - 1.f;
    ret_val = csevl_(&r__1, atn1cs, &ntatn1) - .25f;
    return ret_val;

L20:
    if (y > xmax) {
	xermsg_("SLATEC", "R9ATN1", "NO PRECISION IN ANSWER BECAUSE X IS TOO"
		" BIG", &c__2, &c__2, (ftnlen)6, (ftnlen)6, (ftnlen)43);
    }
    if (y > xbig) {
	xermsg_("SLATEC", "R9ATN1", "ANSWER LT HALF PRECISION BECAUSE X IS T"
		"OO BIG", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)45);
    }

/* Computing 3rd power */
    r__1 = *x;
    ret_val = (atan(*x) - *x) / (r__1 * (r__1 * r__1));
    return ret_val;

} /* r9atn1_ */

