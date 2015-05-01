/* spenc.f -- translated by f2c (version 12.02.01).
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
static integer c__19 = 19;

/* DECK SPENC */
doublereal spenc_(real *x)
{
    /* Initialized data */

    static real spencs[19] = { .1527365598892406f,.08169658058051014f,
	    .00581415714077873f,5.3716198145415e-4f,5.724704675185e-5f,
	    6.67454612164e-6f,8.2764673397e-7f,1.073315673e-7f,
	    1.440077294e-8f,1.98444202e-9f,2.7940058e-10f,4.003991e-11f,
	    5.82346e-12f,8.5767e-13f,1.2768e-13f,1.918e-14f,2.9e-15f,4.4e-16f,
	    6e-17f };
    static real pi26 = 1.644934066848226f;
    static logical first = TRUE_;

    /* System generated locals */
    real ret_val, r__1, r__2;

    /* Local variables */
    static real aln, xbig;
    extern doublereal csevl_(real *, real *, integer *);
    extern integer inits_(real *, integer *, real *);
    extern doublereal r1mach_(integer *);
    static integer nspenc;

/* ***BEGIN PROLOGUE  SPENC */
/* ***PURPOSE  Compute a form of Spence's integral due to K. Mitchell. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C5 */
/* ***TYPE      SINGLE PRECISION (SPENC-S, DSPENC-D) */
/* ***KEYWORDS  FNLIB, SPECIAL FUNCTIONS, SPENCE'S INTEGRAL */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate a form of Spence's function defined by */
/*        integral from 0 to X of  -LOG(1-Y)/Y  DY. */
/* For ABS(X) .LE. 1, the uniformly convergent expansion */
/*        SPENC = sum K=1,infinity  X**K / K**2     is valid. */

/* Spence's function can be used to evaluate much more general integral */
/* forms.  For example, */
/*        integral from 0 to Z of  LOG(A*X+B)/(C*X+D)  DX  = */
/*             LOG(ABS(B-A*D/C))*LOG(ABS(A*(C*X+D)/(A*D-B*C)))/C */
/*             - SPENC (A*(C*Z+D)/(A*D-B*C)) / C. */

/* Ref -- K. Mitchell, Philosophical Magazine, 40, p. 351 (1949). */
/*        Stegun and Abromowitz, AMS 55, p. 1004. */


/* Series for SPEN       on the interval  0.          to  5.00000D-01 */
/*                                        with weighted error   6.82E-17 */
/*                                         log weighted error  16.17 */
/*                               significant figures required  15.22 */
/*                                    decimal places required  16.81 */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CSEVL, INITS, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780201  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  SPENC */
/* ***FIRST EXECUTABLE STATEMENT  SPENC */
    if (first) {
	r__1 = r1mach_(&c__3) * .1f;
	nspenc = inits_(spencs, &c__19, &r__1);
	xbig = 1.f / r1mach_(&c__3);
    }
    first = FALSE_;

    if (*x > 2.f) {
	goto L60;
    }
    if (*x > 1.f) {
	goto L50;
    }
    if (*x > .5f) {
	goto L40;
    }
    if (*x >= 0.f) {
	goto L30;
    }
    if (*x > -1.f) {
	goto L20;
    }

/* HERE IF X .LE. -1.0 */

    aln = log(1.f - *x);
    ret_val = -pi26 - aln * .5f * (log(-(*x)) * 2.f - aln);
    if (*x > -xbig) {
	r__1 = 4.f / (1.f - *x) - 1.f;
	ret_val += (csevl_(&r__1, spencs, &nspenc) + 1.f) / (1.f - *x);
    }
    return ret_val;

/* -1.0 .LT. X .LT. 0.0 */

L20:
/* Computing 2nd power */
    r__1 = log(1.f - *x);
    r__2 = *x * 4.f / (*x - 1.f) - 1.f;
    ret_val = r__1 * r__1 * -.5f - *x * (csevl_(&r__2, spencs, &nspenc) + 1.f)
	     / (*x - 1.f);
    return ret_val;

/* 0.0 .LE. X .LE. 0.5 */

L30:
    r__1 = *x * 4.f - 1.f;
    ret_val = *x * (csevl_(&r__1, spencs, &nspenc) + 1.f);
    return ret_val;

/* 0.5 .LT. X .LE. 1.0 */

L40:
    ret_val = pi26;
    if (*x != 1.f) {
	r__1 = (1.f - *x) * 4.f - 1.f;
	ret_val = pi26 - log(*x) * log(1.f - *x) - (1.f - *x) * (csevl_(&r__1,
		 spencs, &nspenc) + 1.f);
    }
    return ret_val;

/* 1.0 .LT. X .LE. 2.0 */

L50:
/* Computing 2nd power */
    r__1 = *x - 1.f;
    r__2 = (*x - 1.f) * 4.f / *x - 1.f;
    ret_val = pi26 - log(*x) * .5f * log(r__1 * r__1 / *x) + (*x - 1.f) * (
	    csevl_(&r__2, spencs, &nspenc) + 1.f) / *x;
    return ret_val;

/* X .GT. 2.0 */

L60:
/* Computing 2nd power */
    r__1 = log(*x);
    ret_val = pi26 * 2.f - r__1 * r__1 * .5f;
    if (*x < xbig) {
	r__1 = 4.f / *x - 1.f;
	ret_val -= (csevl_(&r__1, spencs, &nspenc) + 1.f) / *x;
    }
    return ret_val;

} /* spenc_ */

