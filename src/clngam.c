/* clngam.f -- translated by f2c (version 12.02.01).
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

/* DECK CLNGAM */
/* Complex */ void clngam_(complex * ret_val, complex *zin)
{
    /* Initialized data */

    static real pi = 3.14159265358979324f;
    static real sq2pil = .91893853320467274f;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1;
    real r__1, r__2;
    doublereal d__1, d__2;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8, q__9, q__10, 
	    q__11, q__12, q__13, q__14, q__15, q__16;

    /* Local variables */
    static integer i__, n;
    static real x, y;
    static complex z__;
    extern doublereal carg_(complex *);
    static complex corr;
    static real cabsz, bound, dxrel;
    extern doublereal r1mach_(integer *);
    extern /* Complex */ void c9lgmc_(complex *, complex *), clnrel_(complex *
	    , complex *);
    static real argsum;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CLNGAM */
/* ***PURPOSE  Compute the logarithm of the absolute value of the Gamma */
/*            function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7A */
/* ***TYPE      COMPLEX (ALNGAM-S, DLNGAM-D, CLNGAM-C) */
/* ***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CLNGAM computes the natural log of the complex valued gamma function */
/* at ZIN, where ZIN is a complex number.  This is a preliminary version, */
/* which is not accurate. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  C9LGMC, CARG, CLNREL, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780401  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  CLNGAM */
/* ***FIRST EXECUTABLE STATEMENT  CLNGAM */
    if (first) {
	n = log(r1mach_(&c__3)) * -.3f;
/* BOUND = N*(0.1*EPS)**(-1/(2*N-1))/(PI*EXP(1)) */
	d__1 = (doublereal) (r1mach_(&c__3) * .1f);
	d__2 = (doublereal) (-1.f / ((n << 1) - 1));
	bound = n * .1171f * pow_dd(&d__1, &d__2);
	dxrel = sqrt(r1mach_(&c__4));
    }
    first = FALSE_;

    z__.r = zin->r, z__.i = zin->i;
    x = zin->r;
    y = r_imag(zin);

    corr.r = 0.f, corr.i = 0.f;
    cabsz = c_abs(&z__);
    if (x >= 0.f && cabsz > bound) {
	goto L50;
    }
    if (x < 0.f && dabs(y) > bound) {
	goto L50;
    }

    if (cabsz < bound) {
	goto L20;
    }

/* USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND */
/* ABS(AIMAG(Y)) SMALL. */

    if (y > 0.f) {
	r_cnjg(&q__1, &z__);
	z__.r = q__1.r, z__.i = q__1.i;
    }
    r__1 = pi * 2.f;
    q__4.r = 0.f, q__4.i = r__1;
    q__3.r = -q__4.r, q__3.i = -q__4.i;
    q__2.r = q__3.r * z__.r - q__3.i * z__.i, q__2.i = q__3.r * z__.i + 
	    q__3.i * z__.r;
    c_exp(&q__1, &q__2);
    corr.r = q__1.r, corr.i = q__1.i;
    if (corr.r == 1.f && r_imag(&corr) == 0.f) {
	xermsg_("SLATEC", "CLNGAM", "Z IS A NEGATIVE INTEGER", &c__3, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)23);
    }

    r__1 = sq2pil + 1.f;
    q__7.r = 0.f, q__7.i = pi;
    q__8.r = z__.r - .5f, q__8.i = z__.i;
    q__6.r = q__7.r * q__8.r - q__7.i * q__8.i, q__6.i = q__7.r * q__8.i + 
	    q__7.i * q__8.r;
    q__5.r = r__1 - q__6.r, q__5.i = -q__6.i;
    q__10.r = -corr.r, q__10.i = -corr.i;
    clnrel_(&q__9, &q__10);
    q__4.r = q__5.r - q__9.r, q__4.i = q__5.i - q__9.i;
    q__12.r = z__.r - .5f, q__12.i = z__.i;
    q__14.r = 1.f - z__.r, q__14.i = -z__.i;
    c_log(&q__13, &q__14);
    q__11.r = q__12.r * q__13.r - q__12.i * q__13.i, q__11.i = q__12.r * 
	    q__13.i + q__12.i * q__13.r;
    q__3.r = q__4.r + q__11.r, q__3.i = q__4.i + q__11.i;
    q__2.r = q__3.r - z__.r, q__2.i = q__3.i - z__.i;
    q__16.r = 1.f - z__.r, q__16.i = -z__.i;
    c9lgmc_(&q__15, &q__16);
    q__1.r = q__2.r - q__15.r, q__1.i = q__2.i - q__15.i;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    if (y > 0.f) {
	r_cnjg(&q__1,  ret_val);
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    return ;

/* USE THE RECURSION RELATION FOR ABS(Z) SMALL. */

L20:
    if (x >= -.5f || dabs(y) > dxrel) {
	goto L30;
    }
    r__2 = x - .5f;
    r__1 = r_int(&r__2);
    q__2.r = z__.r - r__1, q__2.i = z__.i;
    q__1.r = q__2.r / x, q__1.i = q__2.i / x;
    if (c_abs(&q__1) < dxrel) {
	xermsg_("SLATEC", "CLNGAM", "ANSWER LT HALF PRECISION BECAUSE Z TOO "
		"NEAR NEGATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)6, (
		ftnlen)60);
    }

L30:
/* Computing 2nd power */
    r__1 = bound;
/* Computing 2nd power */
    r__2 = y;
    n = sqrt(r__1 * r__1 - r__2 * r__2) - x + 1.f;
    argsum = 0.f;
    corr.r = 1.f, corr.i = 0.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	argsum += carg_(&z__);
	q__1.r = z__.r * corr.r - z__.i * corr.i, q__1.i = z__.r * corr.i + 
		z__.i * corr.r;
	corr.r = q__1.r, corr.i = q__1.i;
	q__1.r = z__.r + 1.f, q__1.i = z__.i;
	z__.r = q__1.r, z__.i = q__1.i;
/* L40: */
    }

    if (corr.r == 0.f && r_imag(&corr) == 0.f) {
	xermsg_("SLATEC", "CLNGAM", "Z IS A NEGATIVE INTEGER", &c__3, &c__2, (
		ftnlen)6, (ftnlen)6, (ftnlen)23);
    }
    r__1 = log(c_abs(&corr));
    q__2.r = r__1, q__2.i = argsum;
    q__1.r = -q__2.r, q__1.i = -q__2.i;
    corr.r = q__1.r, corr.i = q__1.i;

/* USE STIRLING-S APPROXIMATION FOR LARGE Z. */

L50:
    q__6.r = z__.r - .5f, q__6.i = z__.i;
    c_log(&q__7, &z__);
    q__5.r = q__6.r * q__7.r - q__6.i * q__7.i, q__5.i = q__6.r * q__7.i + 
	    q__6.i * q__7.r;
    q__4.r = sq2pil + q__5.r, q__4.i = q__5.i;
    q__3.r = q__4.r - z__.r, q__3.i = q__4.i - z__.i;
    q__2.r = q__3.r + corr.r, q__2.i = q__3.i + corr.i;
    c9lgmc_(&q__8, &z__);
    q__1.r = q__2.r + q__8.r, q__1.i = q__2.i + q__8.i;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;
    return ;

} /* clngam_ */

