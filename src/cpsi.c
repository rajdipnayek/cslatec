/* cpsi.f -- translated by f2c (version 12.02.01).
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
static complex c_b28 = {1.f,0.f};
static complex c_b34 = {.5f,0.f};

/* DECK CPSI */
/* Complex */ void cpsi_(complex * ret_val, complex *zin)
{
    /* Initialized data */

    static real bern[13] = { .083333333333333333f,-.0083333333333333333f,
	    .0039682539682539683f,-.0041666666666666667f,
	    .0075757575757575758f,-.021092796092796093f,.083333333333333333f,
	    -.44325980392156863f,3.0539543302701197f,-26.456212121212121f,
	    281.46014492753623f,-3454.8853937728938f,54827.583333333333f };
    static real pi = 3.141592653589793f;
    static logical first = TRUE_;

    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;
    doublereal d__1, d__2;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Local variables */
    static integer i__, n;
    static real x, y;
    static complex z__;
    static integer ndx;
    static real rbig;
    extern /* Complex */ void ccot_(complex *, complex *);
    static complex corr;
    static real rmin;
    static complex z2inv;
    static real cabsz, bound, dxrel;
    static integer nterm;
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CPSI */
/* ***PURPOSE  Compute the Psi (or Digamma) function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C7C */
/* ***TYPE      COMPLEX (PSI-S, DPSI-D, CPSI-C) */
/* ***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* PSI(X) calculates the psi (or digamma) function of X.  PSI(X) */
/* is the logarithmic derivative of the gamma function of X. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CCOT, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780501  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  CPSI */
/* ***FIRST EXECUTABLE STATEMENT  CPSI */
    if (first) {
	nterm = log(r1mach_(&c__3)) * -.3f;
/* MAYBE BOUND = N*(0.1*EPS)**(-1/(2*N-1)) / (PI*EXP(1)) */
	d__1 = (doublereal) (r1mach_(&c__3) * .1f);
	d__2 = (doublereal) (-1.f / ((nterm << 1) - 1));
	bound = nterm * .1171f * pow_dd(&d__1, &d__2);
	dxrel = sqrt(r1mach_(&c__4));
/* Computing MAX */
	r__1 = log(r1mach_(&c__1)), r__2 = -log(r1mach_(&c__2));
	rmin = exp(dmax(r__1,r__2) + .011f);
	rbig = 1.f / r1mach_(&c__3);
    }
    first = FALSE_;

    z__.r = zin->r, z__.i = zin->i;
    x = z__.r;
    y = r_imag(&z__);
    if (y < 0.f) {
	r_cnjg(&q__1, &z__);
	z__.r = q__1.r, z__.i = q__1.i;
    }

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

    r__1 = -pi;
    q__3.r = pi * z__.r, q__3.i = pi * z__.i;
    ccot_(&q__2, &q__3);
    q__1.r = r__1 * q__2.r, q__1.i = r__1 * q__2.i;
    corr.r = q__1.r, corr.i = q__1.i;
    q__1.r = 1.f - z__.r, q__1.i = -z__.i;
    z__.r = q__1.r, z__.i = q__1.i;
    goto L50;

/* USE THE RECURSION RELATION FOR ABS(Z) SMALL. */

L20:
    if (cabsz < rmin) {
	xermsg_("SLATEC", "CPSI", "CPSI CALLED WITH Z SO NEAR 0 THAT CPSI OV"
		"ERFLOWS", &c__2, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)48);
    }

    if (x >= -.5f || dabs(y) > dxrel) {
	goto L30;
    }
    r__2 = x - .5f;
    r__1 = r_int(&r__2);
    q__2.r = z__.r - r__1, q__2.i = z__.i;
    q__1.r = q__2.r / x, q__1.i = q__2.i / x;
    if (c_abs(&q__1) < dxrel) {
	xermsg_("SLATEC", "CPSI", "ANSWER LT HALF PRECISION BECAUSE Z TOO NE"
		"AR NEGATIVE INTEGER", &c__1, &c__1, (ftnlen)6, (ftnlen)4, (
		ftnlen)60);
    }
    if (y == 0.f && x == r_int(&x)) {
	xermsg_("SLATEC", "CPSI", "Z IS A NEGATIVE INTEGER", &c__3, &c__2, (
		ftnlen)6, (ftnlen)4, (ftnlen)23);
    }

L30:
/* Computing 2nd power */
    r__1 = bound;
/* Computing 2nd power */
    r__2 = y;
    n = sqrt(r__1 * r__1 - r__2 * r__2) - x + 1.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	c_div(&q__2, &c_b28, &z__);
	q__1.r = corr.r - q__2.r, q__1.i = corr.i - q__2.i;
	corr.r = q__1.r, corr.i = q__1.i;
	q__1.r = z__.r + 1.f, q__1.i = z__.i;
	z__.r = q__1.r, z__.i = q__1.i;
/* L40: */
    }

/* NOW EVALUATE THE ASYMPTOTIC SERIES FOR SUITABLY LARGE Z. */

L50:
    if (cabsz > rbig) {
	c_log(&q__2, &z__);
	q__1.r = q__2.r + corr.r, q__1.i = q__2.i + corr.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }
    if (cabsz > rbig) {
	goto L70;
    }

     ret_val->r = 0.f,  ret_val->i = 0.f;
    pow_ci(&q__2, &z__, &c__2);
    c_div(&q__1, &c_b28, &q__2);
    z2inv.r = q__1.r, z2inv.i = q__1.i;
    i__1 = nterm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ndx = nterm + 1 - i__;
	i__2 = ndx - 1;
	q__2.r = z2inv.r *  ret_val->r - z2inv.i *  ret_val->i, q__2.i = 
		z2inv.r *  ret_val->i + z2inv.i *  ret_val->r;
	q__1.r = bern[i__2] + q__2.r, q__1.i = q__2.i;
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
/* L60: */
    }
    c_log(&q__4, &z__);
    c_div(&q__5, &c_b34, &z__);
    q__3.r = q__4.r - q__5.r, q__3.i = q__4.i - q__5.i;
    q__6.r =  ret_val->r * z2inv.r -  ret_val->i * z2inv.i, q__6.i =  
	    ret_val->r * z2inv.i +  ret_val->i * z2inv.r;
    q__2.r = q__3.r - q__6.r, q__2.i = q__3.i - q__6.i;
    q__1.r = q__2.r + corr.r, q__1.i = q__2.i + corr.i;
     ret_val->r = q__1.r,  ret_val->i = q__1.i;

L70:
    if (y < 0.f) {
	r_cnjg(&q__1,  ret_val);
	 ret_val->r = q__1.r,  ret_val->i = q__1.i;
    }

    return ;
} /* cpsi_ */

