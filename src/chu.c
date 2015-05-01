/* chu.f -- translated by f2c (version 12.02.01).
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
static integer c__10 = 10;
static real c_b25 = -1.f;

/* DECK CHU */
doublereal chu_(real *a, real *b, real *x)
{
    /* Initialized data */

    static real pi = 3.14159265358979324f;
    static real eps = 0.f;

    /* System generated locals */
    integer i__1;
    real ret_val, r__1, r__2, r__3;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, m, n;
    static real t, a0, b0, c0, xi, xn, xi1, sum;
    extern doublereal gamr_(real *);
    static real beps;
    extern doublereal poch_(real *, real *);
    static real alnx, pch1i;
    extern doublereal poch1_(real *, real *), r9chu_(real *, real *, real *);
    static real xeps1;
    extern doublereal gamma_(real *);
    static real aintb;
    static integer istrt;
    static real pch1ai;
    extern doublereal r1mach_(integer *);
    static real gamri1, pochai, gamrni, factor;
    extern doublereal exprel_(real *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static real xtoeps;

/* ***BEGIN PROLOGUE  CHU */
/* ***PURPOSE  Compute the logarithmic confluent hypergeometric function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C11 */
/* ***TYPE      SINGLE PRECISION (CHU-S, DCHU-D) */
/* ***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* CHU computes the logarithmic confluent hypergeometric function, */
/* U(A,B,X). */

/* Input Parameters: */
/*       A   real */
/*       B   real */
/*       X   real and positive */

/* This routine is not valid when 1+A-B is close to zero if X is small. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  EXPREL, GAMMA, GAMR, POCH, POCH1, R1MACH, R9CHU, */
/*                    XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  CHU */
/* ***FIRST EXECUTABLE STATEMENT  CHU */
    if (eps == 0.f) {
	eps = r1mach_(&c__3);
    }

    if (*x == 0.f) {
	xermsg_("SLATEC", "CHU", "X IS ZERO SO CHU IS INFINITE", &c__1, &c__2,
		 (ftnlen)6, (ftnlen)3, (ftnlen)28);
    }
    if (*x < 0.f) {
	xermsg_("SLATEC", "CHU", "X IS NEGATIVE, USE CCHU", &c__2, &c__2, (
		ftnlen)6, (ftnlen)3, (ftnlen)23);
    }

/* Computing MAX */
    r__2 = dabs(*a);
/* Computing MAX */
    r__3 = (r__1 = *a + 1.f - *b, dabs(r__1));
    if (dmax(r__2,1.f) * dmax(r__3,1.f) < dabs(*x) * .99f) {
	goto L120;
    }

/* THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL */
/* APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE. */

    if ((r__1 = *a + 1.f - *b, dabs(r__1)) < sqrt(eps)) {
	xermsg_("SLATEC", "CHU", "ALGORITHM IS BAD WHEN 1+A-B IS NEAR ZERO F"
		"OR SMALL X", &c__10, &c__2, (ftnlen)6, (ftnlen)3, (ftnlen)52);
    }

    r__1 = *b + .5f;
    aintb = r_int(&r__1);
    if (*b < 0.f) {
	r__1 = *b - .5f;
	aintb = r_int(&r__1);
    }
    beps = *b - aintb;
    n = aintb;

    alnx = log(*x);
    xtoeps = exp(-beps * alnx);

/* EVALUATE THE FINITE SUM.     ----------------------------------------- */

    if (n >= 1) {
	goto L40;
    }

/* CONSIDER THE CASE B .LT. 1.0 FIRST. */

    sum = 1.f;
    if (n == 0) {
	goto L30;
    }

    t = 1.f;
    m = -n;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi1 = (real) (i__ - 1);
	t = t * (*a + xi1) * *x / ((*b + xi1) * (xi1 + 1.f));
	sum += t;
/* L20: */
    }

L30:
    r__1 = *a + 1.f - *b;
    r__2 = -(*a);
    sum = poch_(&r__1, &r__2) * sum;
    goto L70;

/* NOW CONSIDER THE CASE B .GE. 1.0. */

L40:
    sum = 0.f;
    m = n - 2;
    if (m < 0) {
	goto L70;
    }
    t = 1.f;
    sum = 1.f;
    if (m == 0) {
	goto L60;
    }

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = (real) i__;
	t = t * (*a - *b + xi) * *x / ((1.f - *b + xi) * xi);
	sum += t;
/* L50: */
    }

L60:
    r__1 = *b - 1.f;
    i__1 = 1 - n;
    sum = gamma_(&r__1) * gamr_(a) * pow_ri(x, &i__1) * xtoeps * sum;

/* NOW EVALUATE THE INFINITE SUM.     ----------------------------------- */

L70:
    istrt = 0;
    if (n < 1) {
	istrt = 1 - n;
    }
    xi = (real) istrt;

    r__1 = *a + 1.f - *b;
    factor = pow_ri(&c_b25, &n) * gamr_(&r__1) * pow_ri(x, &istrt);
    if (beps != 0.f) {
	factor = factor * beps * pi / sin(beps * pi);
    }

    pochai = poch_(a, &xi);
    r__1 = xi + 1.f;
    gamri1 = gamr_(&r__1);
    r__1 = aintb + xi;
    gamrni = gamr_(&r__1);
    r__1 = xi - beps;
    r__2 = xi + 1.f - beps;
    b0 = factor * poch_(a, &r__1) * gamrni * gamr_(&r__2);

    if ((r__1 = xtoeps - 1.f, dabs(r__1)) > .5f) {
	goto L90;
    }

/* X**(-BEPS) IS CLOSE TO 1.0, SO WE MUST BE CAREFUL IN EVALUATING */
/* THE DIFFERENCES */

    r__1 = *a + xi;
    r__2 = -beps;
    pch1ai = poch1_(&r__1, &r__2);
    r__1 = xi + 1.f - beps;
    pch1i = poch1_(&r__1, &beps);
    r__1 = *b + xi;
    r__2 = -beps;
    c0 = factor * pochai * gamrni * gamri1 * (-poch1_(&r__1, &r__2) + pch1ai 
	    - pch1i + beps * pch1ai * pch1i);

/* XEPS1 = (1.0 - X**(-BEPS)) / BEPS */
    r__1 = -beps * alnx;
    xeps1 = alnx * exprel_(&r__1);

    ret_val = sum + c0 + xeps1 * b0;
    xn = (real) n;
    for (i__ = 1; i__ <= 1000; ++i__) {
	xi = (real) (istrt + i__);
	xi1 = (real) (istrt + i__ - 1);
	b0 = (*a + xi1 - beps) * b0 * *x / ((xn + xi1) * (xi - beps));
	c0 = (*a + xi1) * c0 * *x / ((*b + xi1) * xi) - ((*a - 1.f) * (xn + 
		xi * 2.f - 1.f) + xi * (xi - beps)) * b0 / (xi * (*b + xi1) * 
		(*a + xi1 - beps));
	t = c0 + xeps1 * b0;
	ret_val += t;
	if (dabs(t) < eps * dabs(ret_val)) {
	    goto L130;
	}
/* L80: */
    }
    xermsg_("SLATEC", "CHU", "NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING "
	    "SERIES", &c__3, &c__2, (ftnlen)6, (ftnlen)3, (ftnlen)52);

/* X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD */
/* FORMULATION IS STABLE. */

L90:
    r__1 = *b + xi;
    a0 = factor * pochai * gamr_(&r__1) * gamri1 / beps;
    b0 = xtoeps * b0 / beps;

    ret_val = sum + a0 - b0;
    for (i__ = 1; i__ <= 1000; ++i__) {
	xi = (real) (istrt + i__);
	xi1 = (real) (istrt + i__ - 1);
	a0 = (*a + xi1) * a0 * *x / ((*b + xi1) * xi);
	b0 = (*a + xi1 - beps) * b0 * *x / ((aintb + xi1) * (xi - beps));
	t = a0 - b0;
	ret_val += t;
	if (dabs(t) < eps * dabs(ret_val)) {
	    goto L130;
	}
/* L100: */
    }
    xermsg_("SLATEC", "CHU", "NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING "
	    "SERIES", &c__3, &c__2, (ftnlen)6, (ftnlen)3, (ftnlen)52);

/* USE LUKE-S RATIONAL APPROX IN THE ASYMPTOTIC REGION. */

L120:
    d__1 = (doublereal) (*x);
    d__2 = (doublereal) (-(*a));
    ret_val = pow_dd(&d__1, &d__2) * r9chu_(a, b, x);

L130:
    return ret_val;
} /* chu_ */

