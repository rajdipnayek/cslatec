/* dchu.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b25 = -1.;

/* DECK DCHU */
doublereal dchu_(doublereal *a, doublereal *b, doublereal *x)
{
    /* Initialized data */

    static doublereal pi = 3.141592653589793238462643383279503;
    static doublereal eps = 0.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static integer i__, m, n;
    static doublereal t, a0, b0, c0, xi, xn, xi1, sum, beps, alnx, pch1i;
    extern doublereal d9chu_(doublereal *, doublereal *, doublereal *);
    static doublereal xeps1;
    extern doublereal dgamr_(doublereal *);
    static doublereal aintb;
    extern doublereal dpoch_(doublereal *, doublereal *), d1mach_(integer *);
    static doublereal pch1ai;
    static integer istrt;
    extern doublereal dpoch1_(doublereal *, doublereal *);
    static doublereal gamri1;
    extern doublereal dgamma_(doublereal *);
    static doublereal pochai, gamrni, factor;
    extern doublereal dexprl_(doublereal *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal xtoeps;

/* ***BEGIN PROLOGUE  DCHU */
/* ***PURPOSE  Compute the logarithmic confluent hypergeometric function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C11 */
/* ***TYPE      DOUBLE PRECISION (CHU-S, DCHU-D) */
/* ***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* DCHU(A,B,X) calculates the double precision logarithmic confluent */
/* hypergeometric function U(A,B,X) for double precision arguments */
/* A, B, and X. */

/* This routine is not valid when 1+A-B is close to zero if X is small. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, D9CHU, DEXPRL, DGAMMA, DGAMR, DPOCH, */
/*                    DPOCH1, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  DCHU */
/* ***FIRST EXECUTABLE STATEMENT  DCHU */
    if (eps == 0.) {
	eps = d1mach_(&c__3);
    }

    if (*x == 0.) {
	xermsg_("SLATEC", "DCHU", "X IS ZERO SO DCHU IS INFINITE", &c__1, &
		c__2, (ftnlen)6, (ftnlen)4, (ftnlen)29);
    }
    if (*x < 0.) {
	xermsg_("SLATEC", "DCHU", "X IS NEGATIVE, USE CCHU", &c__2, &c__2, (
		ftnlen)6, (ftnlen)4, (ftnlen)23);
    }

/* Computing MAX */
    d__2 = abs(*a);
/* Computing MAX */
    d__3 = (d__1 = *a + 1. - *b, abs(d__1));
    if (max(d__2,1.) * max(d__3,1.) < abs(*x) * .99) {
	goto L120;
    }

/* THE ASCENDING SERIES WILL BE USED, BECAUSE THE DESCENDING RATIONAL */
/* APPROXIMATION (WHICH IS BASED ON THE ASYMPTOTIC SERIES) IS UNSTABLE. */

    if ((d__1 = *a + 1. - *b, abs(d__1)) < sqrt(eps)) {
	xermsg_("SLATEC", "DCHU", "ALGORITHMIS BAD WHEN 1+A-B IS NEAR ZERO F"
		"OR SMALL X", &c__10, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)51);
    }

    if (*b >= 0.) {
	d__1 = *b + .5;
	aintb = d_int(&d__1);
    }
    if (*b < 0.) {
	d__1 = *b - .5;
	aintb = d_int(&d__1);
    }
    beps = *b - aintb;
    n = (integer) aintb;

    alnx = log(*x);
    xtoeps = exp(-beps * alnx);

/* EVALUATE THE FINITE SUM.     ----------------------------------------- */

    if (n >= 1) {
	goto L40;
    }

/* CONSIDER THE CASE B .LT. 1.0 FIRST. */

    sum = 1.;
    if (n == 0) {
	goto L30;
    }

    t = 1.;
    m = -n;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi1 = (doublereal) (i__ - 1);
	t = t * (*a + xi1) * *x / ((*b + xi1) * (xi1 + 1.));
	sum += t;
/* L20: */
    }

L30:
    d__1 = *a + 1. - *b;
    d__2 = -(*a);
    sum = dpoch_(&d__1, &d__2) * sum;
    goto L70;

/* NOW CONSIDER THE CASE B .GE. 1.0. */

L40:
    sum = 0.;
    m = n - 2;
    if (m < 0) {
	goto L70;
    }
    t = 1.;
    sum = 1.;
    if (m == 0) {
	goto L60;
    }

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xi = (doublereal) i__;
	t = t * (*a - *b + xi) * *x / ((1. - *b + xi) * xi);
	sum += t;
/* L50: */
    }

L60:
    d__1 = *b - 1.;
    i__1 = 1 - n;
    sum = dgamma_(&d__1) * dgamr_(a) * pow_di(x, &i__1) * xtoeps * sum;

/* NEXT EVALUATE THE INFINITE SUM.     ---------------------------------- */

L70:
    istrt = 0;
    if (n < 1) {
	istrt = 1 - n;
    }
    xi = (doublereal) istrt;

    d__1 = *a + 1. - *b;
    factor = pow_di(&c_b25, &n) * dgamr_(&d__1) * pow_di(x, &istrt);
    if (beps != 0.) {
	factor = factor * beps * pi / sin(beps * pi);
    }

    pochai = dpoch_(a, &xi);
    d__1 = xi + 1.;
    gamri1 = dgamr_(&d__1);
    d__1 = aintb + xi;
    gamrni = dgamr_(&d__1);
    d__1 = xi - beps;
    d__2 = xi + 1. - beps;
    b0 = factor * dpoch_(a, &d__1) * gamrni * dgamr_(&d__2);

    if ((d__1 = xtoeps - 1., abs(d__1)) > .5) {
	goto L90;
    }

/* X**(-BEPS) IS CLOSE TO 1.0D0, SO WE MUST BE CAREFUL IN EVALUATING THE */
/* DIFFERENCES. */

    d__1 = *a + xi;
    d__2 = -beps;
    pch1ai = dpoch1_(&d__1, &d__2);
    d__1 = xi + 1. - beps;
    pch1i = dpoch1_(&d__1, &beps);
    d__1 = *b + xi;
    d__2 = -beps;
    c0 = factor * pochai * gamrni * gamri1 * (-dpoch1_(&d__1, &d__2) + pch1ai 
	    - pch1i + beps * pch1ai * pch1i);

/* XEPS1 = (1.0 - X**(-BEPS))/BEPS = (X**(-BEPS) - 1.0)/(-BEPS) */
    d__1 = -beps * alnx;
    xeps1 = alnx * dexprl_(&d__1);

    ret_val = sum + c0 + xeps1 * b0;
    xn = (doublereal) n;
    for (i__ = 1; i__ <= 1000; ++i__) {
	xi = (doublereal) (istrt + i__);
	xi1 = (doublereal) (istrt + i__ - 1);
	b0 = (*a + xi1 - beps) * b0 * *x / ((xn + xi1) * (xi - beps));
	c0 = (*a + xi1) * c0 * *x / ((*b + xi1) * xi) - ((*a - 1.) * (xn + xi 
		* 2. - 1.) + xi * (xi - beps)) * b0 / (xi * (*b + xi1) * (*a 
		+ xi1 - beps));
	t = c0 + xeps1 * b0;
	ret_val += t;
	if (abs(t) < eps * abs(ret_val)) {
	    goto L130;
	}
/* L80: */
    }
    xermsg_("SLATEC", "DCHU", "NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING"
	    " SERIES", &c__3, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)52);

/* X**(-BEPS) IS VERY DIFFERENT FROM 1.0, SO THE STRAIGHTFORWARD */
/* FORMULATION IS STABLE. */

L90:
    d__1 = *b + xi;
    a0 = factor * pochai * dgamr_(&d__1) * gamri1 / beps;
    b0 = xtoeps * b0 / beps;

    ret_val = sum + a0 - b0;
    for (i__ = 1; i__ <= 1000; ++i__) {
	xi = (doublereal) (istrt + i__);
	xi1 = (doublereal) (istrt + i__ - 1);
	a0 = (*a + xi1) * a0 * *x / ((*b + xi1) * xi);
	b0 = (*a + xi1 - beps) * b0 * *x / ((aintb + xi1) * (xi - beps));
	t = a0 - b0;
	ret_val += t;
	if (abs(t) < eps * abs(ret_val)) {
	    goto L130;
	}
/* L100: */
    }
    xermsg_("SLATEC", "DCHU", "NO CONVERGENCE IN 1000 TERMS OF THE ASCENDING"
	    " SERIES", &c__3, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)52);

/* USE LUKE-S RATIONAL APPROXIMATION IN THE ASYMPTOTIC REGION. */

L120:
    d__1 = -(*a);
    ret_val = pow_dd(x, &d__1) * d9chu_(a, b, x);

L130:
    return ret_val;
} /* dchu_ */

