/* poch.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static real c_b9 = -1.f;

/* DECK POCH */
doublereal poch_(real *a, real *x)
{
    /* Initialized data */

    static real pi = 3.141592653589793238f;

    /* System generated locals */
    integer i__1, i__2;
    real ret_val, r__1, r__2, r__3;

    /* Local variables */
    static real b;
    static integer i__, n;
    static real ax;
    extern doublereal fac_(integer *), cot_(real *);
    static real absa;
    extern doublereal gamr_(real *), gamma_(real *);
    static real alnga, absax, sgnga;
    extern doublereal r9lgmc_(real *);
    extern /* Subroutine */ int algams_(real *, real *, real *);
    static real alngax;
    extern doublereal alnrel_(real *);
    static real sgngax;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  POCH */
/* ***PURPOSE  Evaluate a generalization of Pochhammer's symbol. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1, C7A */
/* ***TYPE      SINGLE PRECISION (POCH-S, DPOCH-D) */
/* ***KEYWORDS  FNLIB, POCHHAMMER, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate a generalization of Pochhammer's symbol */
/* (A)-sub-X = GAMMA(A+X)/GAMMA(A).  For X a non-negative integer, */
/* POCH(A,X) is just Pochhammer's symbol.  A and X are single precision. */
/* This is a preliminary version.  Error handling when POCH(A,X) is */
/* less than half precision is probably incorrect.  Grossly incorrect */
/* arguments are not handled properly. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  ALGAMS, ALNREL, FAC, GAMMA, GAMR, R9LGMC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  POCH */
/* ***FIRST EXECUTABLE STATEMENT  POCH */
    ax = *a + *x;
    if (ax > 0.f) {
	goto L30;
    }
    if (r_int(&ax) != ax) {
	goto L30;
    }

    if (*a > 0.f || r_int(a) != *a) {
	xermsg_("SLATEC", "POCH", "A+X IS NON-POSITIVE INTEGER BUT A IS NOT", 
		&c__2, &c__2, (ftnlen)6, (ftnlen)4, (ftnlen)40);
    }

/* WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS. */

    ret_val = 1.f;
    if (*x == 0.f) {
	return ret_val;
    }

    n = *x;
/* Computing MIN */
    r__1 = *a + *x;
    if (dmin(r__1,*a) < -20.f) {
	goto L20;
    }

    i__1 = -((integer) (*a));
    i__2 = -((integer) (*a)) - n;
    ret_val = pow_ri(&c_b9, &n) * fac_(&i__1) / fac_(&i__2);
    return ret_val;

L20:
    r__1 = *x / (*a - 1.f);
    r__2 = -(*a) + 1.f;
    r__3 = -(*a) - *x + 1.f;
    ret_val = pow_ri(&c_b9, &n) * exp((*a - .5f) * alnrel_(&r__1) + *x * log(
	    -(*a) + 1.f - *x) - *x + r9lgmc_(&r__2) - r9lgmc_(&r__3));
    return ret_val;

/* HERE WE KNOW A+X IS NOT ZERO OR A NEGATIVE INTEGER. */

L30:
    ret_val = 0.f;
    if (*a <= 0.f && r_int(a) == *a) {
	return ret_val;
    }

    n = dabs(*x);
    if ((real) n != *x || n > 20) {
	goto L50;
    }

/* X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE. */

    ret_val = 1.f;
    if (n == 0) {
	return ret_val;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val *= *a + i__ - 1;
/* L40: */
    }
    return ret_val;

L50:
    absax = (r__1 = *a + *x, dabs(r__1));
    absa = dabs(*a);
    if (dmax(absax,absa) > 20.f) {
	goto L60;
    }
    r__1 = *a + *x;
    ret_val = gamma_(&r__1) * gamr_(a);
    return ret_val;

L60:
    if (dabs(*x) > absa * .5f) {
	goto L70;
    }

/* HERE ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS, */
/* A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE */
/* GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) * */
/* SIN(PI*A)/SIN(PI*(A+X)) */

    b = *a;
    if (b < 0.f) {
	b = -(*a) - *x + 1.f;
    }
    r__1 = *x / b;
    r__2 = b + *x;
    ret_val = exp((b - .5f) * alnrel_(&r__1) + *x * log(b + *x) - *x + 
	    r9lgmc_(&r__2) - r9lgmc_(&b));
    if (*a < 0.f && ret_val != 0.f) {
	r__1 = pi * *a;
	ret_val /= cos(pi * *x) + cot_(&r__1) * sin(pi * *x);
    }
    return ret_val;

L70:
    r__1 = *a + *x;
    algams_(&r__1, &alngax, &sgngax);
    algams_(a, &alnga, &sgnga);
    ret_val = sgngax * sgnga * exp(alngax - alnga);

    return ret_val;
} /* poch_ */

