/* dpoch.f -- translated by f2c (version 12.02.01).
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
static doublereal c_b9 = -1.;

/* DECK DPOCH */
doublereal dpoch_(doublereal *a, doublereal *x)
{
    /* Initialized data */

    static doublereal pi = 3.141592653589793238462643383279503;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Local variables */
    static doublereal b;
    static integer i__, n, ia;
    static doublereal ax;
    extern doublereal dfac_(integer *);
    static doublereal absa;
    extern doublereal dcot_(doublereal *);
    static doublereal alnga;
    extern doublereal dgamr_(doublereal *);
    static doublereal absax, sgnga;
    extern doublereal d9lgmc_(doublereal *), dgamma_(doublereal *);
    extern /* Subroutine */ int dlgams_(doublereal *, doublereal *, 
	    doublereal *);
    static doublereal alngax;
    extern doublereal dlnrel_(doublereal *);
    static doublereal sgngax;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DPOCH */
/* ***PURPOSE  Evaluate a generalization of Pochhammer's symbol. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C1, C7A */
/* ***TYPE      DOUBLE PRECISION (POCH-S, DPOCH-D) */
/* ***KEYWORDS  FNLIB, POCHHAMMER, SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate a double precision generalization of Pochhammer's symbol */
/* (A)-sub-X = GAMMA(A+X)/GAMMA(A) for double precision A and X. */
/* For X a non-negative integer, POCH(A,X) is just Pochhammer's symbol. */
/* This is a preliminary version that does not handle wrong arguments */
/* properly and may not properly handle the case when the result is */
/* computed to less than half of double precision. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D9LGMC, DFAC, DGAMMA, DGAMR, DLGAMS, DLNREL, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900727  Added EXTERNAL statement.  (WRB) */
/* ***END PROLOGUE  DPOCH */
/* ***FIRST EXECUTABLE STATEMENT  DPOCH */
    ax = *a + *x;
    if (ax > 0.) {
	goto L30;
    }
    if (d_int(&ax) != ax) {
	goto L30;
    }

    if (*a > 0. || d_int(a) != *a) {
	xermsg_("SLATEC", "DPOCH", "A+X IS NON-POSITIVE INTEGER BUT A IS NOT",
		 &c__2, &c__2, (ftnlen)6, (ftnlen)5, (ftnlen)40);
    }

/* WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS. */

    ret_val = 1.;
    if (*x == 0.) {
	return ret_val;
    }

    n = (integer) (*x);
/* Computing MIN */
    d__1 = *a + *x;
    if (min(d__1,*a) < -20.) {
	goto L20;
    }

    ia = (integer) (*a);
    i__1 = -ia;
    i__2 = -ia - n;
    ret_val = pow_di(&c_b9, &n) * dfac_(&i__1) / dfac_(&i__2);
    return ret_val;

L20:
    d__1 = *x / (*a - 1.);
    d__2 = -(*a) + 1.;
    d__3 = -(*a) - *x + 1.;
    ret_val = pow_di(&c_b9, &n) * exp((*a - .5) * dlnrel_(&d__1) + *x * log(-(
	    *a) + 1. - *x) - *x + d9lgmc_(&d__2) - d9lgmc_(&d__3));
    return ret_val;

/* A+X IS NOT ZERO OR A NEGATIVE INTEGER. */

L30:
    ret_val = 0.;
    if (*a <= 0. && d_int(a) == *a) {
	return ret_val;
    }

    n = (integer) abs(*x);
    if ((doublereal) n != *x || n > 20) {
	goto L50;
    }

/* X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE. */

    ret_val = 1.;
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
    absax = (d__1 = *a + *x, abs(d__1));
    absa = abs(*a);
    if (max(absax,absa) > 20.) {
	goto L60;
    }
    d__1 = *a + *x;
    ret_val = dgamma_(&d__1) * dgamr_(a);
    return ret_val;

L60:
    if (abs(*x) > absa * .5) {
	goto L70;
    }

/* ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS, */
/* A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE */
/* GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) * */
/* SIN(PI*A)/SIN(PI*(A+X)) */

    b = *a;
    if (b < 0.) {
	b = -(*a) - *x + 1.;
    }
    d__1 = *x / b;
    d__2 = b + *x;
    ret_val = exp((b - .5) * dlnrel_(&d__1) + *x * log(b + *x) - *x + d9lgmc_(
	    &d__2) - d9lgmc_(&b));
    if (*a < 0. && ret_val != 0.) {
	d__1 = pi * *a;
	ret_val /= cos(pi * *x) + dcot_(&d__1) * sin(pi * *x);
    }
    return ret_val;

L70:
    d__1 = *a + *x;
    dlgams_(&d__1, &alngax, &sgngax);
    dlgams_(a, &alnga, &sgnga);
    ret_val = sgngax * sgnga * exp(alngax - alnga);

    return ret_val;
} /* dpoch_ */

