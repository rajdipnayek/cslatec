/* d9chu.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK D9CHU */
doublereal d9chu_(doublereal *a, doublereal *b, doublereal *z__)
{
    /* Initialized data */

    static logical first = TRUE_;

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal c2, g1, g2, g3, aa[4], ab, bb[4], bp, ct1, ct2, ct3, 
	    d1z, sab, x2i1, eps, anbn, sqeps;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  D9CHU */
/* ***SUBSIDIARY */
/* ***PURPOSE  Evaluate for large Z  Z**A * U(A,B,Z) where U is the */
/*            logarithmic confluent hypergeometric function. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  C11 */
/* ***TYPE      DOUBLE PRECISION (R9CHU-S, D9CHU-D) */
/* ***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION, */
/*             SPECIAL FUNCTIONS */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/* Evaluate for large Z  Z**A * U(A,B,Z)  where U is the logarithmic */
/* confluent hypergeometric function.  A rational approximation due to Y. */
/* L. Luke is used.  When U is not in the asymptotic region, i.e., when A */
/* or B is large compared with Z, considerable significance loss occurs. */
/* A warning is provided when the computed result is less than half */
/* precision. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  D1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   770801  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900720  Routine changed from user-callable to subsidiary.  (WRB) */
/* ***END PROLOGUE  D9CHU */
/* ***FIRST EXECUTABLE STATEMENT  D9CHU */
    if (first) {
	eps = d1mach_(&c__4) * 4.;
	sqeps = sqrt(d1mach_(&c__4));
    }
    first = FALSE_;

    bp = *a + 1. - *b;
    ab = *a * bp;
    ct2 = (*z__ - ab) * 2.;
    sab = *a + bp;

    bb[0] = 1.;
    aa[0] = 1.;

    ct3 = sab + 1. + ab;
    bb[1] = *z__ * 2. / ct3 + 1.;
    aa[1] = ct2 / ct3 + 1.;

    anbn = ct3 + sab + 3.;
    ct1 = *z__ * 2. / anbn + 1.;
    bb[2] = ct1 * 6. * *z__ / ct3 + 1.;
    aa[2] = ab * 6. / anbn + 1. + ct1 * 3. * ct2 / ct3;

    for (i__ = 4; i__ <= 300; ++i__) {
	x2i1 = (doublereal) ((i__ << 1) - 3);
	ct1 = x2i1 / (x2i1 - 2.);
	anbn = anbn + x2i1 + sab;
	ct2 = (x2i1 - 1.) / anbn;
	c2 = x2i1 * ct2 - 1.;
	d1z = x2i1 * 2. * *z__ / anbn;

	ct3 = sab * ct2;
	g1 = d1z + ct1 * (c2 + ct3);
	g2 = d1z - c2;
	g3 = ct1 * (1. - ct3 - ct2 * 2.);

	bb[3] = g1 * bb[2] + g2 * bb[1] + g3 * bb[0];
	aa[3] = g1 * aa[2] + g2 * aa[1] + g3 * aa[0];
	if ((d__2 = aa[3] * bb[0] - aa[0] * bb[3], abs(d__2)) < eps * (d__1 = 
		bb[3] * bb[0], abs(d__1))) {
	    goto L40;
	}

/* IF OVERFLOWS OR UNDERFLOWS PROVE TO BE A PROBLEM, THE STATEMENTS */
/* BELOW COULD BE ALTERED TO INCORPORATE A DYNAMICALLY ADJUSTED SCALE */
/* FACTOR. */

	for (j = 1; j <= 3; ++j) {
	    aa[j - 1] = aa[j];
	    bb[j - 1] = bb[j];
/* L20: */
	}
/* L30: */
    }
    xermsg_("SLATEC", "D9CHU", "NO CONVERGENCE IN 300 TERMS", &c__2, &c__2, (
	    ftnlen)6, (ftnlen)5, (ftnlen)27);

L40:
    ret_val = aa[3] / bb[3];

    if (ret_val < sqeps || ret_val > 1. / sqeps) {
	xermsg_("SLATEC", "D9CHU", "ANSWER LT HALF PRECISION", &c__2, &c__1, (
		ftnlen)6, (ftnlen)5, (ftnlen)24);
    }

    return ret_val;
} /* d9chu_ */

