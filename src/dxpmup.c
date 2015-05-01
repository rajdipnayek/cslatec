/* dxpmup.f -- translated by f2c (version 12.02.01).
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

static real c_b2 = 1.f;
static integer c_n1 = -1;

/* DECK DXPMUP */
/* Subroutine */ int dxpmup_(doublereal *nu1, doublereal *nu2, integer *mu1, 
	integer *mu2, doublereal *pqa, integer *ipqa, integer *ierror)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__, j, k, l, n, mu;
    static doublereal nu, dmu, prod;
    extern /* Subroutine */ int dxadj_(doublereal *, integer *, integer *);
    static integer iprod;

/* ***BEGIN PROLOGUE  DXPMUP */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute the values of Legendre functions for DXLEGF. */
/*            This subroutine transforms an array of Legendre functions */
/*            of the first kind of negative order stored in array PQA */
/*            into Legendre functions of the first kind of positive */
/*            order stored in array PQA. The original array is destroyed. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      DOUBLE PRECISION (XPMUP-S, DXPMUP-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Smith, John M., (NBS and George Mason University) */
/* ***ROUTINES CALLED  DXADJ */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820728  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXPMUP */
/* ***FIRST EXECUTABLE STATEMENT  DXPMUP */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    nu = *nu1;
    mu = *mu1;
    dmu = (doublereal) mu;
    n = (integer) (*nu2 - *nu1 + .1) + (*mu2 - *mu1) + 1;
    j = 1;
    r__1 = (real) nu;
    if (r_mod(&r__1, &c_b2) != 0.f) {
	goto L210;
    }
L200:
    if (dmu < nu + 1.) {
	goto L210;
    }
    pqa[j] = 0.;
    ipqa[j] = 0;
    ++j;
    if (j > n) {
	return 0;
    }
/*        INCREMENT EITHER MU OR NU AS APPROPRIATE. */
    if (*nu2 - *nu1 > .5) {
	nu += 1.;
    }
    if (*mu2 > *mu1) {
	++mu;
    }
    goto L200;

/*        TRANSFORM P(-MU,NU,X) TO P(MU,NU,X) USING */
/*        P(MU,NU,X)=(NU-MU+1)*(NU-MU+2)*...*(NU+MU)*P(-MU,NU,X)*(-1)**MU */

L210:
    prod = 1.;
    iprod = 0;
    k = mu << 1;
    if (k == 0) {
	goto L222;
    }
    i__1 = k;
    for (l = 1; l <= i__1; ++l) {
	prod *= dmu - nu - l;
/* L220: */
	dxadj_(&prod, &iprod, ierror);
    }
    if (*ierror != 0) {
	return 0;
    }
L222:
    i__1 = n;
    for (i__ = j; i__ <= i__1; ++i__) {
	if (mu == 0) {
	    goto L225;
	}
	pqa[i__] = pqa[i__] * prod * pow_ii(&c_n1, &mu);
	ipqa[i__] += iprod;
	dxadj_(&pqa[i__], &ipqa[i__], ierror);
	if (*ierror != 0) {
	    return 0;
	}
L225:
	if (*nu2 - *nu1 > .5) {
	    goto L230;
	}
	prod = (dmu - nu) * prod * (-dmu - nu - 1.);
	dxadj_(&prod, &iprod, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	++mu;
	dmu += 1.;
	goto L240;
L230:
	prod = prod * (-dmu - nu - 1.) / (dmu - nu - 1.);
	dxadj_(&prod, &iprod, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	nu += 1.;
L240:
	;
    }
    return 0;
} /* dxpmup_ */

