/* dxpnrm.f -- translated by f2c (version 12.02.01).
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

/* DECK DXPNRM */
/* Subroutine */ int dxpnrm_(doublereal *nu1, doublereal *nu2, integer *mu1, 
	integer *mu2, doublereal *pqa, integer *ipqa, integer *ierror)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal c1;
    static integer mu;
    static doublereal nu, dmu, prod;
    extern /* Subroutine */ int dxadj_(doublereal *, integer *, integer *);
    static integer iprod;

/* ***BEGIN PROLOGUE  DXPNRM */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute the values of Legendre functions for DXLEGF. */
/*            This subroutine transforms an array of Legendre functions */
/*            of the first kind of negative order stored in array PQA */
/*            into normalized Legendre polynomials stored in array PQA. */
/*            The original array is destroyed. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      DOUBLE PRECISION (XPNRM-S, DXPNRM-D) */
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
/* ***END PROLOGUE  DXPNRM */
/* ***FIRST EXECUTABLE STATEMENT  DXPNRM */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    l = (integer) (*mu2 - *mu1 + (*nu2 - *nu1 + 1.5));
    mu = *mu1;
    dmu = (doublereal) (*mu1);
    nu = *nu1;

/*         IF MU .GT.NU, NORM P =0. */

    j = 1;
L500:
    if (dmu <= nu) {
	goto L505;
    }
    pqa[j] = 0.;
    ipqa[j] = 0;
    ++j;
    if (j > l) {
	return 0;
    }

/*        INCREMENT EITHER MU OR NU AS APPROPRIATE. */

    if (*mu2 > *mu1) {
	dmu += 1.;
    }
    if (*nu2 - *nu1 > .5) {
	nu += 1.;
    }
    goto L500;

/*         TRANSFORM P(-MU,NU,X) INTO NORMALIZED P(MU,NU,X) USING */
/*              NORM P(MU,NU,X)= */
/*                 SQRT((NU+.5)*FACTORIAL(NU+MU)/FACTORIAL(NU-MU)) */
/*                              *P(-MU,NU,X) */

L505:
    prod = 1.;
    iprod = 0;
    k = mu << 1;
    if (k <= 0) {
	goto L520;
    }
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	prod *= sqrt(nu + dmu + 1. - i__);
/* L510: */
	dxadj_(&prod, &iprod, ierror);
    }
    if (*ierror != 0) {
	return 0;
    }
L520:
    i__1 = l;
    for (i__ = j; i__ <= i__1; ++i__) {
	c1 = prod * sqrt(nu + .5);
	pqa[i__] *= c1;
	ipqa[i__] += iprod;
	dxadj_(&pqa[i__], &ipqa[i__], ierror);
	if (*ierror != 0) {
	    return 0;
	}
	if (*nu2 - *nu1 > .5) {
	    goto L530;
	}
	if (dmu >= nu) {
	    goto L525;
	}
	prod = sqrt(nu + dmu + 1.) * prod;
	if (nu > dmu) {
	    prod *= sqrt(nu - dmu);
	}
	dxadj_(&prod, &iprod, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	++mu;
	dmu += 1.;
	goto L540;
L525:
	prod = 0.;
	iprod = 0;
	++mu;
	dmu += 1.;
	goto L540;
L530:
	prod = sqrt(nu + dmu + 1.) * prod;
	if (nu != dmu - 1.) {
	    prod /= sqrt(nu - dmu + 1.);
	}
	dxadj_(&prod, &iprod, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	nu += 1.;
L540:
	;
    }
    return 0;
} /* dxpnrm_ */

