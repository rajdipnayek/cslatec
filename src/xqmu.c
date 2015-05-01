/* xqmu.f -- translated by f2c (version 12.02.01).
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

/* DECK XQMU */
/* Subroutine */ int xqmu_(real *nu1, real *nu2, integer *mu1, integer *mu2, 
	real *theta, real *x, real *sx, integer *id, real *pqa, integer *ipqa,
	 integer *ierror)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static integer k;
    static real x1, x2, pq;
    static integer mu;
    static real nu, pq1, pq2, dmu;
    static integer ipq, ipq1, ipq2;
    extern /* Subroutine */ int xadd_(real *, integer *, real *, integer *, 
	    real *, integer *, integer *), xadj_(real *, integer *, integer *)
	    , xpqnu_(real *, real *, integer *, real *, integer *, real *, 
	    integer *, integer *);

/* ***BEGIN PROLOGUE  XQMU */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute the values of Legendre functions for XLEGF. */
/*            Method: forward mu-wise recurrence for Q(MU,NU,X) for fixed */
/*            nu to obtain Q(MU1,NU,X), Q(MU1+1,NU,X), ..., Q(MU2,NU,X). */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      SINGLE PRECISION (XQMU-S, DXQMU-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Smith, John M., (NBS and George Mason University) */
/* ***ROUTINES CALLED  XADD, XADJ, XPQNU */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820728  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  XQMU */
/* ***FIRST EXECUTABLE STATEMENT  XQMU */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    mu = 0;

/*        CALL XPQNU TO OBTAIN Q(0.,NU1,X) */

    xpqnu_(nu1, nu2, &mu, theta, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }
    pq2 = pqa[1];
    ipq2 = ipqa[1];
    mu = 1;

/*        CALL XPQNU TO OBTAIN Q(1.,NU1,X) */

    xpqnu_(nu1, nu2, &mu, theta, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }
    nu = *nu1;
    k = 0;
    mu = 1;
    dmu = 1.f;
    pq1 = pqa[1];
    ipq1 = ipqa[1];
    if (*mu1 > 0) {
	goto L310;
    }
    ++k;
    pqa[k] = pq2;
    ipqa[k] = ipq2;
    if (*mu2 < 1) {
	goto L330;
    }
L310:
    if (*mu1 > 1) {
	goto L320;
    }
    ++k;
    pqa[k] = pq1;
    ipqa[k] = ipq1;
    if (*mu2 <= 1) {
	goto L330;
    }
L320:

/*        FORWARD RECURRENCE IN MU TO OBTAIN */
/*                  Q(MU1,NU,X),Q(MU1+1,NU,X),....,Q(MU2,NU,X) USING */
/*             Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X) */
/*                               -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X) */

    x1 = dmu * -2.f * *x * *sx * pq1;
    x2 = (nu + dmu) * (nu - dmu + 1.f) * pq2;
    r__1 = -x2;
    xadd_(&x1, &ipq1, &r__1, &ipq2, &pq, &ipq, ierror);
    if (*ierror != 0) {
	return 0;
    }
    xadj_(&pq, &ipq, ierror);
    if (*ierror != 0) {
	return 0;
    }
    pq2 = pq1;
    ipq2 = ipq1;
    pq1 = pq;
    ipq1 = ipq;
    ++mu;
    dmu += 1.f;
    if (mu < *mu1) {
	goto L320;
    }
    ++k;
    pqa[k] = pq;
    ipqa[k] = ipq;
    if (*mu2 > mu) {
	goto L320;
    }
L330:
    return 0;
} /* xqmu_ */

