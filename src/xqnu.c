/* xqnu.f -- translated by f2c (version 12.02.01).
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

/* DECK XQNU */
/* Subroutine */ int xqnu_(real *nu1, real *nu2, integer *mu1, real *theta, 
	real *x, real *sx, integer *id, real *pqa, integer *ipqa, integer *
	ierror)
{
    /* System generated locals */
    real r__1;

    /* Local variables */
    static integer k;
    static real x1, x2, pq;
    static integer mu;
    static real nu, pq1, pq2, dmu;
    static integer ipq, ipq1, ipq2;
    static real pql1, pql2;
    extern /* Subroutine */ int xadd_(real *, integer *, real *, integer *, 
	    real *, integer *, integer *), xadj_(real *, integer *, integer *)
	    ;
    static integer ipql1, ipql2;
    extern /* Subroutine */ int xpqnu_(real *, real *, integer *, real *, 
	    integer *, real *, integer *, integer *);

/* ***BEGIN PROLOGUE  XQNU */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute the values of Legendre functions for XLEGF. */
/*            Method: backward nu-wise recurrence for Q(MU,NU,X) for */
/*            fixed mu to obtain Q(MU1,NU1,X), Q(MU1,NU1+1,X), ..., */
/*            Q(MU1,NU2,X). */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      SINGLE PRECISION (XQNU-S, DXQNU-D) */
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
/* ***END PROLOGUE  XQNU */
/* ***FIRST EXECUTABLE STATEMENT  XQNU */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    k = 0;
    pq2 = 0.f;
    ipq2 = 0;
    pql2 = 0.f;
    ipql2 = 0;
    if (*mu1 == 1) {
	goto L290;
    }
    mu = 0;

/*        CALL XPQNU TO OBTAIN Q(0.,NU2,X) AND Q(0.,NU2-1,X) */

    xpqnu_(nu1, nu2, &mu, theta, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }
    if (*mu1 == 0) {
	return 0;
    }
    k = *nu2 - *nu1 + 1.5f;
    pq2 = pqa[k];
    ipq2 = ipqa[k];
    pql2 = pqa[k - 1];
    ipql2 = ipqa[k - 1];
L290:
    mu = 1;

/*        CALL XPQNU TO OBTAIN Q(1.,NU2,X) AND Q(1.,NU2-1,X) */

    xpqnu_(nu1, nu2, &mu, theta, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }
    if (*mu1 == 1) {
	return 0;
    }
    nu = *nu2;
    pq1 = pqa[k];
    ipq1 = ipqa[k];
    pql1 = pqa[k - 1];
    ipql1 = ipqa[k - 1];
L300:
    mu = 1;
    dmu = 1.f;
L320:

/*        FORWARD RECURRENCE IN MU TO OBTAIN Q(MU1,NU2,X) AND */
/*              Q(MU1,NU2-1,X) USING */
/*              Q(MU+1,NU,X)=-2.*MU*X*SQRT(1./(1.-X**2))*Q(MU,NU,X) */
/*                   -(NU+MU)*(NU-MU+1.)*Q(MU-1,NU,X) */

/*              FIRST FOR NU=NU2 */

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
    pqa[k] = pq;
    ipqa[k] = ipq;
    if (k == 1) {
	return 0;
    }
    if (nu < *nu2) {
	goto L340;
    }

/*              THEN FOR NU=NU2-1 */

    nu += -1.f;
    pq2 = pql2;
    ipq2 = ipql2;
    pq1 = pql1;
    ipq1 = ipql1;
    --k;
    goto L300;

/*         BACKWARD RECURRENCE IN NU TO OBTAIN */
/*              Q(MU1,NU1,X),Q(MU1,NU1+1,X),....,Q(MU1,NU2,X) */
/*              USING */
/*              (NU-MU+1.)*Q(MU,NU+1,X)= */
/*                       (2.*NU+1.)*X*Q(MU,NU,X)-(NU+MU)*Q(MU,NU-1,X) */

L340:
    pq1 = pqa[k];
    ipq1 = ipqa[k];
    pq2 = pqa[k + 1];
    ipq2 = ipqa[k + 1];
L350:
    if (nu <= *nu1) {
	return 0;
    }
    --k;
    x1 = (nu * 2.f + 1.f) * *x * pq1 / (nu + dmu);
    x2 = -(nu - dmu + 1.f) * pq2 / (nu + dmu);
    xadd_(&x1, &ipq1, &x2, &ipq2, &pq, &ipq, ierror);
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
    pqa[k] = pq;
    ipqa[k] = ipq;
    nu += -1.f;
    goto L350;
} /* xqnu_ */

