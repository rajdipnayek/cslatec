/* xpmu.f -- translated by f2c (version 12.02.01).
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

/* DECK XPMU */
/* Subroutine */ int xpmu_(real *nu1, real *nu2, integer *mu1, integer *mu2, 
	real *theta, real *x, real *sx, integer *id, real *pqa, integer *ipqa,
	 integer *ierror)
{
    static integer j, n;
    static real p0, x1, x2;
    static integer mu, ip0;
    extern /* Subroutine */ int xadd_(real *, integer *, real *, integer *, 
	    real *, integer *, integer *), xadj_(real *, integer *, integer *)
	    , xpqnu_(real *, real *, integer *, real *, integer *, real *, 
	    integer *, integer *);

/* ***BEGIN PROLOGUE  XPMU */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute the values of Legendre functions for XLEGF. */
/*            Method: backward mu-wise recurrence for P(-MU,NU,X) for */
/*            fixed nu to obtain P(-MU2,NU1,X), P(-(MU2-1),NU1,X), ..., */
/*            P(-MU1,NU1,X) and store in ascending mu order. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      SINGLE PRECISION (XPMU-S, DXPMU-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Smith, John M., (NBS and George Mason University) */
/* ***ROUTINES CALLED  XADD, XADJ, XPQNU */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820728  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  XPMU */

/*        CALL XPQNU TO OBTAIN P(-MU2,NU,X) */

/* ***FIRST EXECUTABLE STATEMENT  XPMU */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    xpqnu_(nu1, nu2, mu2, theta, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }
    p0 = pqa[1];
    ip0 = ipqa[1];
    mu = *mu2 - 1;

/*        CALL XPQNU TO OBTAIN P(-MU2-1,NU,X) */

    xpqnu_(nu1, nu2, &mu, theta, id, &pqa[1], &ipqa[1], ierror);
    if (*ierror != 0) {
	return 0;
    }
    n = *mu2 - *mu1 + 1;
    pqa[n] = p0;
    ipqa[n] = ip0;
    if (n == 1) {
	goto L300;
    }
    pqa[n - 1] = pqa[1];
    ipqa[n - 1] = ipqa[1];
    if (n == 2) {
	goto L300;
    }
    j = n - 2;
L290:

/*        BACKWARD RECURRENCE IN MU TO OBTAIN */
/*              P(-MU2,NU1,X),P(-(MU2-1),NU1,X),....P(-MU1,NU1,X) */
/*              USING */
/*              (NU-MU)*(NU+MU+1.)*P(-(MU+1),NU,X)= */
/*                2.*MU*X*SQRT((1./(1.-X**2))*P(-MU,NU,X)-P(-(MU-1),NU,X) */

    x1 = mu * 2.f * *x * *sx * pqa[j + 1];
    x2 = -(*nu1 - mu) * (*nu1 + mu + 1.f) * pqa[j + 2];
    xadd_(&x1, &ipqa[j + 1], &x2, &ipqa[j + 2], &pqa[j], &ipqa[j], ierror);
    if (*ierror != 0) {
	return 0;
    }
    xadj_(&pqa[j], &ipqa[j], ierror);
    if (*ierror != 0) {
	return 0;
    }
    if (j == 1) {
	goto L300;
    }
    --j;
    --mu;
    goto L290;
L300:
    return 0;
} /* xpmu_ */

