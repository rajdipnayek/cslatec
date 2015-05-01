/* xpqnu.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    integer nbitsf;
} xblk1_;

#define xblk1_1 xblk1_

/* Table of constant values */

static real c_b2 = 1.f;

/* DECK XPQNU */
/* Subroutine */ int xpqnu_(real *nu1, real *nu2, integer *mu, real *theta, 
	integer *id, real *pqa, integer *ipqa, integer *ierror)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static real a;
    static integer i__, j, k;
    static real r__, w, x, y, z__;
    static integer j0;
    static real x1, x2;
    static integer ia;
    static real di;
    static integer if__;
    static real pq, nu, xs, pq1, pq2;
    static integer ix1;
    static real dmu;
    static integer ipq, ixs, ipq1, ipq2;
    extern /* Subroutine */ int xadd_(real *, integer *, real *, integer *, 
	    real *, integer *, integer *), xadj_(real *, integer *, integer *)
	    ;
    static real flok;
    extern doublereal xpsi_(real *, integer *, integer *);
    static integer ipsik, ipsix;
    static real factmu;

/* ***BEGIN PROLOGUE  XPQNU */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute the values of Legendre functions for XLEGF. */
/*            This subroutine calculates initial values of P or Q using */
/*            power series, then performs forward nu-wise recurrence to */
/*            obtain P(-MU,NU,X), Q(0,NU,X), or Q(1,NU,X). The nu-wise */
/*            recurrence is stable for P for all mu and for Q for mu=0,1. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      SINGLE PRECISION (XPQNU-S, DXPQNU-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Smith, John M., (NBS and George Mason University) */
/* ***ROUTINES CALLED  XADD, XADJ, XPSI */
/* ***COMMON BLOCKS    XBLK1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820728  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  XPQNU */

/*        J0, IPSIK, AND IPSIX ARE INITIALIZED IN THIS SUBROUTINE. */
/*        J0 IS THE NUMBER OF TERMS USED IN SERIES EXPANSION */
/*        IN SUBROUTINE XPQNU. */
/*        IPSIK, IPSIX ARE VALUES OF K AND X RESPECTIVELY */
/*        USED IN THE CALCULATION OF THE XPSI FUNCTION. */

/* ***FIRST EXECUTABLE STATEMENT  XPQNU */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    j0 = xblk1_1.nbitsf;
    ipsik = xblk1_1.nbitsf / 10 + 1;
    ipsix = ipsik * 5;
    ipq = 0;
/*        FIND NU IN INTERVAL [-.5,.5) IF ID=2  ( CALCULATION OF Q ) */
    nu = r_mod(nu1, &c_b2);
    if (nu >= .5f) {
	nu += -1.f;
    }
/*        FIND NU IN INTERVAL (-1.5,-.5] IF ID=1,3, OR 4  ( CALC. OF P ) */
    if (*id != 2 && nu > -.5f) {
	nu += -1.f;
    }
/*        CALCULATE MU FACTORIAL */
    k = *mu;
    dmu = (real) (*mu);
    if (*mu <= 0) {
	goto L60;
    }
    factmu = 1.f;
    if__ = 0;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	factmu *= i__;
/* L50: */
	xadj_(&factmu, &if__, ierror);
    }
    if (*ierror != 0) {
	return 0;
    }
L60:
    if (k == 0) {
	factmu = 1.f;
    }
    if (k == 0) {
	if__ = 0;
    }

/*        X=COS(THETA) */
/*        Y=SIN(THETA/2)**2=(1-X)/2=.5-.5*X */
/*        R=TAN(THETA/2)=SQRT((1-X)/(1+X) */

    x = cos(*theta);
/* Computing 2nd power */
    r__1 = sin(*theta / 2.f);
    y = r__1 * r__1;
    r__ = tan(*theta / 2.f);

/*        USE ASCENDING SERIES TO CALCULATE TWO VALUES OF P OR Q */
/*        FOR USE AS STARTING VALUES IN RECURRENCE RELATION. */

    pq2 = 0.f;
    for (j = 1; j <= 2; ++j) {
	ipq1 = 0;
	if (*id == 2) {
	    goto L80;
	}

/*        SERIES FOR P ( ID = 1, 3, OR 4 ) */
/*        P(-MU,NU,X)=1./FACTORIAL(MU)*SQRT(((1.-X)/(1.+X))**MU) */
/*                *SUM(FROM 0 TO J0-1)A(J)*(.5-.5*X)**J */

	ipq = 0;
	pq = 1.f;
	a = 1.f;
	ia = 0;
	i__1 = j0;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    di = (real) i__;
	    a = a * y * (di - 2.f - nu) * (di - 1.f + nu) / ((di - 1.f + dmu) 
		    * (di - 1.f));
	    xadj_(&a, &ia, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
	    if (a == 0.f) {
		goto L66;
	    }
	    xadd_(&pq, &ipq, &a, &ia, &pq, &ipq, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
/* L65: */
	}
L66:
	if (*mu <= 0) {
	    goto L90;
	}
	x2 = r__;
	x1 = pq;
	k = *mu;
	i__1 = k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x1 *= x2;
/* L77: */
	    xadj_(&x1, &ipq, ierror);
	}
	if (*ierror != 0) {
	    return 0;
	}
	pq = x1 / factmu;
	ipq -= if__;
	xadj_(&pq, &ipq, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	goto L90;

/*        Z=-LN(R)=.5*LN((1+X)/(1-X)) */

L80:
	z__ = -log(r__);
	r__1 = nu + 1.f;
	w = xpsi_(&r__1, &ipsik, &ipsix);
	xs = 1.f / sin(*theta);

/*        SERIES SUMMATION FOR Q ( ID = 2 ) */
/*        Q(0,NU,X)=SUM(FROM 0 TO J0-1)((.5*LN((1+X)/(1-X)) */
/*    +XPSI(J+1,IPSIK,IPSIX)-XPSI(NU+1,IPSIK,IPSIX)))*A(J)*(.5-.5*X)**J */

/*        Q(1,NU,X)=-SQRT(1./(1.-X**2))+SQRT((1-X)/(1+X)) */
/*             *SUM(FROM 0 T0 J0-1)(-NU*(NU+1)/2*LN((1+X)/(1-X)) */
/*                 +(J-NU)*(J+NU+1)/(2*(J+1))+NU*(NU+1)* */
/*     (XPSI(NU+1,IPSIK,IPSIX)-XPSI(J+1,IPSIK,IPSIX))*A(J)*(.5-.5*X)**J */

/*        NOTE, IN THIS LOOP K=J+1 */

	pq = 0.f;
	ipq = 0;
	ia = 0;
	a = 1.f;
	i__1 = j0;
	for (k = 1; k <= i__1; ++k) {
	    flok = (real) k;
	    if (k == 1) {
		goto L81;
	    }
	    a = a * y * (flok - 2.f - nu) * (flok - 1.f + nu) / ((flok - 1.f 
		    + dmu) * (flok - 1.f));
	    xadj_(&a, &ia, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
L81:
	    if (*mu >= 1) {
		goto L83;
	    }
	    x1 = (xpsi_(&flok, &ipsik, &ipsix) - w + z__) * a;
	    ix1 = ia;
	    xadd_(&pq, &ipq, &x1, &ix1, &pq, &ipq, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
	    goto L85;
L83:
	    x1 = (nu * (nu + 1.f) * (z__ - w + xpsi_(&flok, &ipsik, &ipsix)) 
		    + (nu - flok + 1.f) * (nu + flok) / (flok * 2.f)) * a;
	    ix1 = ia;
	    xadd_(&pq, &ipq, &x1, &ix1, &pq, &ipq, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
L85:
	    ;
	}
	if (*mu >= 1) {
	    pq = -r__ * pq;
	}
	ixs = 0;
	if (*mu >= 1) {
	    r__1 = -xs;
	    xadd_(&pq, &ipq, &r__1, &ixs, &pq, &ipq, ierror);
	}
	if (*ierror != 0) {
	    return 0;
	}
	if (j == 2) {
	    *mu = -(*mu);
	}
	if (j == 2) {
	    dmu = -dmu;
	}
L90:
	if (j == 1) {
	    pq2 = pq;
	}
	if (j == 1) {
	    ipq2 = ipq;
	}
	nu += 1.f;
/* L100: */
    }
    k = 0;
    if (nu - 1.5f < *nu1) {
	goto L120;
    }
    ++k;
    pqa[k] = pq2;
    ipqa[k] = ipq2;
    if (nu > *nu2 + .5f) {
	return 0;
    }
L120:
    pq1 = pq;
    ipq1 = ipq;
    if (nu < *nu1 + .5f) {
	goto L130;
    }
    ++k;
    pqa[k] = pq;
    ipqa[k] = ipq;
    if (nu > *nu2 + .5f) {
	return 0;
    }

/*        FORWARD NU-WISE RECURRENCE FOR F(MU,NU,X) FOR FIXED MU */
/*        USING */
/*        (NU+MU+1)*F(MU,NU,X)=(2.*NU+1)*F(MU,NU,X)-(NU-MU)*F(MU,NU-1,X) */
/*        WHERE F(MU,NU,X) MAY BE P(-MU,NU,X) OR IF MU IS REPLACED */
/*        BY -MU THEN F(MU,NU,X) MAY BE Q(MU,NU,X). */
/*        NOTE, IN THIS LOOP, NU=NU+1 */

L130:
    x1 = (nu * 2.f - 1.f) / (nu + dmu) * x * pq1;
    x2 = (nu - 1.f - dmu) / (nu + dmu) * pq2;
    r__1 = -x2;
    xadd_(&x1, &ipq1, &r__1, &ipq2, &pq, &ipq, ierror);
    if (*ierror != 0) {
	return 0;
    }
    xadj_(&pq, &ipq, ierror);
    if (*ierror != 0) {
	return 0;
    }
    nu += 1.f;
    pq2 = pq1;
    ipq2 = ipq1;
    goto L120;

} /* xpqnu_ */

