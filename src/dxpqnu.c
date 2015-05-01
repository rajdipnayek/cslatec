/* dxpqnu.f -- translated by f2c (version 12.02.01).
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
} dxblk1_;

#define dxblk1_1 dxblk1_

/* Table of constant values */

static doublereal c_b2 = 1.;

/* DECK DXPQNU */
/* Subroutine */ int dxpqnu_(doublereal *nu1, doublereal *nu2, integer *mu, 
	doublereal *theta, integer *id, doublereal *pqa, integer *ipqa, 
	integer *ierror)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal a;
    static integer i__, j, k;
    static doublereal r__, w, x, y, z__;
    static integer j0;
    static doublereal x1, x2;
    static integer ia;
    static doublereal di;
    static integer if__;
    static doublereal pq, nu, xs, pq1, pq2;
    static integer ix1;
    static doublereal dmu;
    static integer ipq, ixs, ipq1, ipq2;
    static doublereal flok;
    extern /* Subroutine */ int dxadd_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), dxadj_(doublereal 
	    *, integer *, integer *);
    static integer ipsik;
    extern doublereal dxpsi_(doublereal *, integer *, integer *);
    static integer ipsix;
    static doublereal factmu;

/* ***BEGIN PROLOGUE  DXPQNU */
/* ***SUBSIDIARY */
/* ***PURPOSE  To compute the values of Legendre functions for DXLEGF. */
/*            This subroutine calculates initial values of P or Q using */
/*            power series, then performs forward nu-wise recurrence to */
/*            obtain P(-MU,NU,X), Q(0,NU,X), or Q(1,NU,X). The nu-wise */
/*            recurrence is stable for P for all mu and for Q for mu=0,1. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  C3A2, C9 */
/* ***TYPE      DOUBLE PRECISION (XPQNU-S, DXPQNU-D) */
/* ***KEYWORDS  LEGENDRE FUNCTIONS */
/* ***AUTHOR  Smith, John M., (NBS and George Mason University) */
/* ***ROUTINES CALLED  DXADD, DXADJ, DXPSI */
/* ***COMMON BLOCKS    DXBLK1 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820728  DATE WRITTEN */
/*   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS) */
/*   901019  Revisions to prologue.  (DWL and WRB) */
/*   901106  Changed all specific intrinsics to generic.  (WRB) */
/*           Corrected order of sections in prologue and added TYPE */
/*           section.  (WRB) */
/*   920127  Revised PURPOSE section of prologue.  (DWL) */
/* ***END PROLOGUE  DXPQNU */

/*        J0, IPSIK, AND IPSIX ARE INITIALIZED IN THIS SUBROUTINE. */
/*        J0 IS THE NUMBER OF TERMS USED IN SERIES EXPANSION */
/*        IN SUBROUTINE DXPQNU. */
/*        IPSIK, IPSIX ARE VALUES OF K AND X RESPECTIVELY */
/*        USED IN THE CALCULATION OF THE DXPSI FUNCTION. */

/* ***FIRST EXECUTABLE STATEMENT  DXPQNU */
    /* Parameter adjustments */
    --ipqa;
    --pqa;

    /* Function Body */
    *ierror = 0;
    j0 = dxblk1_1.nbitsf;
    ipsik = dxblk1_1.nbitsf / 10 + 1;
    ipsix = ipsik * 5;
    ipq = 0;
/*        FIND NU IN INTERVAL [-.5,.5) IF ID=2  ( CALCULATION OF Q ) */
    nu = d_mod(nu1, &c_b2);
    if (nu >= .5) {
	nu += -1.;
    }
/*        FIND NU IN INTERVAL (-1.5,-.5] IF ID=1,3, OR 4  ( CALC. OF P ) */
    if (*id != 2 && nu > -.5) {
	nu += -1.;
    }
/*        CALCULATE MU FACTORIAL */
    k = *mu;
    dmu = (doublereal) (*mu);
    if (*mu <= 0) {
	goto L60;
    }
    factmu = 1.;
    if__ = 0;
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	factmu *= i__;
/* L50: */
	dxadj_(&factmu, &if__, ierror);
    }
    if (*ierror != 0) {
	return 0;
    }
L60:
    if (k == 0) {
	factmu = 1.;
    }
    if (k == 0) {
	if__ = 0;
    }

/*        X=COS(THETA) */
/*        Y=SIN(THETA/2)**2=(1-X)/2=.5-.5*X */
/*        R=TAN(THETA/2)=SQRT((1-X)/(1+X) */

    x = cos(*theta);
/* Computing 2nd power */
    d__1 = sin(*theta / 2.);
    y = d__1 * d__1;
    r__ = tan(*theta / 2.);

/*        USE ASCENDING SERIES TO CALCULATE TWO VALUES OF P OR Q */
/*        FOR USE AS STARTING VALUES IN RECURRENCE RELATION. */

    pq2 = 0.;
    for (j = 1; j <= 2; ++j) {
	ipq1 = 0;
	if (*id == 2) {
	    goto L80;
	}

/*        SERIES FOR P ( ID = 1, 3, OR 4 ) */
/*        P(-MU,NU,X)=1./FACTORIAL(MU)*SQRT(((1.-X)/(1.+X))**MU) */
/*                *SUM(FROM 0 TO J0-1)A(J)*(.5-.5*X)**J */

	ipq = 0;
	pq = 1.;
	a = 1.;
	ia = 0;
	i__1 = j0;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    di = (doublereal) i__;
	    a = a * y * (di - 2. - nu) * (di - 1. + nu) / ((di - 1. + dmu) * (
		    di - 1.));
	    dxadj_(&a, &ia, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
	    if (a == 0.) {
		goto L66;
	    }
	    dxadd_(&pq, &ipq, &a, &ia, &pq, &ipq, ierror);
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
	    dxadj_(&x1, &ipq, ierror);
	}
	if (*ierror != 0) {
	    return 0;
	}
	pq = x1 / factmu;
	ipq -= if__;
	dxadj_(&pq, &ipq, ierror);
	if (*ierror != 0) {
	    return 0;
	}
	goto L90;

/*        Z=-LN(R)=.5*LN((1+X)/(1-X)) */

L80:
	z__ = -log(r__);
	d__1 = nu + 1.;
	w = dxpsi_(&d__1, &ipsik, &ipsix);
	xs = 1. / sin(*theta);

/*        SERIES SUMMATION FOR Q ( ID = 2 ) */
/*        Q(0,NU,X)=SUM(FROM 0 TO J0-1)((.5*LN((1+X)/(1-X)) */
/*    +DXPSI(J+1,IPSIK,IPSIX)-DXPSI(NU+1,IPSIK,IPSIX)))*A(J)*(.5-.5*X)**J */

/*        Q(1,NU,X)=-SQRT(1./(1.-X**2))+SQRT((1-X)/(1+X)) */
/*             *SUM(FROM 0 T0 J0-1)(-NU*(NU+1)/2*LN((1+X)/(1-X)) */
/*                 +(J-NU)*(J+NU+1)/(2*(J+1))+NU*(NU+1)* */
/*     (DXPSI(NU+1,IPSIK,IPSIX)-DXPSI(J+1,IPSIK,IPSIX))*A(J)*(.5-.5*X)**J */

/*        NOTE, IN THIS LOOP K=J+1 */

	pq = 0.;
	ipq = 0;
	ia = 0;
	a = 1.;
	i__1 = j0;
	for (k = 1; k <= i__1; ++k) {
	    flok = (doublereal) k;
	    if (k == 1) {
		goto L81;
	    }
	    a = a * y * (flok - 2. - nu) * (flok - 1. + nu) / ((flok - 1. + 
		    dmu) * (flok - 1.));
	    dxadj_(&a, &ia, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
L81:
	    if (*mu >= 1) {
		goto L83;
	    }
	    x1 = (dxpsi_(&flok, &ipsik, &ipsix) - w + z__) * a;
	    ix1 = ia;
	    dxadd_(&pq, &ipq, &x1, &ix1, &pq, &ipq, ierror);
	    if (*ierror != 0) {
		return 0;
	    }
	    goto L85;
L83:
	    x1 = (nu * (nu + 1.) * (z__ - w + dxpsi_(&flok, &ipsik, &ipsix)) 
		    + (nu - flok + 1.) * (nu + flok) / (flok * 2.)) * a;
	    ix1 = ia;
	    dxadd_(&pq, &ipq, &x1, &ix1, &pq, &ipq, ierror);
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
	    d__1 = -xs;
	    dxadd_(&pq, &ipq, &d__1, &ixs, &pq, &ipq, ierror);
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
	nu += 1.;
/* L100: */
    }
    k = 0;
    if (nu - 1.5 < *nu1) {
	goto L120;
    }
    ++k;
    pqa[k] = pq2;
    ipqa[k] = ipq2;
    if (nu > *nu2 + .5) {
	return 0;
    }
L120:
    pq1 = pq;
    ipq1 = ipq;
    if (nu < *nu1 + .5) {
	goto L130;
    }
    ++k;
    pqa[k] = pq;
    ipqa[k] = ipq;
    if (nu > *nu2 + .5) {
	return 0;
    }

/*        FORWARD NU-WISE RECURRENCE FOR F(MU,NU,X) FOR FIXED MU */
/*        USING */
/*        (NU+MU+1)*F(MU,NU,X)=(2.*NU+1)*F(MU,NU,X)-(NU-MU)*F(MU,NU-1,X) */
/*        WHERE F(MU,NU,X) MAY BE P(-MU,NU,X) OR IF MU IS REPLACED */
/*        BY -MU THEN F(MU,NU,X) MAY BE Q(MU,NU,X). */
/*        NOTE, IN THIS LOOP, NU=NU+1 */

L130:
    x1 = (nu * 2. - 1.) / (nu + dmu) * x * pq1;
    x2 = (nu - 1. - dmu) / (nu + dmu) * pq2;
    d__1 = -x2;
    dxadd_(&x1, &ipq1, &d__1, &ipq2, &pq, &ipq, ierror);
    if (*ierror != 0) {
	return 0;
    }
    dxadj_(&pq, &ipq, ierror);
    if (*ierror != 0) {
	return 0;
    }
    nu += 1.;
    pq2 = pq1;
    ipq2 = ipq1;
    goto L120;

} /* dxpqnu_ */

