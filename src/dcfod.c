/* dcfod.f -- translated by f2c (version 12.02.01).
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

/* DECK DCFOD */
/* Subroutine */ int dcfod_(integer *meth, doublereal *elco, doublereal *
	tesco)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, ib;
    static doublereal pc[12];
    static integer nq;
    static doublereal fnq;
    static integer nqm1, nqp1;
    static doublereal ragq, pint, xpin, fnqm1, agamq, rqfac, tsign, rq1fac;

/* ***BEGIN PROLOGUE  DCFOD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DDEBDF */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (CFOD-S, DCFOD-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   DCFOD defines coefficients needed in the integrator package DDEBDF */

/* ***SEE ALSO  DDEBDF */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820301  DATE WRITTEN */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DCFOD */


/*     ------------------------------------------------------------------ */
/*      DCFOD  IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS */
/*      NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS */
/*      GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED. */
/*      THE MAXIMUM ORDER ASSUMED HERE IS 12 IF METH = 1 AND 5 IF METH = */
/*      2.  (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.) */
/*      DCFOD  IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM, */
/*      AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED. */

/*      THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS. */
/*      THE COEFFICIENTS EL(I), 1 .LE. I .LE. NQ+1, FOR THE METHOD OF */
/*      ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A */
/*      GENERATING POLYNOMIAL, I.E., */
/*          L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ. */
/*      FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY */
/*          DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) = */
/*      0.  FOR THE BDF METHODS, L(X) IS GIVEN BY */
/*          L(X) = (X+1)*(X+2)* ... *(X+NQ)/K, */
/*      WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ). */

/*      THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE */
/*      LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER. */
/*      AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP */
/*      SIZE AT ORDER NQ - 1 IF K = 1, AT ORDER NQ IF K = 2, AND AT ORDER */
/*      NQ + 1 IF K = 3. */
/*     ------------------------------------------------------------------ */

/* ***FIRST EXECUTABLE STATEMENT  DCFOD */
    /* Parameter adjustments */
    tesco -= 4;
    elco -= 14;

    /* Function Body */
    switch (*meth) {
	case 1:  goto L10;
	case 2:  goto L60;
    }

L10:
    elco[14] = 1.;
    elco[15] = 1.;
    tesco[4] = 0.;
    tesco[5] = 2.;
    tesco[7] = 1.;
    tesco[39] = 0.;
    pc[0] = 1.;
    rqfac = 1.;
    for (nq = 2; nq <= 12; ++nq) {
/*           ------------------------------------------------------------ */
/*            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE */
/*                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ-1). */
/*            INITIALLY, P(X) = 1. */
/*           ------------------------------------------------------------ */
	rq1fac = rqfac;
	rqfac /= nq;
	nqm1 = nq - 1;
	fnqm1 = (doublereal) nqm1;
	nqp1 = nq + 1;
/*           FORM COEFFICIENTS OF P(X)*(X+NQ-1). */
/*           ---------------------------------- */
	pc[nq - 1] = 0.;
	i__1 = nqm1;
	for (ib = 1; ib <= i__1; ++ib) {
	    i__ = nqp1 - ib;
	    pc[i__ - 1] = pc[i__ - 2] + fnqm1 * pc[i__ - 1];
/* L20: */
	}
	pc[0] = fnqm1 * pc[0];
/*           COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X). */
/*           ----------------------- */
	pint = pc[0];
	xpin = pc[0] / 2.;
	tsign = 1.;
	i__1 = nq;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    tsign = -tsign;
	    pint += tsign * pc[i__ - 1] / i__;
	    xpin += tsign * pc[i__ - 1] / (i__ + 1);
/* L30: */
	}
/*           STORE COEFFICIENTS IN ELCO AND TESCO. */
/*           -------------------------------- */
	elco[nq * 13 + 1] = pint * rq1fac;
	elco[nq * 13 + 2] = 1.;
	i__1 = nq;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    elco[i__ + 1 + nq * 13] = rq1fac * pc[i__ - 1] / i__;
/* L40: */
	}
	agamq = rqfac * xpin;
	ragq = 1. / agamq;
	tesco[nq * 3 + 2] = ragq;
	if (nq < 12) {
	    tesco[nqp1 * 3 + 1] = ragq * rqfac / nqp1;
	}
	tesco[nqm1 * 3 + 3] = ragq;
/* L50: */
    }
    goto L100;

L60:
    pc[0] = 1.;
    rq1fac = 1.;
    for (nq = 1; nq <= 5; ++nq) {
/*           ------------------------------------------------------------ */
/*            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE */
/*                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ). */
/*            INITIALLY, P(X) = 1. */
/*           ------------------------------------------------------------ */
	fnq = (doublereal) nq;
	nqp1 = nq + 1;
/*           FORM COEFFICIENTS OF P(X)*(X+NQ). */
/*           ------------------------------------ */
	pc[nqp1 - 1] = 0.;
	i__1 = nq;
	for (ib = 1; ib <= i__1; ++ib) {
	    i__ = nq + 2 - ib;
	    pc[i__ - 1] = pc[i__ - 2] + fnq * pc[i__ - 1];
/* L70: */
	}
	pc[0] = fnq * pc[0];
/*           STORE COEFFICIENTS IN ELCO AND TESCO. */
/*           -------------------------------- */
	i__1 = nqp1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    elco[i__ + nq * 13] = pc[i__ - 1] / pc[1];
/* L80: */
	}
	elco[nq * 13 + 2] = 1.;
	tesco[nq * 3 + 1] = rq1fac;
	tesco[nq * 3 + 2] = nqp1 / elco[nq * 13 + 1];
	tesco[nq * 3 + 3] = (nq + 2) / elco[nq * 13 + 1];
	rq1fac /= fnq;
/* L90: */
    }
L100:
    return 0;
/*     ----------------------- END OF SUBROUTINE DCFOD */
/*     ----------------------- */
} /* dcfod_ */

