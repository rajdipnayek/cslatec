/* polyvl.f -- translated by f2c (version 12.02.01).
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

/* DECK POLYVL */
/* Subroutine */ int polyvl_(integer *nder, real *xx, real *yfit, real *yp, 
	integer *n, real *x, real *c__, real *work, integer *ierr)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, m, mm;
    static real xk;
    static integer im1, km1;
    static real fac;
    static integer ndr;
    static real pone, ptwo;
    static integer km1pi, npkm1, nmkp1, km2pn;
    static real pione;
    static integer izero;
    static real pitwo;
    static integer km2pni;

/* ***BEGIN PROLOGUE  POLYVL */
/* ***PURPOSE  Calculate the value of a polynomial and its first NDER */
/*            derivatives where the polynomial was produced by a previous */
/*            call to POLINT. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3 */
/* ***TYPE      SINGLE PRECISION (POLYVL-S, DPOLVL-D) */
/* ***KEYWORDS  POLYNOMIAL EVALUATION */
/* ***AUTHOR  Huddleston, R. E., (SNLL) */
/* ***DESCRIPTION */

/*     Written by Robert E. Huddleston, Sandia Laboratories, Livermore */

/*     Abstract - */
/*        Subroutine POLYVL calculates the value of the polynomial and */
/*     its first NDER derivatives where the polynomial was produced by */
/*     a previous call to POLINT. */
/*        The variable N and the arrays X and C must not be altered */
/*     between the call to POLINT and the call to POLYVL. */

/*     ******  Dimensioning Information ******* */

/*     YP   must be dimensioned by at least NDER */
/*     X    must be dimensioned by at least N (see the abstract ) */
/*     C    must be dimensioned by at least N (see the abstract ) */
/*     WORK must be dimensioned by at least 2*N if NDER is .GT. 0. */

/*     *** Note *** */
/*       If NDER=0, neither YP nor WORK need to be dimensioned variables. */
/*       If NDER=1, YP does not need to be a dimensioned variable. */


/*     *****  Input parameters */

/*     NDER - the number of derivatives to be evaluated */

/*     XX   - the argument at which the polynomial and its derivatives */
/*            are to be evaluated. */

/*     N    - ***** */
/*            *       N, X, and C must not be altered between the call */
/*     X    - *       to POLINT and the call to POLYVL. */
/*     C    - ***** */


/*     *****  Output Parameters */

/*     YFIT - the value of the polynomial at XX */

/*     YP   - the derivatives of the polynomial at XX.  The derivative of */
/*            order J at XX is stored in  YP(J) , J = 1,...,NDER. */

/*     IERR - Output error flag with the following possible values. */
/*          = 1  indicates normal execution */

/*     ***** Storage Parameters */

/*     WORK  = this is an array to provide internal working storage for */
/*             POLYVL.  It must be dimensioned by at least 2*N if NDER is */
/*             .GT. 0.  If NDER=0, WORK does not need to be a dimensioned */
/*             variable. */

/* ***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston, */
/*                 Curve fitting by polynomials in one variable, Report */
/*                 SLA-74-0270, Sandia Laboratories, June 1974. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   740601  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  POLYVL */
/* ***FIRST EXECUTABLE STATEMENT  POLYVL */
    /* Parameter adjustments */
    --work;
    --c__;
    --x;
    --yp;

    /* Function Body */
    *ierr = 1;
    if (*nder > 0) {
	goto L10020;
    }

/*     *****   CODING FOR THE CASE NDER = 0 */

    pione = 1.f;
    pone = c__[1];
    *yfit = pone;
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	pitwo = (*xx - x[k - 1]) * pione;
	pione = pitwo;
	ptwo = pone + pitwo * c__[k];
	pone = ptwo;
/* L10010: */
    }
    *yfit = ptwo;
    return 0;

/*     *****   END OF NDER = 0 CASE */

L10020:
    if (*n > 1) {
	goto L10040;
    }
    *yfit = c__[1];

/*     *****  CODING FOR THE CASE  N=1 AND NDER .GT. 0 */

    i__1 = *nder;
    for (k = 1; k <= i__1; ++k) {
	yp[k] = 0.f;
/* L10030: */
    }
    return 0;

/*     *****  END OF THE CASE  N = 1 AND  NDER .GT. 0 */

L10040:
    if (*nder < *n) {
	goto L10050;
    }

/*     *****  SET FLAGS FOR NUMBER OF DERIVATIVES AND FOR DERIVATIVES */
/*            IN EXCESS OF THE DEGREE (N-1) OF THE POLYNOMIAL. */

    izero = 1;
    ndr = *n - 1;
    goto L10060;
L10050:
    izero = 0;
    ndr = *nder;
L10060:
    m = ndr + 1;
    mm = m;

/*     *****  START OF THE CASE NDER .GT. 0  AND N .GT. 1 */
/*     *****  THE POLYNOMIAL AND ITS DERIVATIVES WILL BE EVALUATED AT XX */

    i__1 = ndr;
    for (k = 1; k <= i__1; ++k) {
	yp[k] = c__[k + 1];
/* L10070: */
    }

/*     *****  THE FOLLOWING SECTION OF CODE IS EASIER TO READ IF ONE */
/*            BREAKS WORK INTO TWO ARRAYS W AND V. THE CODE WOULD THEN */
/*            READ */
/*                W(1) = 1. */
/*                PONE = C(1) */
/*               *DO   K = 2,N */
/*               *   V(K-1) =  XX - X(K-1) */
/*               *   W(K)   =  V(K-1)*W(K-1) */
/*               *   PTWO   =  PONE + W(K)*C(K) */
/*               *   PONE   =  PWO */

/*               YFIT = PTWO */

    work[1] = 1.f;
    pone = c__[1];
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	km1 = k - 1;
	npkm1 = *n + k - 1;
	work[npkm1] = *xx - x[km1];
	work[k] = work[npkm1] * work[km1];
	ptwo = pone + work[k] * c__[k];
	pone = ptwo;
/* L10080: */
    }
    *yfit = ptwo;

/*     ** AT THIS POINT THE POLYNOMIAL HAS BEEN EVALUATED AND INFORMATION */
/*        FOR THE DERIVATIVE EVALUATIONS HAVE BEEN STORED IN THE ARRAY */
/*        WORK */
    if (*n == 2) {
	goto L10110;
    }
    if (m == *n) {
	mm = ndr;
    }

/*     ***** EVALUATE THE DERIVATIVES AT XX */

/*                  ******  DO K=2,MM   (FOR MOST CASES, MM = NDER + 1) */
/*                  *  ******  DO I=2,N-K+1 */
/*                  *  *       W(I) = V(K-2+I)*W(I-1) + W(I) */
/*                  *  *       YP(K-1) = YP(K-1) + W(I)*C(K-1+I) */
/*                  ******  CONTINUE */

    i__1 = mm;
    for (k = 2; k <= i__1; ++k) {
	nmkp1 = *n - k + 1;
	km1 = k - 1;
	km2pn = k - 2 + *n;
	i__2 = nmkp1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    km2pni = km2pn + i__;
	    im1 = i__ - 1;
	    km1pi = km1 + i__;
	    work[i__] = work[km2pni] * work[im1] + work[i__];
	    yp[km1] += work[i__] * c__[km1pi];
/* L10090: */
	}
    }
    if (ndr == 1) {
	goto L10110;
    }
    fac = 1.f;
    i__2 = ndr;
    for (k = 2; k <= i__2; ++k) {
	xk = (real) k;
	fac = xk * fac;
	yp[k] = fac * yp[k];
/* L10100: */
    }

/*     ***** END OF DERIVATIVE EVALUATIONS */

L10110:
    if (izero == 0) {
	return 0;
    }

/*     *****  SET EXCESS DERIVATIVES TO ZERO. */

    i__2 = *nder;
    for (k = *n; k <= i__2; ++k) {
	yp[k] = 0.f;
/* L10120: */
    }
    return 0;
} /* polyvl_ */

