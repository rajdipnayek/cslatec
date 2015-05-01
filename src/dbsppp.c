/* dbsppp.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static integer c__1 = 1;

/* DECK DBSPPP */
/* Subroutine */ int dbsppp_(doublereal *t, doublereal *a, integer *n, 
	integer *k, integer *ldc, doublereal *c__, doublereal *xi, integer *
	lxi, doublereal *work)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;

    /* Local variables */
    static integer nk, inev, ileft;
    extern /* Subroutine */ int dbspdr_(doublereal *, doublereal *, integer *,
	     integer *, integer *, doublereal *), dbspev_(doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *), xermsg_(char *, char *, 
	    char *, integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBSPPP */
/* ***PURPOSE  Convert the B-representation of a B-spline to the piecewise */
/*            polynomial (PP) form. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      DOUBLE PRECISION (BSPPP-S, DBSPPP-D) */
/* ***KEYWORDS  B-SPLINE, PIECEWISE POLYNOMIAL */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract    **** a double precision routine **** */
/*         DBSPPP is the BSPLPP routine of the reference. */

/*         DBSPPP converts the B-representation (T,A,N,K) to the */
/*         piecewise polynomial (PP) form (C,XI,LXI,K) for use with */
/*         DPPVAL.  Here XI(*), the break point array of length LXI, is */
/*         the knot array T(*) with multiplicities removed.  The columns */
/*         of the matrix C(I,J) contain the right Taylor derivatives */
/*         for the polynomial expansion about XI(J) for the intervals */
/*         XI(J) .LE. X .LE. XI(J+1), I=1,K, J=1,LXI.  Function DPPVAL */
/*         makes this evaluation at a specified point X in */
/*         XI(1) .LE. X .LE. XI(LXI+1) */

/*     Description of Arguments */

/*         Input      T,A are double precision */
/*          T       - knot vector of length N+K */
/*          A       - B-spline coefficient vector of length N */
/*          N       - number of B-spline coefficients */
/*                    N = sum of knot multiplicities-K */
/*          K       - order of the B-spline, K .GE. 1 */
/*          LDC     - leading dimension of C, LDC .GE. K */

/*         Output     C,XI,WORK are double precision */
/*          C       - matrix of dimension at least (K,LXI) containing */
/*                    right derivatives at break points */
/*          XI      - XI break point vector of length LXI+1 */
/*          LXI     - number of break points, LXI .LE. N-K+1 */
/*          WORK    - work vector of length K*(N+3) */

/*     Error Conditions */
/*         Improper input is a fatal error */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  DBSPDR, DBSPEV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DBSPPP */

/*     DIMENSION T(N+K),XI(LXI+1),C(LDC,*) */
/*     HERE, * = THE FINAL VALUE OF THE OUTPUT PARAMETER LXI. */
/* ***FIRST EXECUTABLE STATEMENT  DBSPPP */
    /* Parameter adjustments */
    --t;
    --a;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --xi;
    --work;

    /* Function Body */
    if (*k < 1) {
	goto L100;
    }
    if (*n < *k) {
	goto L105;
    }
    if (*ldc < *k) {
	goto L110;
    }
    dbspdr_(&t[1], &a[1], n, k, k, &work[1]);
    *lxi = 0;
    xi[1] = t[*k];
    inev = 1;
    nk = *n * *k + 1;
    i__1 = *n;
    for (ileft = *k; ileft <= i__1; ++ileft) {
	if (t[ileft + 1] == t[ileft]) {
	    goto L10;
	}
	++(*lxi);
	xi[*lxi + 1] = t[ileft + 1];
	dbspev_(&t[1], &work[1], n, k, k, &xi[*lxi], &inev, &c__[*lxi * 
		c_dim1 + 1], &work[nk]);
L10:
	;
    }
    return 0;
L100:
    xermsg_("SLATEC", "DBSPPP", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "DBSPPP", "N DOES NOT SATISFY N.GE.K", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L110:
    xermsg_("SLATEC", "DBSPPP", "LDC DOES NOT SATISFY LDC.GE.K", &c__2, &c__1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)29);
    return 0;
} /* dbsppp_ */

