/* bfqad.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BFQAD */
/* Subroutine */ int bfqad_(E_fp f, real *t, real *bcoef, integer *n, integer 
	*k, integer *id, real *x1, real *x2, real *tol, real *quad, integer *
	ierr, real *work)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real a, b, q, aa, bb, ta, tb;
    static integer il1, il2, np1;
    static real ans;
    static integer ilo, npk, iflg, left, inbv;
    static real wtol;
    extern /* Subroutine */ int bsgq8_(E_fp, real *, real *, integer *, 
	    integer *, integer *, real *, real *, integer *, real *, real *, 
	    integer *, real *);
    static integer mflag;
    extern /* Subroutine */ int intrv_(real *, integer *, real *, integer *, 
	    integer *, integer *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BFQAD */
/* ***PURPOSE  Compute the integral of a product of a function and a */
/*            derivative of a B-spline. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A2A1, E3, K6 */
/* ***TYPE      SINGLE PRECISION (BFQAD-S, DBFQAD-D) */
/* ***KEYWORDS  INTEGRAL OF B-SPLINE, QUADRATURE */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BFQAD computes the integral on (X1,X2) of a product of a */
/*         function F and the ID-th derivative of a K-th order B-spline, */
/*         using the B-representation (T,BCOEF,N,K).  (X1,X2) must be */
/*         a subinterval of T(K) .LE. X .le. T(N+1).  An integration */
/*         routine BSGQ8 (a modification */
/*         of GAUS8), integrates the product on sub- */
/*         intervals of (X1,X2) formed by included (distinct) knots. */

/*     Description of Arguments */
/*         Input */
/*           F      - external function of one argument for the */
/*                    integrand BF(X)=F(X)*BVALU(T,BCOEF,N,K,ID,X,INBV, */
/*                    WORK) */
/*           T      - knot array of length N+K */
/*           BCOEF  - coefficient array of length N */
/*           N      - length of coefficient array */
/*           K      - order of B-spline, K .GE. 1 */
/*           ID     - order of the spline derivative, 0 .LE. ID .LE. K-1 */
/*                    ID=0 gives the spline function */
/*           X1,X2  - end points of quadrature interval in */
/*                    T(K) .LE. X .LE. T(N+1) */
/*           TOL    - desired accuracy for the quadrature, suggest */
/*                    10.*STOL .LT. TOL .LE. 0.1 where STOL is the single */
/*                    precision unit roundoff for the machine = R1MACH(4) */

/*         Output */
/*           QUAD   - integral of BF(X) on (X1,X2) */
/*           IERR   - a status code */
/*                    IERR=1  normal return */
/*                         2  some quadrature on (X1,X2) does not meet */
/*                            the requested tolerance. */
/*           WORK   - work vector of length 3*K */

/*     Error Conditions */
/*         X1 or X2 not in T(K) .LE. X .LE. T(N+1) is a fatal error. */
/*         TOL not greater than the single precision unit roundoff or */
/*         less than 0.1 is a fatal error. */
/*         Some quadrature fails to meet the requested tolerance. */

/* ***REFERENCES  D. E. Amos, Quadrature subroutines for splines and */
/*                 B-splines, Report SAND79-1825, Sandia Laboratories, */
/*                 December 1979. */
/* ***ROUTINES CALLED  BSGQ8, INTRV, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BFQAD */


/* ***FIRST EXECUTABLE STATEMENT  BFQAD */
    /* Parameter adjustments */
    --work;
    --bcoef;
    --t;

    /* Function Body */
    *ierr = 1;
    *quad = 0.f;
    if (*k < 1) {
	goto L100;
    }
    if (*n < *k) {
	goto L105;
    }
    if (*id < 0 || *id >= *k) {
	goto L110;
    }
    wtol = r1mach_(&c__4);
    if (*tol < wtol || *tol > .1f) {
	goto L30;
    }
    aa = dmin(*x1,*x2);
    bb = dmax(*x1,*x2);
    if (aa < t[*k]) {
	goto L20;
    }
    np1 = *n + 1;
    if (bb > t[np1]) {
	goto L20;
    }
    if (aa == bb) {
	return 0;
    }
    npk = *n + *k;

    ilo = 1;
    intrv_(&t[1], &npk, &aa, &ilo, &il1, &mflag);
    intrv_(&t[1], &npk, &bb, &ilo, &il2, &mflag);
    if (il2 >= np1) {
	il2 = *n;
    }
    inbv = 1;
    q = 0.f;
    i__1 = il2;
    for (left = il1; left <= i__1; ++left) {
	ta = t[left];
	tb = t[left + 1];
	if (ta == tb) {
	    goto L10;
	}
	a = dmax(aa,ta);
	b = dmin(bb,tb);
	bsgq8_((E_fp)f, &t[1], &bcoef[1], n, k, id, &a, &b, &inbv, tol, &ans, 
		&iflg, &work[1]);
	if (iflg > 1) {
	    *ierr = 2;
	}
	q += ans;
L10:
	;
    }
    if (*x1 > *x2) {
	q = -q;
    }
    *quad = q;
    return 0;


L20:
    xermsg_("SLATEC", "BFQAD", "X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE"
	    ".T(N+1)", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)51);
    return 0;
L30:
    xermsg_("SLATEC", "BFQAD", "TOL IS LESS THAN THE SINGLE PRECISION TOLERA"
	    "NCE OR GREATER THAN 0.1", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (
	    ftnlen)67);
    return 0;
L100:
    xermsg_("SLATEC", "BFQAD", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "BFQAD", "N DOES NOT SATISFY N.GE.K", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;
L110:
    xermsg_("SLATEC", "BFQAD", "ID DOES NOT SATISFY 0 .LE. ID .LT. K", &c__2, 
	    &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)36);
    return 0;
} /* bfqad_ */

