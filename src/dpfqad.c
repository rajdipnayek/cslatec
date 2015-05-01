/* dpfqad.f -- translated by f2c (version 12.02.01).
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

/* DECK DPFQAD */
/* Subroutine */ int dpfqad_(D_fp f, integer *ldc, doublereal *c__, 
	doublereal *xi, integer *lxi, integer *k, integer *id, doublereal *x1,
	 doublereal *x2, doublereal *tol, doublereal *quad, integer *ierr)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;

    /* Local variables */
    static doublereal a, b, q, aa, bb, ta, tb;
    static integer mf1, mf2, il1, il2;
    static doublereal ans;
    static integer ilo, iflg, left;
    static doublereal wtol;
    static integer inppv;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int dppgq8_(D_fp, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    xermsg_(char *, char *, char *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen), dintrv_(doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *);

/* ***BEGIN PROLOGUE  DPFQAD */
/* ***PURPOSE  Compute the integral on (X1,X2) of a product of a */
/*            function F and the ID-th derivative of a B-spline, */
/*            (PP-representation). */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A2A1, E3, K6 */
/* ***TYPE      DOUBLE PRECISION (PFQAD-S, DPFQAD-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract    **** a double precision routine **** */
/*         DPFQAD computes the integral on (X1,X2) of a product of a */
/*         function F and the ID-th derivative of a B-spline, using the */
/*         PP-representation (C,XI,LXI,K).  (X1,X2) is normally a sub- */
/*         interval of XI(1) .LE. X .LE. XI(LXI+1).  An integration */
/*         routine, DPPGQ8 (a modification of GAUS8), integrates the */
/*         product on subintervals of (X1,X2) formed by the included */
/*         break points.  Integration outside of (XI(1),XI(LXI+1)) is */
/*         permitted provided F is defined. */

/*         The maximum number of significant digits obtainable in */
/*         DBSQAD is the smaller of 18 and the number of digits */
/*         carried in double precision arithmetic. */

/*     Description of arguments */
/*         Input      F,C,XI,X1,X2,TOL are double precision */
/*           F      - external function of one argument for the */
/*                    integrand PF(X)=F(X)*DPPVAL(LDC,C,XI,LXI,K,ID,X, */
/*                    INPPV) */
/*           LDC    - leading dimension of matrix C, LDC .GE. K */
/*           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI */
/*           XI(*)  - break point array of length LXI+1 */
/*           LXI    - number of polynomial pieces */
/*           K      - order of B-spline, K .GE. 1 */
/*           ID     - order of the spline derivative, 0 .LE. ID .LE. K-1 */
/*                    ID=0 gives the spline function */
/*           X1,X2  - end points of quadrature interval, normally in */
/*                    XI(1) .LE. X .LE. XI(LXI+1) */
/*           TOL    - desired accuracy for the quadrature, suggest */
/*                    10.*DTOL .LT. TOL .LE. 0.1 where DTOL is the */
/*                    maximum of 1.0D-18 and double precision unit */
/*                    roundoff for the machine = D1MACH(4) */

/*         Output     QUAD is double precision */
/*           QUAD   - integral of PF(X) on (X1,X2) */
/*           IERR   - a status code */
/*                    IERR=1 normal return */
/*                         2 some quadrature does not meet the */
/*                           requested tolerance */

/*     Error Conditions */
/*         Improper input is a fatal error. */
/*         Some quadrature does not meet the requested tolerance. */

/* ***REFERENCES  D. E. Amos, Quadrature subroutines for splines and */
/*                 B-splines, Report SAND79-1825, Sandia Laboratories, */
/*                 December 1979. */
/* ***ROUTINES CALLED  D1MACH, DINTRV, DPPGQ8, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPFQAD */


/* ***FIRST EXECUTABLE STATEMENT  DPFQAD */
    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --xi;

    /* Function Body */
    *ierr = 1;
    *quad = 0.;
    if (*k < 1) {
	goto L100;
    }
    if (*ldc < *k) {
	goto L105;
    }
    if (*id < 0 || *id >= *k) {
	goto L110;
    }
    if (*lxi < 1) {
	goto L115;
    }
    wtol = d1mach_(&c__4);
    wtol = max(wtol,1e-18);
    if (*tol < wtol || *tol > .1) {
	goto L20;
    }
    aa = min(*x1,*x2);
    bb = max(*x1,*x2);
    if (aa == bb) {
	return 0;
    }
    ilo = 1;
    dintrv_(&xi[1], lxi, &aa, &ilo, &il1, &mf1);
    dintrv_(&xi[1], lxi, &bb, &ilo, &il2, &mf2);
    q = 0.;
    inppv = 1;
    i__1 = il2;
    for (left = il1; left <= i__1; ++left) {
	ta = xi[left];
	a = max(aa,ta);
	if (left == 1) {
	    a = aa;
	}
	tb = bb;
	if (left < *lxi) {
	    tb = xi[left + 1];
	}
	b = min(bb,tb);
	dppgq8_((D_fp)f, ldc, &c__[c_offset], &xi[1], lxi, k, id, &a, &b, &
		inppv, tol, &ans, &iflg);
	if (iflg > 1) {
	    *ierr = 2;
	}
	q += ans;
/* L10: */
    }
    if (*x1 > *x2) {
	q = -q;
    }
    *quad = q;
    return 0;

L20:
    xermsg_("SLATEC", "DPFQAD", "TOL IS LESS DTOL OR GREATER THAN 0.1", &c__2,
	     &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)36);
    return 0;
L100:
    xermsg_("SLATEC", "DPFQAD", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "DPFQAD", "LDC DOES NOT SATISFY LDC.GE.K", &c__2, &c__1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)29);
    return 0;
L110:
    xermsg_("SLATEC", "DPFQAD", "ID DOES NOT SATISFY 0.LE.ID.LT.K", &c__2, &
	    c__1, (ftnlen)6, (ftnlen)6, (ftnlen)32);
    return 0;
L115:
    xermsg_("SLATEC", "DPFQAD", "LXI DOES NOT SATISFY LXI.GE.1", &c__2, &c__1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)29);
    return 0;
} /* dpfqad_ */

