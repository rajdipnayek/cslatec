/* dppqad.f -- translated by f2c (version 12.02.01).
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

/* DECK DPPQAD */
/* Subroutine */ int dppqad_(integer *ldc, doublereal *c__, doublereal *xi, 
	integer *lxi, integer *k, doublereal *x1, doublereal *x2, doublereal *
	pquad)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static doublereal a;
    static integer i__;
    static doublereal q, s, x, aa, bb;
    static integer ii, il, im;
    static doublereal ta, tb, dx, ss[2];
    static integer mf1, mf2, il1, il2;
    static doublereal flk;
    static integer ilo, left;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dintrv_(doublereal *, integer 
	    *, doublereal *, integer *, integer *, integer *);

/* ***BEGIN PROLOGUE  DPPQAD */
/* ***PURPOSE  Compute the integral on (X1,X2) of a K-th order B-spline */
/*            using the piecewise polynomial (PP) representation. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A2A1, E3, K6 */
/* ***TYPE      DOUBLE PRECISION (PPQAD-S, DPPQAD-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract    **** a double precision routine **** */
/*         DPPQAD computes the integral on (X1,X2) of a K-th order */
/*         B-spline using the piecewise polynomial representation */
/*         (C,XI,LXI,K).  Here the Taylor expansion about the left */
/*         end point XI(J) of the J-th interval is integrated and */
/*         evaluated on subintervals of (X1,X2) which are formed by */
/*         included break points.  Integration outside (XI(1),XI(LXI+1)) */
/*         is permitted. */

/*     Description of Arguments */
/*         Input      C,XI,X1,X2 are double precision */
/*           LDC    - leading dimension of matrix C, LDC .GE. K */
/*           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI */
/*           XI(*)  - break point array of length LXI+1 */
/*           LXI    - number of polynomial pieces */
/*           K      - order of B-spline, K .GE. 1 */
/*           X1,X2  - end points of quadrature interval, normally in */
/*                    XI(1) .LE. X .LE. XI(LXI+1) */

/*         Output     PQUAD is double precision */
/*           PQUAD  - integral of the PP representation over (X1,X2) */

/*     Error Conditions */
/*         Improper input is a fatal error */

/* ***REFERENCES  D. E. Amos, Quadrature subroutines for splines and */
/*                 B-splines, Report SAND79-1825, Sandia Laboratories, */
/*                 December 1979. */
/* ***ROUTINES CALLED  DINTRV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   890911  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPPQAD */


/* ***FIRST EXECUTABLE STATEMENT  DPPQAD */
    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --xi;

    /* Function Body */
    *pquad = 0.;
    if (*k < 1) {
	goto L100;
    }
    if (*lxi < 1) {
	goto L105;
    }
    if (*ldc < *k) {
	goto L110;
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
	x = min(bb,tb);
	for (ii = 1; ii <= 2; ++ii) {
	    ss[ii - 1] = 0.;
	    dx = x - xi[left];
	    if (dx == 0.) {
		goto L20;
	    }
	    s = c__[*k + left * c_dim1];
	    flk = (doublereal) (*k);
	    im = *k - 1;
	    il = im;
	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s = s * dx / flk + c__[im + left * c_dim1];
		--im;
		flk += -1.;
/* L10: */
	    }
	    ss[ii - 1] = s * dx;
L20:
	    x = a;
/* L30: */
	}
	q += ss[0] - ss[1];
/* L40: */
    }
    if (*x1 > *x2) {
	q = -q;
    }
    *pquad = q;
    return 0;


L100:
    xermsg_("SLATEC", "DPPQAD", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)6, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "DPPQAD", "LXI DOES NOT SATISFY LXI.GE.1", &c__2, &c__1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)29);
    return 0;
L110:
    xermsg_("SLATEC", "DPPQAD", "LDC DOES NOT SATISFY LDC.GE.K", &c__2, &c__1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)29);
    return 0;
} /* dppqad_ */

