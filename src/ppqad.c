/* ppqad.f -- translated by f2c (version 12.02.01).
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

/* DECK PPQAD */
/* Subroutine */ int ppqad_(integer *ldc, real *c__, real *xi, integer *lxi, 
	integer *k, real *x1, real *x2, real *pquad)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static real a;
    static integer i__;
    static real q, s, x, aa, bb;
    static integer ii, il, im;
    static real ta, tb, dx, ss[2];
    static integer mf1, mf2, il1, il2;
    static real flk;
    static integer ilo, left;
    extern /* Subroutine */ int intrv_(real *, integer *, real *, integer *, 
	    integer *, integer *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  PPQAD */
/* ***PURPOSE  Compute the integral on (X1,X2) of a K-th order B-spline */
/*            using the piecewise polynomial (PP) representation. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A2A1, E3, K6 */
/* ***TYPE      SINGLE PRECISION (PPQAD-S, DPPQAD-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         PPQAD computes the integral on (X1,X2) of a K-th order */
/*         B-spline using the piecewise polynomial representation */
/*         (C,XI,LXI,K).  Here the Taylor expansion about the left */
/*         end point XI(J) of the J-th interval is integrated and */
/*         evaluated on subintervals of (X1,X2) which are formed by */
/*         included break points.  Integration outside (XI(1),XI(LXI+1)) */
/*         is permitted. */

/*     Description of Arguments */
/*         Input */
/*           LDC    - leading dimension of matrix C, LDC .GE. K */
/*           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI */
/*           XI(*)  - break point array of length LXI+1 */
/*           LXI    - number of polynomial pieces */
/*           K      - order of B-spline, K .GE. 1 */
/*           X1,X2  - end points of quadrature interval, normally in */
/*                    XI(1) .LE. X .LE. XI(LXI+1) */

/*         Output */
/*           PQUAD  - integral of the PP representation over (X1,X2) */

/*     Error Conditions */
/*         Improper input is a fatal error */

/* ***REFERENCES  D. E. Amos, Quadrature subroutines for splines and */
/*                 B-splines, Report SAND79-1825, Sandia Laboratories, */
/*                 December 1979. */
/* ***ROUTINES CALLED  INTRV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  PPQAD */


/* ***FIRST EXECUTABLE STATEMENT  PPQAD */
    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --xi;

    /* Function Body */
    *pquad = 0.f;
    if (*k < 1) {
	goto L100;
    }
    if (*lxi < 1) {
	goto L105;
    }
    if (*ldc < *k) {
	goto L110;
    }
    aa = dmin(*x1,*x2);
    bb = dmax(*x1,*x2);
    if (aa == bb) {
	return 0;
    }
    ilo = 1;
    intrv_(&xi[1], lxi, &aa, &ilo, &il1, &mf1);
    intrv_(&xi[1], lxi, &bb, &ilo, &il2, &mf2);
    q = 0.f;
    i__1 = il2;
    for (left = il1; left <= i__1; ++left) {
	ta = xi[left];
	a = dmax(aa,ta);
	if (left == 1) {
	    a = aa;
	}
	tb = bb;
	if (left < *lxi) {
	    tb = xi[left + 1];
	}
	x = dmin(bb,tb);
	for (ii = 1; ii <= 2; ++ii) {
	    ss[ii - 1] = 0.f;
	    dx = x - xi[left];
	    if (dx == 0.f) {
		goto L20;
	    }
	    s = c__[*k + left * c_dim1];
	    flk = (real) (*k);
	    im = *k - 1;
	    il = im;
	    i__2 = il;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		s = s * dx / flk + c__[im + left * c_dim1];
		--im;
		flk += -1.f;
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
    xermsg_("SLATEC", "PPQAD", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "PPQAD", "LXI DOES NOT SATISFY LXI.GE.1", &c__2, &c__1, 
	    (ftnlen)6, (ftnlen)5, (ftnlen)29);
    return 0;
L110:
    xermsg_("SLATEC", "PPQAD", "LDC DOES NOT SATISFY LDC.GE.K", &c__2, &c__1, 
	    (ftnlen)6, (ftnlen)5, (ftnlen)29);
    return 0;
} /* ppqad_ */

