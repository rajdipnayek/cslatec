/* bsqad.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;
static integer c__2 = 2;
static integer c__1 = 1;

/* DECK BSQAD */
/* Subroutine */ int bsqad_(real *t, real *bcoef, integer *n, integer *k, 
	real *x1, real *x2, real *bquad, real *work)
{
    /* Initialized data */

    static real gpts[9] = { .577350269189625764f,.238619186083196909f,
	    .661209386466264514f,.932469514203152028f,.148874338981631211f,
	    .433395394129247191f,.679409568299024406f,.865063366688984511f,
	    .97390652851717172f };
    static real gwts[9] = { 1.f,.467913934572691047f,.360761573048138608f,
	    .171324492379170345f,.29552422471475287f,.269266719309996355f,
	    .219086362515982044f,.149451349150580593f,.0666713443086881376f };

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static real a, b;
    static integer i__, m;
    static real q, c1, y1, y2, aa, bb;
    static integer jf, mf;
    static real ta, tb, gx;
    static integer il1, il2, np1;
    static real bma, bpa;
    static integer ilo, npk;
    static real sum[5];
    static integer left, inbv, mflag;
    extern doublereal bvalu_(real *, real *, integer *, integer *, integer *, 
	    real *, integer *, real *);
    extern /* Subroutine */ int intrv_(real *, integer *, real *, integer *, 
	    integer *, integer *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BSQAD */
/* ***PURPOSE  Compute the integral of a K-th order B-spline using the */
/*            B-representation. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  H2A2A1, E3, K6 */
/* ***TYPE      SINGLE PRECISION (BSQAD-S, DBSQAD-D) */
/* ***KEYWORDS  INTEGRAL OF B-SPLINES, QUADRATURE */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*         BSQAD computes the integral on (X1,X2) of a K-th order */
/*         B-spline using the B-representation (T,BCOEF,N,K).  Orders */
/*         K as high as 20 are permitted by applying a 2, 6, or 10 */
/*         point Gauss formula on subintervals of (X1,X2) which are */
/*         formed by included (distinct) knots. */

/*         If orders K greater than 20 are needed, use BFQAD with */
/*         F(X) = 1. */

/*     Description of Arguments */
/*         Input */
/*           T      - knot array of length N+K */
/*           BCOEF  - B-spline coefficient array of length N */
/*           N      - length of coefficient array */
/*           K      - order of B-spline, 1 .LE. K .LE. 20 */
/*           X1,X2  - end points of quadrature interval in */
/*                    T(K) .LE. X .LE. T(N+1) */

/*         Output */
/*           BQUAD  - integral of the B-spline over (X1,X2) */
/*           WORK   - work vector of length 3*K */

/*     Error Conditions */
/*         Improper input is a fatal error */

/* ***REFERENCES  D. E. Amos, Quadrature subroutines for splines and */
/*                 B-splines, Report SAND79-1825, Sandia Laboratories, */
/*                 December 1979. */
/* ***ROUTINES CALLED  BVALU, INTRV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BSQAD */


    /* Parameter adjustments */
    --work;
    --bcoef;
    --t;

    /* Function Body */

/* ***FIRST EXECUTABLE STATEMENT  BSQAD */
    *bquad = 0.f;
    if (*k < 1 || *k > 20) {
	goto L65;
    }
    if (*n < *k) {
	goto L70;
    }
    aa = dmin(*x1,*x2);
    bb = dmax(*x1,*x2);
    if (aa < t[*k]) {
	goto L60;
    }
    np1 = *n + 1;
    if (bb > t[np1]) {
	goto L60;
    }
    if (aa == bb) {
	return 0;
    }
    npk = *n + *k;
/*     SELECTION OF 2, 6, OR 10 POINT GAUSS FORMULA */
    jf = 0;
    mf = 1;
    if (*k <= 4) {
	goto L10;
    }
    jf = 1;
    mf = 3;
    if (*k <= 12) {
	goto L10;
    }
    jf = 4;
    mf = 5;
L10:

    i__1 = mf;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum[i__ - 1] = 0.f;
/* L20: */
    }
    ilo = 1;
    inbv = 1;
    intrv_(&t[1], &npk, &aa, &ilo, &il1, &mflag);
    intrv_(&t[1], &npk, &bb, &ilo, &il2, &mflag);
    if (il2 >= np1) {
	il2 = *n;
    }
    i__1 = il2;
    for (left = il1; left <= i__1; ++left) {
	ta = t[left];
	tb = t[left + 1];
	if (ta == tb) {
	    goto L40;
	}
	a = dmax(aa,ta);
	b = dmin(bb,tb);
	bma = (b - a) * .5f;
	bpa = (b + a) * .5f;
	i__2 = mf;
	for (m = 1; m <= i__2; ++m) {
	    c1 = bma * gpts[jf + m - 1];
	    gx = -c1 + bpa;
	    y2 = bvalu_(&t[1], &bcoef[1], n, k, &c__0, &gx, &inbv, &work[1]);
	    gx = c1 + bpa;
	    y1 = bvalu_(&t[1], &bcoef[1], n, k, &c__0, &gx, &inbv, &work[1]);
	    sum[m - 1] += (y1 + y2) * bma;
/* L30: */
	}
L40:
	;
    }
    q = 0.f;
    i__1 = mf;
    for (m = 1; m <= i__1; ++m) {
	q += gwts[jf + m - 1] * sum[m - 1];
/* L50: */
    }
    if (*x1 > *x2) {
	q = -q;
    }
    *bquad = q;
    return 0;


L60:
    xermsg_("SLATEC", "BSQAD", "X1 OR X2 OR BOTH DO NOT SATISFY T(K).LE.X.LE"
	    ".T(N+1)", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)51);
    return 0;
L65:
    xermsg_("SLATEC", "BSQAD", "K DOES NOT SATISFY 1.LE.K.LE.20", &c__2, &
	    c__1, (ftnlen)6, (ftnlen)5, (ftnlen)31);
    return 0;
L70:
    xermsg_("SLATEC", "BSQAD", "N DOES NOT SATISFY N.GE.K", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;
} /* bsqad_ */

