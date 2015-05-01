/* bspev.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK BSPEV */
/* Subroutine */ int bspev_(real *t, real *ad, integer *n, integer *k, 
	integer *nderiv, real *x, integer *inev, real *svalue, real *work)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, l, id, jj, ll, kp1;
    static real sum;
    static integer left, kp1mn, mflag;
    extern /* Subroutine */ int bspvn_(real *, integer *, integer *, integer *
	    , real *, integer *, real *, real *, integer *);
    static integer iwork;
    extern /* Subroutine */ int intrv_(real *, integer *, real *, integer *, 
	    integer *, integer *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BSPEV */
/* ***PURPOSE  Calculate the value of the spline and its derivatives from */
/*            the B-representation. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      SINGLE PRECISION (BSPEV-S, DBSPEV-D) */
/* ***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract */
/*         BSPEV is the BSPLEV routine of the reference. */

/*         BSPEV calculates the value of the spline and its derivatives */
/*         at X from the B-representation (T,A,N,K) and returns them */
/*         in SVALUE(I),I=1,NDERIV, T(K) .LE. X .LE. T(N+1).  AD(I) can */
/*         be the B-spline coefficients A(I), I=1,N if NDERIV=1.  Other- */
/*         wise AD must be computed before hand by a call to BSPDR (T,A, */
/*         N,K,NDERIV,AD).  If X=T(I),I=K,N, right limiting values are */
/*         obtained. */

/*         To compute left derivatives or left limiting values at a */
/*         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1. */

/*         BSPEV calls INTRV, BSPVN */

/*     Description of Arguments */
/*         Input */
/*          T       - knot vector of length N+K */
/*          AD      - vector of length (2*N-NDERIV+1)*NDERIV/2 containing */
/*                    the difference table from BSPDR. */
/*          N       - number of B-spline coefficients */
/*                    N = sum of knot multiplicities-K */
/*          K       - order of the B-spline, K .GE. 1 */
/*          NDERIV  - number of derivatives, 1 .LE. NDERIV .LE. K. */
/*                    NDERIV=1 gives the zero-th derivative = function */
/*                    value */
/*          X       - argument, T(K) .LE. X .LE. T(N+1) */
/*          INEV    - an initialization parameter which must be set */
/*                    to 1 the first time BSPEV is called. */

/*         Output */
/*          INEV    - INEV contains information for efficient process- */
/*                    ing after the initial call and INEV must not */
/*                    be changed by the user.  Distinct splines require */
/*                    distinct INEV parameters. */
/*          SVALUE  - vector of length NDERIV containing the spline */
/*                    value in SVALUE(1) and the NDERIV-1 derivatives */
/*                    in the remaining components. */
/*          WORK    - work vector of length 3*K */

/*     Error Conditions */
/*         Improper input is a fatal error. */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  BSPVN, INTRV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BSPEV */

/*     DIMENSION T(N+K) */
/* ***FIRST EXECUTABLE STATEMENT  BSPEV */
    /* Parameter adjustments */
    --work;
    --svalue;
    --ad;
    --t;

    /* Function Body */
    if (*k < 1) {
	goto L100;
    }
    if (*n < *k) {
	goto L105;
    }
    if (*nderiv < 1 || *nderiv > *k) {
	goto L115;
    }
    id = *nderiv;
    i__1 = *n + 1;
    intrv_(&t[1], &i__1, x, inev, &i__, &mflag);
    if (*x < t[*k]) {
	goto L110;
    }
    if (mflag == 0) {
	goto L30;
    }
    if (*x > t[i__]) {
	goto L110;
    }
L20:
    if (i__ == *k) {
	goto L120;
    }
    --i__;
    if (*x == t[i__]) {
	goto L20;
    }

/* *I* HAS BEEN FOUND IN (K,N) SO THAT T(I) .LE. X .LT. T(I+1) */
/*     (OR .LE. T(I+1), IF T(I) .LT. T(I+1) = T(N+1) ). */
L30:
    kp1mn = *k + 1 - id;
    kp1 = *k + 1;
    bspvn_(&t[1], &kp1mn, k, &c__1, x, &i__, &work[1], &work[kp1], &iwork);
    jj = (*n + *n - id + 2) * (id - 1) / 2;
/*     ADIF(LEFTPL,ID) = AD(LEFTPL-ID+1 + (2*N-ID+2)*(ID-1)/2) */
/*     LEFTPL = LEFT + L */
L40:
    left = i__ - kp1mn;
    sum = 0.f;
    ll = left + jj + 2 - id;
    i__1 = kp1mn;
    for (l = 1; l <= i__1; ++l) {
	sum += work[l] * ad[ll];
	++ll;
/* L50: */
    }
    svalue[id] = sum;
    --id;
    if (id == 0) {
	goto L60;
    }
    jj -= *n - id + 1;
    ++kp1mn;
    bspvn_(&t[1], &kp1mn, k, &c__2, x, &i__, &work[1], &work[kp1], &iwork);
    goto L40;

L60:
    return 0;


L100:
    xermsg_("SLATEC", "BSPEV", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;
L105:
    xermsg_("SLATEC", "BSPEV", "N DOES NOT SATISFY N.GE.K", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;
L110:
    xermsg_("SLATEC", "BSPEV", "X IS NOT IN T(K).LE.X.LE.T(N+1)", &c__2, &
	    c__1, (ftnlen)6, (ftnlen)5, (ftnlen)31);
    return 0;
L115:
    xermsg_("SLATEC", "BSPEV", "NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)40);
    return 0;
L120:
    xermsg_("SLATEC", "BSPEV", "A LEFT LIMITING VALUE CANNOT BE OBTAINED AT "
	    "T(K)", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)48);
    return 0;
} /* bspev_ */

