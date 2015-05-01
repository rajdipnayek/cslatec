/* bvalu.f -- translated by f2c (version 12.02.01).
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

/* DECK BVALU */
doublereal bvalu_(real *t, real *a, integer *n, integer *k, integer *ideriv, 
	real *x, integer *inbv, real *work)
{
    /* System generated locals */
    integer i__1, i__2;
    real ret_val;

    /* Local variables */
    static integer i__, j, j1, j2, jj, km1, ip1, ihi, imk, kmj, ipj, ilo, kpk;
    static real fkmj;
    static integer ip1mj, mflag, imkpj;
    extern /* Subroutine */ int intrv_(real *, integer *, real *, integer *, 
	    integer *, integer *);
    static integer iderp1, kmider, ihmkmj;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  BVALU */
/* ***PURPOSE  Evaluate the B-representation of a B-spline at X for the */
/*            function value or any of its derivatives. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      SINGLE PRECISION (BVALU-S, DBVALU-D) */
/* ***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract */
/*         BVALU is the BVALUE function of the reference. */

/*         BVALU evaluates the B-representation (T,A,N,K) of a B-spline */
/*         at X for the function value on IDERIV = 0 or any of its */
/*         derivatives on IDERIV = 1,2,...,K-1.  Right limiting values */
/*         (right derivatives) are returned except at the right end */
/*         point X=T(N+1) where left limiting values are computed.  The */
/*         spline is defined on T(K) .LE. X .LE. T(N+1).  BVALU returns */
/*         a fatal error message when X is outside of this interval. */

/*         To compute left derivatives or left limiting values at a */
/*         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1. */

/*         BVALU calls INTRV */

/*     Description of Arguments */
/*         Input */
/*          T       - knot vector of length N+K */
/*          A       - B-spline coefficient vector of length N */
/*          N       - number of B-spline coefficients */
/*                    N = sum of knot multiplicities-K */
/*          K       - order of the B-spline, K .GE. 1 */
/*          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1 */
/*                    IDERIV=0 returns the B-spline value */
/*          X       - argument, T(K) .LE. X .LE. T(N+1) */
/*          INBV    - an initialization parameter which must be set */
/*                    to 1 the first time BVALU is called. */

/*         Output */
/*          INBV    - INBV contains information for efficient process- */
/*                    ing after the initial call and INBV must not */
/*                    be changed by the user.  Distinct splines require */
/*                    distinct INBV parameters. */
/*          WORK    - work vector of length 3*K. */
/*          BVALU   - value of the IDERIV-th derivative at X */

/*     Error Conditions */
/*         An improper input is a fatal error */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  INTRV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BVALU */

/*     DIMENSION T(N+K), WORK(3*K) */
/* ***FIRST EXECUTABLE STATEMENT  BVALU */
    /* Parameter adjustments */
    --work;
    --a;
    --t;

    /* Function Body */
    ret_val = 0.f;
    if (*k < 1) {
	goto L102;
    }
    if (*n < *k) {
	goto L101;
    }
    if (*ideriv < 0 || *ideriv >= *k) {
	goto L110;
    }
    kmider = *k - *ideriv;

/* *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1) */
/*     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)). */
    km1 = *k - 1;
    i__1 = *n + 1;
    intrv_(&t[1], &i__1, x, inbv, &i__, &mflag);
    if (*x < t[*k]) {
	goto L120;
    }
    if (mflag == 0) {
	goto L20;
    }
    if (*x > t[i__]) {
	goto L130;
    }
L10:
    if (i__ == *k) {
	goto L140;
    }
    --i__;
    if (*x == t[i__]) {
	goto L10;
    }

/* *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES */
/*     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K */

L20:
    imk = i__ - *k;
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	imkpj = imk + j;
	work[j] = a[imkpj];
/* L30: */
    }
    if (*ideriv == 0) {
	goto L60;
    }
    i__1 = *ideriv;
    for (j = 1; j <= i__1; ++j) {
	kmj = *k - j;
	fkmj = (real) kmj;
	i__2 = kmj;
	for (jj = 1; jj <= i__2; ++jj) {
	    ihi = i__ + jj;
	    ihmkmj = ihi - kmj;
	    work[jj] = (work[jj + 1] - work[jj]) / (t[ihi] - t[ihmkmj]) * 
		    fkmj;
/* L40: */
	}
/* L50: */
    }

/* *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE, */
/*     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV). */
L60:
    if (*ideriv == km1) {
	goto L100;
    }
    ip1 = i__ + 1;
    kpk = *k + *k;
    j1 = *k + 1;
    j2 = kpk + 1;
    i__1 = kmider;
    for (j = 1; j <= i__1; ++j) {
	ipj = i__ + j;
	work[j1] = t[ipj] - *x;
	ip1mj = ip1 - j;
	work[j2] = *x - t[ip1mj];
	++j1;
	++j2;
/* L70: */
    }
    iderp1 = *ideriv + 1;
    i__1 = km1;
    for (j = iderp1; j <= i__1; ++j) {
	kmj = *k - j;
	ilo = kmj;
	i__2 = kmj;
	for (jj = 1; jj <= i__2; ++jj) {
	    work[jj] = (work[jj + 1] * work[kpk + ilo] + work[jj] * work[*k + 
		    jj]) / (work[kpk + ilo] + work[*k + jj]);
	    --ilo;
/* L80: */
	}
/* L90: */
    }
L100:
    ret_val = work[1];
    return ret_val;


L101:
    xermsg_("SLATEC", "BVALU", "N DOES NOT SATISFY N.GE.K", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return ret_val;
L102:
    xermsg_("SLATEC", "BVALU", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return ret_val;
L110:
    xermsg_("SLATEC", "BVALU", "IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)40);
    return ret_val;
L120:
    xermsg_("SLATEC", "BVALU", "X IS N0T GREATER THAN OR EQUAL TO T(K)", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)38);
    return ret_val;
L130:
    xermsg_("SLATEC", "BVALU", "X IS NOT LESS THAN OR EQUAL TO T(N+1)", &c__2,
	     &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)37);
    return ret_val;
L140:
    xermsg_("SLATEC", "BVALU", "A LEFT LIMITING VALUE CANNOT BE OBTAINED AT "
	    "T(K)", &c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)48);
    return ret_val;
} /* bvalu_ */

