/* bspvd.f -- translated by f2c (version 12.02.01).
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

/* DECK BSPVD */
/* Subroutine */ int bspvd_(real *t, integer *k, integer *nderiv, real *x, 
	integer *ileft, integer *ldvnik, real *vnikx, real *work)
{
    /* System generated locals */
    integer vnikx_dim1, vnikx_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, l, m;
    static real v;
    static integer jj, jm, kp1, kmd;
    static real fkmd;
    static integer jlow, mhigh, ipkmd;
    extern /* Subroutine */ int bspvn_(real *, integer *, integer *, integer *
	    , real *, integer *, real *, real *, integer *);
    static integer iwork, jp1mid;
    static real factor;
    static integer ideriv;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer ldummy;

/* ***BEGIN PROLOGUE  BSPVD */
/* ***PURPOSE  Calculate the value and all derivatives of order less than */
/*            NDERIV of all basis functions which do not vanish at X. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  E3, K6 */
/* ***TYPE      SINGLE PRECISION (BSPVD-S, DBSPVD-D) */
/* ***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     Written by Carl de Boor and modified by D. E. Amos */

/*     Abstract */
/*         BSPVD is the BSPLVD routine of the reference. */

/*         BSPVD calculates the value and all derivatives of order */
/*         less than NDERIV of all basis functions which do not */
/*         (possibly) vanish at X.  ILEFT is input such that */
/*         T(ILEFT) .LE. X .LT. T(ILEFT+1).  A call to INTRV(T,N+1,X, */
/*         ILO,ILEFT,MFLAG) will produce the proper ILEFT.  The output of */
/*         BSPVD is a matrix VNIKX(I,J) of dimension at least (K,NDERIV) */
/*         whose columns contain the K nonzero basis functions and */
/*         their NDERIV-1 right derivatives at X, I=1,K, J=1,NDERIV. */
/*         These basis functions have indices ILEFT-K+I, I=1,K, */
/*         K .LE. ILEFT .LE. N. The nonzero part of the I-th basis */
/*         function lies in (T(I),T(I+K)), I=1,N. */

/*         If X=T(ILEFT+1) then VNIKX contains left limiting values */
/*         (left derivatives) at T(ILEFT+1).  In particular, ILEFT = N */
/*         produces left limiting values at the right end point */
/*         X=T(N+1). To obtain left limiting values at T(I), I=K+1,N+1, */
/*         set X= next lower distinct knot, call INTRV to get ILEFT, */
/*         set X=T(I), and then call BSPVD. */

/*     Description of Arguments */
/*         Input */
/*          T       - knot vector of length N+K, where */
/*                    N = number of B-spline basis functions */
/*                    N = sum of knot multiplicities-K */
/*          K       - order of the B-spline, K .GE. 1 */
/*          NDERIV  - number of derivatives = NDERIV-1, */
/*                    1 .LE. NDERIV .LE. K */
/*          X       - argument of basis functions, */
/*                    T(K) .LE. X .LE. T(N+1) */
/*          ILEFT   - largest integer such that */
/*                    T(ILEFT) .LE. X .LT. T(ILEFT+1) */
/*          LDVNIK  - leading dimension of matrix VNIKX */

/*         Output */
/*          VNIKX   - matrix of dimension at least (K,NDERIV) contain- */
/*                    ing the nonzero basis functions at X and their */
/*                    derivatives columnwise. */
/*          WORK    - a work vector of length (K+1)*(K+2)/2 */

/*     Error Conditions */
/*         Improper input is a fatal error */

/* ***REFERENCES  Carl de Boor, Package for calculating with B-splines, */
/*                 SIAM Journal on Numerical Analysis 14, 3 (June 1977), */
/*                 pp. 441-472. */
/* ***ROUTINES CALLED  BSPVN, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  BSPVD */

/*     DIMENSION T(ILEFT+K), WORK((K+1)*(K+2)/2) */
/*     A(I,J) = WORK(I+J*(J+1)/2),  I=1,J+1  J=1,K-1 */
/*     A(I,K) = W0RK(I+K*(K-1)/2)  I=1.K */
/*     WORK(1) AND WORK((K+1)*(K+2)/2) ARE NOT USED. */
/* ***FIRST EXECUTABLE STATEMENT  BSPVD */
    /* Parameter adjustments */
    --t;
    vnikx_dim1 = *ldvnik;
    vnikx_offset = 1 + vnikx_dim1;
    vnikx -= vnikx_offset;
    --work;

    /* Function Body */
    if (*k < 1) {
	goto L200;
    }
    if (*nderiv < 1 || *nderiv > *k) {
	goto L205;
    }
    if (*ldvnik < *k) {
	goto L210;
    }
    ideriv = *nderiv;
    kp1 = *k + 1;
    jj = kp1 - ideriv;
    bspvn_(&t[1], &jj, k, &c__1, x, ileft, &vnikx[vnikx_offset], &work[1], &
	    iwork);
    if (ideriv == 1) {
	goto L100;
    }
    mhigh = ideriv;
    i__1 = mhigh;
    for (m = 2; m <= i__1; ++m) {
	jp1mid = 1;
	i__2 = *k;
	for (j = ideriv; j <= i__2; ++j) {
	    vnikx[j + ideriv * vnikx_dim1] = vnikx[jp1mid + vnikx_dim1];
	    ++jp1mid;
/* L10: */
	}
	--ideriv;
	jj = kp1 - ideriv;
	bspvn_(&t[1], &jj, k, &c__2, x, ileft, &vnikx[vnikx_offset], &work[1],
		 &iwork);
/* L20: */
    }

    jm = kp1 * (kp1 + 1) / 2;
    i__1 = jm;
    for (l = 1; l <= i__1; ++l) {
	work[l] = 0.f;
/* L30: */
    }
/*     A(I,I) = WORK(I*(I+3)/2) = 1.0       I = 1,K */
    l = 2;
    j = 0;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j += l;
	work[j] = 1.f;
	++l;
/* L40: */
    }
    kmd = *k;
    i__1 = mhigh;
    for (m = 2; m <= i__1; ++m) {
	--kmd;
	fkmd = (real) kmd;
	i__ = *ileft;
	j = *k;
	jj = j * (j + 1) / 2;
	jm = jj - j;
	i__2 = kmd;
	for (ldummy = 1; ldummy <= i__2; ++ldummy) {
	    ipkmd = i__ + kmd;
	    factor = fkmd / (t[ipkmd] - t[i__]);
	    i__3 = j;
	    for (l = 1; l <= i__3; ++l) {
		work[l + jj] = (work[l + jj] - work[l + jm]) * factor;
/* L50: */
	    }
	    --i__;
	    --j;
	    jj = jm;
	    jm -= j;
/* L60: */
	}

	i__2 = *k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v = 0.f;
	    jlow = max(i__,m);
	    jj = jlow * (jlow + 1) / 2;
	    i__3 = *k;
	    for (j = jlow; j <= i__3; ++j) {
		v = work[i__ + jj] * vnikx[j + m * vnikx_dim1] + v;
		jj = jj + j + 1;
/* L70: */
	    }
	    vnikx[i__ + m * vnikx_dim1] = v;
/* L80: */
	}
/* L90: */
    }
L100:
    return 0;


L200:
    xermsg_("SLATEC", "BSPVD", "K DOES NOT SATISFY K.GE.1", &c__2, &c__1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)25);
    return 0;
L205:
    xermsg_("SLATEC", "BSPVD", "NDERIV DOES NOT SATISFY 1.LE.NDERIV.LE.K", &
	    c__2, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)40);
    return 0;
L210:
    xermsg_("SLATEC", "BSPVD", "LDVNIK DOES NOT SATISFY LDVNIK.GE.K", &c__2, &
	    c__1, (ftnlen)6, (ftnlen)5, (ftnlen)35);
    return 0;
} /* bspvd_ */

