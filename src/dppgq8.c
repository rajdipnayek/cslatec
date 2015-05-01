/* dppgq8.f -- translated by f2c (version 12.02.01).
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

static integer c__14 = 14;
static integer c__5 = 5;
static doublereal c_b6 = 1.;
static doublereal c_b8 = 2.;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c_n1 = -1;

/* DECK DPPGQ8 */
/* Subroutine */ int dppgq8_(D_fp fun, integer *ldc, doublereal *c__, 
	doublereal *xi, integer *lxi, integer *kk, integer *id, doublereal *a,
	 doublereal *b, integer *inppv, doublereal *err, doublereal *ans, 
	integer *ierr)
{
    /* Initialized data */

    static doublereal x1 = .183434642495649805;
    static integer nlmn = 1;
    static integer kmx = 5000;
    static integer kml = 6;
    static doublereal x2 = .525532409916328986;
    static doublereal x3 = .79666647741362674;
    static doublereal x4 = .960289856497536232;
    static doublereal w1 = .362683783378361983;
    static doublereal w2 = .313706645877887287;
    static doublereal w3 = .222381034453374471;
    static doublereal w4 = .101228536290376259;
    static doublereal sq2 = 1.41421356;

    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14, d__15, d__16, d__17, d__18;

    /* Local variables */
    static integer k, l;
    static doublereal aa[60], ae, be, cc, ee, ef, hh[60], gl, gr[60];
    static integer lr[60];
    static doublereal vl[60], vr;
    static integer nib;
    static doublereal glr;
    static integer lmn;
    static doublereal eps, est, tol;
    static integer lmx, mxl;
    static doublereal area, anib;
    static integer nlmx, nbits;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);
    extern doublereal dppval_(integer *, doublereal *, doublereal *, integer *
	    , integer *, integer *, doublereal *, integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DPPGQ8 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DPFQAD */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (PPGQ8-S, DPPGQ8-D) */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract    **** A DOUBLE PRECISION routine **** */

/*        DPPGQ8, a modification of GAUS8, integrates the */
/*        product of FUN(X) by the ID-th derivative of a spline */
/*        DPPVAL(LDC,C,XI,LXI,KK,ID,X,INPPV)  between limits A and B. */

/*     Description of Arguments */

/*      Input-- FUN,C,XI,A,B,ERR are DOUBLE PRECISION */
/*        FUN - Name of external function of one argument which */
/*              multiplies DPPVAL. */
/*        LDC - Leading dimension of matrix C, LDC .GE. KK */
/*        C   - Matrix of Taylor derivatives of dimension at least */
/*              (K,LXI) */
/*        XI  - Breakpoint vector of length LXI+1 */
/*        LXI - Number of polynomial pieces */
/*        KK  - Order of the spline, KK .GE. 1 */
/*        ID  - Order of the spline derivative, 0 .LE. ID .LE. KK-1 */
/*        A   - Lower limit of integral */
/*        B   - Upper limit of integral (may be less than A) */
/*        INPPV- Initialization parameter for DPPVAL */
/*        ERR - Is a requested pseudorelative error tolerance.  Normally */
/*              pick a value of ABS(ERR) .LT. 1D-3.  ANS will normally */
/*              have no more error than ABS(ERR) times the integral of */
/*              the absolute value of FUN(X)*DPPVAL(LDC,C,XI,LXI,KK,ID,X, */
/*              INPPV). */


/*      Output-- ERR,ANS are DOUBLE PRECISION */
/*        ERR - Will be an estimate of the absolute error in ANS if the */
/*              input value of ERR was negative.  (ERR Is unchanged if */
/*              the input value of ERR was nonnegative.)  The estimated */
/*              error is solely for information to the user and should */
/*              not be used as a correction to the computed integral. */
/*        ANS - Computed value of integral */
/*        IERR- A status code */
/*            --Normal Codes */
/*               1 ANS most likely meets requested error tolerance, */
/*                 or A=B. */
/*              -1 A and B are too nearly equal to allow normal */
/*                 integration.  ANS is set to zero. */
/*            --Abnormal Code */
/*               2 ANS probably does not meet requested error tolerance. */

/* ***SEE ALSO  DPFQAD */
/* ***ROUTINES CALLED  D1MACH, DPPVAL, I1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890911  Removed unnecessary intrinsics.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  DPPGQ8 */

    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --xi;

    /* Function Body */

/*     INITIALIZE */

/* ***FIRST EXECUTABLE STATEMENT  DPPGQ8 */
    k = i1mach_(&c__14);
    anib = d1mach_(&c__5) * k / .30102;
    nbits = (integer) anib;
/* Computing MIN */
    i__1 = nbits * 5 / 8;
    nlmx = min(i__1,60);
    *ans = 0.;
    *ierr = 1;
    be = 0.;
    if (*a == *b) {
	goto L140;
    }
    lmx = nlmx;
    lmn = nlmn;
    if (*b == 0.) {
	goto L10;
    }
    if (d_sign(&c_b6, b) * *a <= 0.) {
	goto L10;
    }
    cc = (d__1 = 1. - *a / *b, abs(d__1));
    if (cc > .1) {
	goto L10;
    }
    if (cc <= 0.) {
	goto L140;
    }
    anib = .5 - log(cc) / .69314718;
    nib = (integer) anib;
/* Computing MIN */
    i__1 = nlmx, i__2 = nbits - nib - 7;
    lmx = min(i__1,i__2);
    if (lmx < 1) {
	goto L130;
    }
    lmn = min(lmn,lmx);
L10:
/* Computing MAX */
    i__1 = 5 - nbits;
    d__1 = abs(*err), d__2 = pow_di(&c_b8, &i__1);
    tol = max(d__1,d__2) / 2.;
    if (*err == 0.) {
	tol = sqrt(d1mach_(&c__4));
    }
    eps = tol;
    hh[0] = (*b - *a) / 4.;
    aa[0] = *a;
    lr[0] = 1;
    l = 1;
    d__1 = aa[l - 1] + hh[l - 1] * 2.;
    d__2 = hh[l - 1] * 2.;
    d__3 = d__1 - x1 * d__2;
    d__4 = d__1 - x1 * d__2;
    d__5 = d__1 + x1 * d__2;
    d__6 = d__1 + x1 * d__2;
    d__7 = d__1 - x2 * d__2;
    d__8 = d__1 - x2 * d__2;
    d__9 = d__1 + x2 * d__2;
    d__10 = d__1 + x2 * d__2;
    d__11 = d__1 - x3 * d__2;
    d__12 = d__1 - x3 * d__2;
    d__13 = d__1 + x3 * d__2;
    d__14 = d__1 + x3 * d__2;
    d__15 = d__1 - x4 * d__2;
    d__16 = d__1 - x4 * d__2;
    d__17 = d__1 + x4 * d__2;
    d__18 = d__1 + x4 * d__2;
    est = d__2 * (w1 * ((*fun)(&d__3) * dppval_(ldc, &c__[c_offset], &xi[1], 
	    lxi, kk, id, &d__4, inppv) + (*fun)(&d__5) * dppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &d__6, inppv)) + w2 * ((*fun)(&
	    d__7) * dppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &d__8, 
	    inppv) + (*fun)(&d__9) * dppval_(ldc, &c__[c_offset], &xi[1], lxi,
	     kk, id, &d__10, inppv)) + (w3 * ((*fun)(&d__11) * dppval_(ldc, &
	    c__[c_offset], &xi[1], lxi, kk, id, &d__12, inppv) + (*fun)(&
	    d__13) * dppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &d__14,
	     inppv)) + w4 * ((*fun)(&d__15) * dppval_(ldc, &c__[c_offset], &
	    xi[1], lxi, kk, id, &d__16, inppv) + (*fun)(&d__17) * dppval_(ldc,
	     &c__[c_offset], &xi[1], lxi, kk, id, &d__18, inppv))));
    k = 8;
    area = abs(est);
    ef = .5;
    mxl = 0;

/*     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC. */

L20:
    d__1 = aa[l - 1] + hh[l - 1];
    d__2 = d__1 - x1 * hh[l - 1];
    d__3 = d__1 - x1 * hh[l - 1];
    d__4 = d__1 + x1 * hh[l - 1];
    d__5 = d__1 + x1 * hh[l - 1];
    d__6 = d__1 - x2 * hh[l - 1];
    d__7 = d__1 - x2 * hh[l - 1];
    d__8 = d__1 + x2 * hh[l - 1];
    d__9 = d__1 + x2 * hh[l - 1];
    d__10 = d__1 - x3 * hh[l - 1];
    d__11 = d__1 - x3 * hh[l - 1];
    d__12 = d__1 + x3 * hh[l - 1];
    d__13 = d__1 + x3 * hh[l - 1];
    d__14 = d__1 - x4 * hh[l - 1];
    d__15 = d__1 - x4 * hh[l - 1];
    d__16 = d__1 + x4 * hh[l - 1];
    d__17 = d__1 + x4 * hh[l - 1];
    gl = hh[l - 1] * (w1 * ((*fun)(&d__2) * dppval_(ldc, &c__[c_offset], &xi[
	    1], lxi, kk, id, &d__3, inppv) + (*fun)(&d__4) * dppval_(ldc, &
	    c__[c_offset], &xi[1], lxi, kk, id, &d__5, inppv)) + w2 * ((*fun)(
	    &d__6) * dppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &d__7, 
	    inppv) + (*fun)(&d__8) * dppval_(ldc, &c__[c_offset], &xi[1], lxi,
	     kk, id, &d__9, inppv)) + (w3 * ((*fun)(&d__10) * dppval_(ldc, &
	    c__[c_offset], &xi[1], lxi, kk, id, &d__11, inppv) + (*fun)(&
	    d__12) * dppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &d__13,
	     inppv)) + w4 * ((*fun)(&d__14) * dppval_(ldc, &c__[c_offset], &
	    xi[1], lxi, kk, id, &d__15, inppv) + (*fun)(&d__16) * dppval_(ldc,
	     &c__[c_offset], &xi[1], lxi, kk, id, &d__17, inppv))));
    d__1 = aa[l - 1] + hh[l - 1] * 3.;
    d__2 = d__1 - x1 * hh[l - 1];
    d__3 = d__1 - x1 * hh[l - 1];
    d__4 = d__1 + x1 * hh[l - 1];
    d__5 = d__1 + x1 * hh[l - 1];
    d__6 = d__1 - x2 * hh[l - 1];
    d__7 = d__1 - x2 * hh[l - 1];
    d__8 = d__1 + x2 * hh[l - 1];
    d__9 = d__1 + x2 * hh[l - 1];
    d__10 = d__1 - x3 * hh[l - 1];
    d__11 = d__1 - x3 * hh[l - 1];
    d__12 = d__1 + x3 * hh[l - 1];
    d__13 = d__1 + x3 * hh[l - 1];
    d__14 = d__1 - x4 * hh[l - 1];
    d__15 = d__1 - x4 * hh[l - 1];
    d__16 = d__1 + x4 * hh[l - 1];
    d__17 = d__1 + x4 * hh[l - 1];
    gr[l - 1] = hh[l - 1] * (w1 * ((*fun)(&d__2) * dppval_(ldc, &c__[c_offset]
	    , &xi[1], lxi, kk, id, &d__3, inppv) + (*fun)(&d__4) * dppval_(
	    ldc, &c__[c_offset], &xi[1], lxi, kk, id, &d__5, inppv)) + w2 * ((
	    *fun)(&d__6) * dppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &
	    d__7, inppv) + (*fun)(&d__8) * dppval_(ldc, &c__[c_offset], &xi[1]
	    , lxi, kk, id, &d__9, inppv)) + (w3 * ((*fun)(&d__10) * dppval_(
	    ldc, &c__[c_offset], &xi[1], lxi, kk, id, &d__11, inppv) + (*fun)(
	    &d__12) * dppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &
	    d__13, inppv)) + w4 * ((*fun)(&d__14) * dppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &d__15, inppv) + (*fun)(&d__16) * 
	    dppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &d__17, inppv)))
	    );
    k += 16;
    area += abs(gl) + (d__1 = gr[l - 1], abs(d__1)) - abs(est);
    glr = gl + gr[l - 1];
    ee = (d__1 = est - glr, abs(d__1)) * ef;
/* Computing MAX */
    d__1 = eps * area, d__2 = tol * abs(glr);
    ae = max(d__1,d__2);
    if (ee - ae <= 0.) {
	goto L40;
    } else {
	goto L50;
    }
L30:
    mxl = 1;
L40:
    be += est - glr;
    if (lr[l - 1] <= 0) {
	goto L60;
    } else {
	goto L80;
    }

/*     CONSIDER THE LEFT HALF OF THIS LEVEL */

L50:
    if (k > kmx) {
	lmx = kml;
    }
    if (l >= lmx) {
	goto L30;
    }
    ++l;
    eps *= .5;
    ef /= sq2;
    hh[l - 1] = hh[l - 2] * .5;
    lr[l - 1] = -1;
    aa[l - 1] = aa[l - 2];
    est = gl;
    goto L20;

/*     PROCEED TO RIGHT HALF AT THIS LEVEL */

L60:
    vl[l - 1] = glr;
L70:
    est = gr[l - 2];
    lr[l - 1] = 1;
    aa[l - 1] += hh[l - 1] * 4.;
    goto L20;

/*     RETURN ONE LEVEL */

L80:
    vr = glr;
L90:
    if (l <= 1) {
	goto L120;
    }
    --l;
    eps *= 2.;
    ef *= sq2;
    if (lr[l - 1] <= 0) {
	goto L100;
    } else {
	goto L110;
    }
L100:
    vl[l - 1] = vl[l] + vr;
    goto L70;
L110:
    vr = vl[l] + vr;
    goto L90;

/*      EXIT */

L120:
    *ans = vr;
    if (mxl == 0 || abs(be) <= tol * 2. * area) {
	goto L140;
    }
    *ierr = 2;
    xermsg_("SLATEC", "DPPGQ8", "ANS IS PROBABLY INSUFFICIENTLY ACCURATE.", &
	    c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    goto L140;
L130:
    *ierr = -1;
    xermsg_("SLATEC", "DPPGQ8", "A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMA"
	    "L INTEGRATION.  ANSWER IS SET TO ZERO, AND IERR=-1.", &c__1, &
	    c_n1, (ftnlen)6, (ftnlen)6, (ftnlen)94);
L140:
    if (*err < 0.) {
	*err = be;
    }
    return 0;
} /* dppgq8_ */

