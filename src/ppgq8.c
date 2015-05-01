/* ppgq8.f -- translated by f2c (version 12.02.01).
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

static integer c__11 = 11;
static integer c__5 = 5;
static real c_b6 = 1.f;
static real c_b8 = 2.f;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c_n1 = -1;

/* DECK PPGQ8 */
/* Subroutine */ int ppgq8_(E_fp fun, integer *ldc, real *c__, real *xi, 
	integer *lxi, integer *kk, integer *id, real *a, real *b, integer *
	inppv, real *err, real *ans, integer *ierr)
{
    /* Initialized data */

    static real x1 = .183434642495649805f;
    static integer nlmn = 1;
    static integer kmx = 5000;
    static integer kml = 6;
    static real x2 = .525532409916328986f;
    static real x3 = .79666647741362674f;
    static real x4 = .960289856497536232f;
    static real w1 = .362683783378361983f;
    static real w2 = .313706645877887287f;
    static real w3 = .222381034453374471f;
    static real w4 = .101228536290376259f;
    static real sq2 = 1.41421356f;

    /* System generated locals */
    integer c_dim1, c_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11, 
	    r__12, r__13, r__14, r__15, r__16, r__17, r__18;

    /* Local variables */
    static integer k, l;
    static real aa[30], ae, be, cc, ee, ef, hh[30], gl, gr[30];
    static integer lr[30];
    static real vl[30], vr;
    static integer nib;
    static real glr;
    static integer lmn;
    static real eps, est, tol;
    static integer lmx, mxl;
    static real area, anib;
    static integer nlmx, nbits;
    extern doublereal ppval_(integer *, real *, real *, integer *, integer *, 
	    integer *, real *, integer *);
    extern integer i1mach_(integer *);
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  PPGQ8 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to PFQAD */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PPGQ8-S, DPPGQ8-D) */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        PPGQ8, a modification of GAUS8, integrates the */
/*        product of FUN(X) by the ID-th derivative of a spline */
/*        PPVAL(LDC,C,XI,LXI,KK,ID,X,INPPV)  between limits A and B. */

/*     Description of arguments */

/*        INPUT-- */
/*        FUN - Name of external function of one argument which */
/*              multiplies PPVAL. */
/*        LDC - Leading dimension of matrix C, LDC.GE.KK */
/*        C   - Matrix of Taylor derivatives of dimension at least */
/*              (K,LXI) */
/*        XI  - Breakpoint vector of length LXI+1 */
/*        LXI - Number of polynomial pieces */
/*        KK  - Order of the spline, KK.GE.1 */
/*        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1 */
/*        A   - Lower limit of integral */
/*        B   - Upper limit of integral (may be less than A) */
/*        INPPV- Initialization parameter for PPVAL */
/*        ERR - Is a requested pseudorelative error tolerance.  Normally */
/*              pick a value of ABS(ERR).LT.1E-3.  ANS will normally */
/*              have no more error than ABS(ERR) times the integral of */
/*              the absolute value of FUN(X)*PPVAL(LDC,C,XI,LXI,KK,ID,X, */
/*              INPPV). */

/*        OUTPUT-- */
/*        ERR - Will be an estimate of the absolute error in ANS if the */
/*              input value of ERR was negative.  (ERR is unchanged if */
/*              the input value of ERR was nonnegative.)  The estimated */
/*              error is solely for information to the user and should */
/*              not be used as a correction to the computed integral. */
/*        ANS - Computed value of integral */
/*        IERR- A status code */
/*            --Normal codes */
/*               1 ANS most likely meets requested error tolerance, */
/*                 or A=B. */
/*              -1 A and B ARE too nearly equal to allow normal */
/*                 integration.  ANS is set to zero. */
/*            --Abnormal code */
/*               2 ANS probably does not meet requested error tolerance. */

/* ***SEE ALSO  PFQAD */
/* ***ROUTINES CALLED  I1MACH, PPVAL, R1MACH, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800901  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910408  Updated the AUTHOR section.  (WRB) */
/* ***END PROLOGUE  PPGQ8 */

    /* Parameter adjustments */
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --xi;

    /* Function Body */

/*     INITIALIZE */

/* ***FIRST EXECUTABLE STATEMENT  PPGQ8 */
    k = i1mach_(&c__11);
    anib = r1mach_(&c__5) * k / .30102f;
    nbits = (integer) anib;
    nlmx = nbits * 5 / 8;
    *ans = 0.f;
    *ierr = 1;
    be = 0.f;
    if (*a == *b) {
	goto L140;
    }
    lmx = nlmx;
    lmn = nlmn;
    if (*b == 0.f) {
	goto L10;
    }
    if (r_sign(&c_b6, b) * *a <= 0.f) {
	goto L10;
    }
    cc = (r__1 = 1.f - *a / *b, dabs(r__1));
    if (cc > .1f) {
	goto L10;
    }
    if (cc <= 0.f) {
	goto L140;
    }
    anib = .5f - log(cc) / .69314718f;
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
    r__1 = dabs(*err), r__2 = pow_ri(&c_b8, &i__1);
    tol = dmax(r__1,r__2) / 2.f;
    if (*err == 0.f) {
	tol = sqrt(r1mach_(&c__4));
    }
    eps = tol;
    hh[0] = (*b - *a) / 4.f;
    aa[0] = *a;
    lr[0] = 1;
    l = 1;
    r__1 = aa[l - 1] + hh[l - 1] * 2.f;
    r__2 = hh[l - 1] * 2.f;
    r__3 = r__1 - x1 * r__2;
    r__4 = r__1 - x1 * r__2;
    r__5 = r__1 + x1 * r__2;
    r__6 = r__1 + x1 * r__2;
    r__7 = r__1 - x2 * r__2;
    r__8 = r__1 - x2 * r__2;
    r__9 = r__1 + x2 * r__2;
    r__10 = r__1 + x2 * r__2;
    r__11 = r__1 - x3 * r__2;
    r__12 = r__1 - x3 * r__2;
    r__13 = r__1 + x3 * r__2;
    r__14 = r__1 + x3 * r__2;
    r__15 = r__1 - x4 * r__2;
    r__16 = r__1 - x4 * r__2;
    r__17 = r__1 + x4 * r__2;
    r__18 = r__1 + x4 * r__2;
    est = r__2 * (w1 * ((*fun)(&r__3) * ppval_(ldc, &c__[c_offset], &xi[1], 
	    lxi, kk, id, &r__4, inppv) + (*fun)(&r__5) * ppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &r__6, inppv)) + w2 * ((*fun)(&
	    r__7) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &r__8, 
	    inppv) + (*fun)(&r__9) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, 
	    kk, id, &r__10, inppv)) + (w3 * ((*fun)(&r__11) * ppval_(ldc, &
	    c__[c_offset], &xi[1], lxi, kk, id, &r__12, inppv) + (*fun)(&
	    r__13) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &r__14, 
	    inppv)) + w4 * ((*fun)(&r__15) * ppval_(ldc, &c__[c_offset], &xi[
	    1], lxi, kk, id, &r__16, inppv) + (*fun)(&r__17) * ppval_(ldc, &
	    c__[c_offset], &xi[1], lxi, kk, id, &r__18, inppv))));
    k = 8;
    area = dabs(est);
    ef = .5f;
    mxl = 0;

/*     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC. */

L20:
    r__1 = aa[l - 1] + hh[l - 1];
    r__2 = r__1 - x1 * hh[l - 1];
    r__3 = r__1 - x1 * hh[l - 1];
    r__4 = r__1 + x1 * hh[l - 1];
    r__5 = r__1 + x1 * hh[l - 1];
    r__6 = r__1 - x2 * hh[l - 1];
    r__7 = r__1 - x2 * hh[l - 1];
    r__8 = r__1 + x2 * hh[l - 1];
    r__9 = r__1 + x2 * hh[l - 1];
    r__10 = r__1 - x3 * hh[l - 1];
    r__11 = r__1 - x3 * hh[l - 1];
    r__12 = r__1 + x3 * hh[l - 1];
    r__13 = r__1 + x3 * hh[l - 1];
    r__14 = r__1 - x4 * hh[l - 1];
    r__15 = r__1 - x4 * hh[l - 1];
    r__16 = r__1 + x4 * hh[l - 1];
    r__17 = r__1 + x4 * hh[l - 1];
    gl = hh[l - 1] * (w1 * ((*fun)(&r__2) * ppval_(ldc, &c__[c_offset], &xi[1]
	    , lxi, kk, id, &r__3, inppv) + (*fun)(&r__4) * ppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &r__5, inppv)) + w2 * ((*fun)(&
	    r__6) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &r__7, 
	    inppv) + (*fun)(&r__8) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, 
	    kk, id, &r__9, inppv)) + (w3 * ((*fun)(&r__10) * ppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &r__11, inppv) + (*fun)(&r__12) * 
	    ppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &r__13, inppv)) 
	    + w4 * ((*fun)(&r__14) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, 
	    kk, id, &r__15, inppv) + (*fun)(&r__16) * ppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &r__17, inppv))));
    r__1 = aa[l - 1] + hh[l - 1] * 3.f;
    r__2 = r__1 - x1 * hh[l - 1];
    r__3 = r__1 - x1 * hh[l - 1];
    r__4 = r__1 + x1 * hh[l - 1];
    r__5 = r__1 + x1 * hh[l - 1];
    r__6 = r__1 - x2 * hh[l - 1];
    r__7 = r__1 - x2 * hh[l - 1];
    r__8 = r__1 + x2 * hh[l - 1];
    r__9 = r__1 + x2 * hh[l - 1];
    r__10 = r__1 - x3 * hh[l - 1];
    r__11 = r__1 - x3 * hh[l - 1];
    r__12 = r__1 + x3 * hh[l - 1];
    r__13 = r__1 + x3 * hh[l - 1];
    r__14 = r__1 - x4 * hh[l - 1];
    r__15 = r__1 - x4 * hh[l - 1];
    r__16 = r__1 + x4 * hh[l - 1];
    r__17 = r__1 + x4 * hh[l - 1];
    gr[l - 1] = hh[l - 1] * (w1 * ((*fun)(&r__2) * ppval_(ldc, &c__[c_offset],
	     &xi[1], lxi, kk, id, &r__3, inppv) + (*fun)(&r__4) * ppval_(ldc, 
	    &c__[c_offset], &xi[1], lxi, kk, id, &r__5, inppv)) + w2 * ((*fun)
	    (&r__6) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &r__7, 
	    inppv) + (*fun)(&r__8) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, 
	    kk, id, &r__9, inppv)) + (w3 * ((*fun)(&r__10) * ppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &r__11, inppv) + (*fun)(&r__12) * 
	    ppval_(ldc, &c__[c_offset], &xi[1], lxi, kk, id, &r__13, inppv)) 
	    + w4 * ((*fun)(&r__14) * ppval_(ldc, &c__[c_offset], &xi[1], lxi, 
	    kk, id, &r__15, inppv) + (*fun)(&r__16) * ppval_(ldc, &c__[
	    c_offset], &xi[1], lxi, kk, id, &r__17, inppv))));
    k += 16;
    area += dabs(gl) + (r__1 = gr[l - 1], dabs(r__1)) - dabs(est);
    glr = gl + gr[l - 1];
    ee = (r__1 = est - glr, dabs(r__1)) * ef;
/* Computing MAX */
    r__1 = eps * area, r__2 = tol * dabs(glr);
    ae = dmax(r__1,r__2);
    if (ee - ae <= 0.f) {
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
    eps *= .5f;
    ef /= sq2;
    hh[l - 1] = hh[l - 2] * .5f;
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
    aa[l - 1] += hh[l - 1] * 4.f;
    goto L20;

/*     RETURN ONE LEVEL */

L80:
    vr = glr;
L90:
    if (l <= 1) {
	goto L120;
    }
    --l;
    eps *= 2.f;
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
    if (mxl == 0 || dabs(be) <= tol * 2.f * area) {
	goto L140;
    }
    *ierr = 2;
    xermsg_("SLATEC", "PPGQ8", "ANS IS PROBABLY INSUFFICIENTLY ACCURATE.", &
	    c__3, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)40);
    goto L140;
L130:
    *ierr = -1;
    xermsg_("SLATEC", "PPGQ8", "A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL"
	    " INTEGRATION. ANS IS SET TO ZERO AND IERR TO -1.", &c__1, &c_n1, (
	    ftnlen)6, (ftnlen)5, (ftnlen)92);
L140:
    if (*err < 0.f) {
	*err = be;
    }
    return 0;
} /* ppgq8_ */

