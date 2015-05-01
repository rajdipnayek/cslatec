/* dbsgq8.f -- translated by f2c (version 12.02.01).
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

/* DECK DBSGQ8 */
/* Subroutine */ int dbsgq8_(D_fp fun, doublereal *xt, doublereal *bc, 
	integer *n, integer *kk, integer *id, doublereal *a, doublereal *b, 
	integer *inbv, doublereal *err, doublereal *ans, integer *ierr, 
	doublereal *work)
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
    integer i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10, 
	    d__11, d__12, d__13, d__14, d__15, d__16, d__17, d__18;

    /* Local variables */
    static doublereal c__;
    static integer k, l;
    static doublereal aa[60], ae, ce, ee, ef, hh[60], gl, gr[60];
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
    extern doublereal dbvalu_(doublereal *, doublereal *, integer *, integer *
	    , integer *, doublereal *, integer *, doublereal *);
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DBSGQ8 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBFQAD */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (BSGQ8-S, DBSGQ8-D) */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract    **** A DOUBLE PRECISION routine **** */

/*        DBSGQ8, a modification of GAUS8, integrates the */
/*        product of FUN(X) by the ID-th derivative of a spline */
/*        DBVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B. */

/*     Description of Arguments */

/*        INPUT-- FUN,XT,BC,A,B,ERR are DOUBLE PRECISION */
/*        FUN - Name of external function of one argument which */
/*              multiplies DBVALU. */
/*        XT  - Knot array for DBVALU */
/*        BC  - B-coefficient array for DBVALU */
/*        N   - Number of B-coefficients for DBVALU */
/*        KK  - Order of the spline, KK.GE.1 */
/*        ID  - Order of the spline derivative, 0.LE.ID.LE.KK-1 */
/*        A   - Lower limit of integral */
/*        B   - Upper limit of integral (may be less than A) */
/*        INBV- Initialization parameter for DBVALU */
/*        ERR - Is a requested pseudorelative error tolerance.  Normally */
/*              pick a value of ABS(ERR).LT.1D-3.  ANS will normally */
/*              have no more error than ABS(ERR) times the integral of */
/*              the absolute value of FUN(X)*DBVALU(XT,BC,N,KK,X,ID, */
/*              INBV,WORK). */


/*        OUTPUT-- ERR,ANS,WORK are DOUBLE PRECISION */
/*        ERR - Will be an estimate of the absolute error in ANS if the */
/*              input value of ERR was negative.  (ERR is unchanged if */
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
/*        WORK- Work vector of length 3*K for DBVALU */

/* ***SEE ALSO  DBFQAD */
/* ***ROUTINES CALLED  D1MACH, DBVALU, I1MACH, XERMSG */
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
/* ***END PROLOGUE  DBSGQ8 */

    /* Parameter adjustments */
    --work;
    --bc;
    --xt;

    /* Function Body */

/*     INITIALIZE */

/* ***FIRST EXECUTABLE STATEMENT  DBSGQ8 */
    k = i1mach_(&c__14);
    anib = d1mach_(&c__5) * k / .30102;
    nbits = (integer) anib;
/* Computing MIN */
    i__1 = nbits * 5 / 8;
    nlmx = min(i__1,60);
    *ans = 0.;
    *ierr = 1;
    ce = 0.;
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
    c__ = (d__1 = 1. - *a / *b, abs(d__1));
    if (c__ > .1) {
	goto L10;
    }
    if (c__ <= 0.) {
	goto L140;
    }
    anib = .5 - log(c__) / .69314718;
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
    est = d__2 * (w1 * ((*fun)(&d__3) * dbvalu_(&xt[1], &bc[1], n, kk, id, &
	    d__4, inbv, &work[1]) + (*fun)(&d__5) * dbvalu_(&xt[1], &bc[1], n,
	     kk, id, &d__6, inbv, &work[1])) + w2 * ((*fun)(&d__7) * dbvalu_(&
	    xt[1], &bc[1], n, kk, id, &d__8, inbv, &work[1]) + (*fun)(&d__9) *
	     dbvalu_(&xt[1], &bc[1], n, kk, id, &d__10, inbv, &work[1])) + (
	    w3 * ((*fun)(&d__11) * dbvalu_(&xt[1], &bc[1], n, kk, id, &d__12, 
	    inbv, &work[1]) + (*fun)(&d__13) * dbvalu_(&xt[1], &bc[1], n, kk, 
	    id, &d__14, inbv, &work[1])) + w4 * ((*fun)(&d__15) * dbvalu_(&xt[
	    1], &bc[1], n, kk, id, &d__16, inbv, &work[1]) + (*fun)(&d__17) * 
	    dbvalu_(&xt[1], &bc[1], n, kk, id, &d__18, inbv, &work[1]))));
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
    gl = hh[l - 1] * (w1 * ((*fun)(&d__2) * dbvalu_(&xt[1], &bc[1], n, kk, id,
	     &d__3, inbv, &work[1]) + (*fun)(&d__4) * dbvalu_(&xt[1], &bc[1], 
	    n, kk, id, &d__5, inbv, &work[1])) + w2 * ((*fun)(&d__6) * 
	    dbvalu_(&xt[1], &bc[1], n, kk, id, &d__7, inbv, &work[1]) + (*fun)
	    (&d__8) * dbvalu_(&xt[1], &bc[1], n, kk, id, &d__9, inbv, &work[1]
	    )) + (w3 * ((*fun)(&d__10) * dbvalu_(&xt[1], &bc[1], n, kk, id, &
	    d__11, inbv, &work[1]) + (*fun)(&d__12) * dbvalu_(&xt[1], &bc[1], 
	    n, kk, id, &d__13, inbv, &work[1])) + w4 * ((*fun)(&d__14) * 
	    dbvalu_(&xt[1], &bc[1], n, kk, id, &d__15, inbv, &work[1]) + (*
	    fun)(&d__16) * dbvalu_(&xt[1], &bc[1], n, kk, id, &d__17, inbv, &
	    work[1]))));
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
    gr[l - 1] = hh[l - 1] * (w1 * ((*fun)(&d__2) * dbvalu_(&xt[1], &bc[1], n, 
	    kk, id, &d__3, inbv, &work[1]) + (*fun)(&d__4) * dbvalu_(&xt[1], &
	    bc[1], n, kk, id, &d__5, inbv, &work[1])) + w2 * ((*fun)(&d__6) * 
	    dbvalu_(&xt[1], &bc[1], n, kk, id, &d__7, inbv, &work[1]) + (*fun)
	    (&d__8) * dbvalu_(&xt[1], &bc[1], n, kk, id, &d__9, inbv, &work[1]
	    )) + (w3 * ((*fun)(&d__10) * dbvalu_(&xt[1], &bc[1], n, kk, id, &
	    d__11, inbv, &work[1]) + (*fun)(&d__12) * dbvalu_(&xt[1], &bc[1], 
	    n, kk, id, &d__13, inbv, &work[1])) + w4 * ((*fun)(&d__14) * 
	    dbvalu_(&xt[1], &bc[1], n, kk, id, &d__15, inbv, &work[1]) + (*
	    fun)(&d__16) * dbvalu_(&xt[1], &bc[1], n, kk, id, &d__17, inbv, &
	    work[1]))));
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
    ce += est - glr;
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
    if (mxl == 0 || abs(ce) <= tol * 2. * area) {
	goto L140;
    }
    *ierr = 2;
    xermsg_("SLATEC", "DBSGQ8", "ANS IS PROBABLY INSUFFICIENTLY ACCURATE.", &
	    c__3, &c__1, (ftnlen)6, (ftnlen)6, (ftnlen)40);
    goto L140;
L130:
    *ierr = -1;
    xermsg_("SLATEC", "DBSGQ8", "A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMA"
	    "L INTEGRATION.  ANS IS SET TO ZERO AND IERR TO -1.", &c__1, &c_n1,
	     (ftnlen)6, (ftnlen)6, (ftnlen)93);
L140:
    if (*err < 0.) {
	*err = ce;
    }
    return 0;
} /* dbsgq8_ */

