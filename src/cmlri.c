/* cmlri.f -- translated by f2c (version 12.02.01).
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

/* DECK CMLRI */
/* Subroutine */ int cmlri_(complex *z__, real *fnu, integer *kode, integer *
	n, complex *y, integer *nz, real *tol)
{
    /* Initialized data */

    static complex czero = {0.f,0.f};
    static complex cone = {1.f,0.f};
    static complex ctwo = {2.f,0.f};

    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2, r__3;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Local variables */
    static integer i__, k, m;
    static real x;
    static complex p1, p2;
    static real ak, bk;
    static complex ck;
    static real ap, at;
    static integer kk, km;
    static real az;
    static complex pt, rz;
    static real ack, fnf, fkk;
    static integer iaz;
    static real rho;
    static integer inu;
    static complex sum;
    static real tst, rho2, flam, fkap, scle, tfnf;
    static integer idum, ifnu;
    extern doublereal gamln_(real *, integer *);
    static integer itime;
    static complex cnorm;
    extern doublereal r1mach_(integer *);

/* ***BEGIN PROLOGUE  CMLRI */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBESI and CBESK */
/* ***LIBRARY   SLATEC */
/* ***TYPE      ALL (CMLRI-A, ZMLRI-A) */
/* ***AUTHOR  Amos, D. E., (SNL) */
/* ***DESCRIPTION */

/*     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z).GE.0.0 BY THE */
/*     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES. */

/* ***SEE ALSO  CBESI, CBESK */
/* ***ROUTINES CALLED  GAMLN, R1MACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   830501  DATE WRITTEN */
/*   910415  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  CMLRI */
    /* Parameter adjustments */
    --y;

    /* Function Body */
    scle = 1e3f * r1mach_(&c__1) / *tol;
/* ***FIRST EXECUTABLE STATEMENT  CMLRI */
    *nz = 0;
    az = c_abs(z__);
    x = z__->r;
    iaz = az;
    ifnu = *fnu;
    inu = ifnu + *n - 1;
    at = iaz + 1.f;
    q__2.r = at, q__2.i = 0.f;
    c_div(&q__1, &q__2, z__);
    ck.r = q__1.r, ck.i = q__1.i;
    c_div(&q__1, &ctwo, z__);
    rz.r = q__1.r, rz.i = q__1.i;
    p1.r = czero.r, p1.i = czero.i;
    p2.r = cone.r, p2.i = cone.i;
    ack = (at + 1.f) / az;
    rho = ack + sqrt(ack * ack - 1.f);
    rho2 = rho * rho;
    tst = (rho2 + rho2) / ((rho2 - 1.f) * (rho - 1.f));
    tst /= *tol;
/* ----------------------------------------------------------------------- */
/*     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES */
/* ----------------------------------------------------------------------- */
    ak = at;
    for (i__ = 1; i__ <= 80; ++i__) {
	pt.r = p2.r, pt.i = p2.i;
	q__2.r = ck.r * p2.r - ck.i * p2.i, q__2.i = ck.r * p2.i + ck.i * 
		p2.r;
	q__1.r = p1.r - q__2.r, q__1.i = p1.i - q__2.i;
	p2.r = q__1.r, p2.i = q__1.i;
	p1.r = pt.r, p1.i = pt.i;
	q__1.r = ck.r + rz.r, q__1.i = ck.i + rz.i;
	ck.r = q__1.r, ck.i = q__1.i;
	ap = c_abs(&p2);
	if (ap > tst * ak * ak) {
	    goto L20;
	}
	ak += 1.f;
/* L10: */
    }
    goto L110;
L20:
    ++i__;
    k = 0;
    if (inu < iaz) {
	goto L40;
    }
/* ----------------------------------------------------------------------- */
/*     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS */
/* ----------------------------------------------------------------------- */
    p1.r = czero.r, p1.i = czero.i;
    p2.r = cone.r, p2.i = cone.i;
    at = inu + 1.f;
    q__2.r = at, q__2.i = 0.f;
    c_div(&q__1, &q__2, z__);
    ck.r = q__1.r, ck.i = q__1.i;
    ack = at / az;
    tst = sqrt(ack / *tol);
    itime = 1;
    for (k = 1; k <= 80; ++k) {
	pt.r = p2.r, pt.i = p2.i;
	q__2.r = ck.r * p2.r - ck.i * p2.i, q__2.i = ck.r * p2.i + ck.i * 
		p2.r;
	q__1.r = p1.r - q__2.r, q__1.i = p1.i - q__2.i;
	p2.r = q__1.r, p2.i = q__1.i;
	p1.r = pt.r, p1.i = pt.i;
	q__1.r = ck.r + rz.r, q__1.i = ck.i + rz.i;
	ck.r = q__1.r, ck.i = q__1.i;
	ap = c_abs(&p2);
	if (ap < tst) {
	    goto L30;
	}
	if (itime == 2) {
	    goto L40;
	}
	ack = c_abs(&ck);
	flam = ack + sqrt(ack * ack - 1.f);
	fkap = ap / c_abs(&p1);
	rho = dmin(flam,fkap);
	tst *= sqrt(rho / (rho * rho - 1.f));
	itime = 2;
L30:
	;
    }
    goto L110;
L40:
/* ----------------------------------------------------------------------- */
/*     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION */
/* ----------------------------------------------------------------------- */
    ++k;
/* Computing MAX */
    i__1 = i__ + iaz, i__2 = k + inu;
    kk = max(i__1,i__2);
    fkk = (real) kk;
    p1.r = czero.r, p1.i = czero.i;
/* ----------------------------------------------------------------------- */
/*     SCALE P2 AND SUM BY SCLE */
/* ----------------------------------------------------------------------- */
    q__1.r = scle, q__1.i = 0.f;
    p2.r = q__1.r, p2.i = q__1.i;
    fnf = *fnu - ifnu;
    tfnf = fnf + fnf;
    r__1 = fkk + tfnf + 1.f;
    r__2 = fkk + 1.f;
    r__3 = tfnf + 1.f;
    bk = gamln_(&r__1, &idum) - gamln_(&r__2, &idum) - gamln_(&r__3, &idum);
    bk = exp(bk);
    sum.r = czero.r, sum.i = czero.i;
    km = kk - inu;
    i__1 = km;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pt.r = p2.r, pt.i = p2.i;
	r__1 = fkk + fnf;
	q__4.r = r__1, q__4.i = 0.f;
	q__3.r = q__4.r * rz.r - q__4.i * rz.i, q__3.i = q__4.r * rz.i + 
		q__4.i * rz.r;
	q__2.r = q__3.r * p2.r - q__3.i * p2.i, q__2.i = q__3.r * p2.i + 
		q__3.i * p2.r;
	q__1.r = p1.r + q__2.r, q__1.i = p1.i + q__2.i;
	p2.r = q__1.r, p2.i = q__1.i;
	p1.r = pt.r, p1.i = pt.i;
	ak = 1.f - tfnf / (fkk + tfnf);
	ack = bk * ak;
	r__1 = ack + bk;
	q__3.r = r__1, q__3.i = 0.f;
	q__2.r = q__3.r * p1.r - q__3.i * p1.i, q__2.i = q__3.r * p1.i + 
		q__3.i * p1.r;
	q__1.r = sum.r + q__2.r, q__1.i = sum.i + q__2.i;
	sum.r = q__1.r, sum.i = q__1.i;
	bk = ack;
	fkk += -1.f;
/* L50: */
    }
    i__1 = *n;
    y[i__1].r = p2.r, y[i__1].i = p2.i;
    if (*n == 1) {
	goto L70;
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	pt.r = p2.r, pt.i = p2.i;
	r__1 = fkk + fnf;
	q__4.r = r__1, q__4.i = 0.f;
	q__3.r = q__4.r * rz.r - q__4.i * rz.i, q__3.i = q__4.r * rz.i + 
		q__4.i * rz.r;
	q__2.r = q__3.r * p2.r - q__3.i * p2.i, q__2.i = q__3.r * p2.i + 
		q__3.i * p2.r;
	q__1.r = p1.r + q__2.r, q__1.i = p1.i + q__2.i;
	p2.r = q__1.r, p2.i = q__1.i;
	p1.r = pt.r, p1.i = pt.i;
	ak = 1.f - tfnf / (fkk + tfnf);
	ack = bk * ak;
	r__1 = ack + bk;
	q__3.r = r__1, q__3.i = 0.f;
	q__2.r = q__3.r * p1.r - q__3.i * p1.i, q__2.i = q__3.r * p1.i + 
		q__3.i * p1.r;
	q__1.r = sum.r + q__2.r, q__1.i = sum.i + q__2.i;
	sum.r = q__1.r, sum.i = q__1.i;
	bk = ack;
	fkk += -1.f;
	m = *n - i__ + 1;
	i__2 = m;
	y[i__2].r = p2.r, y[i__2].i = p2.i;
/* L60: */
    }
L70:
    if (ifnu <= 0) {
	goto L90;
    }
    i__1 = ifnu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pt.r = p2.r, pt.i = p2.i;
	r__1 = fkk + fnf;
	q__4.r = r__1, q__4.i = 0.f;
	q__3.r = q__4.r * rz.r - q__4.i * rz.i, q__3.i = q__4.r * rz.i + 
		q__4.i * rz.r;
	q__2.r = q__3.r * p2.r - q__3.i * p2.i, q__2.i = q__3.r * p2.i + 
		q__3.i * p2.r;
	q__1.r = p1.r + q__2.r, q__1.i = p1.i + q__2.i;
	p2.r = q__1.r, p2.i = q__1.i;
	p1.r = pt.r, p1.i = pt.i;
	ak = 1.f - tfnf / (fkk + tfnf);
	ack = bk * ak;
	r__1 = ack + bk;
	q__3.r = r__1, q__3.i = 0.f;
	q__2.r = q__3.r * p1.r - q__3.i * p1.i, q__2.i = q__3.r * p1.i + 
		q__3.i * p1.r;
	q__1.r = sum.r + q__2.r, q__1.i = sum.i + q__2.i;
	sum.r = q__1.r, sum.i = q__1.i;
	bk = ack;
	fkk += -1.f;
/* L80: */
    }
L90:
    pt.r = z__->r, pt.i = z__->i;
    if (*kode == 2) {
	q__2.r = x, q__2.i = 0.f;
	q__1.r = pt.r - q__2.r, q__1.i = pt.i - q__2.i;
	pt.r = q__1.r, pt.i = q__1.i;
    }
    q__4.r = fnf, q__4.i = 0.f;
    q__3.r = -q__4.r, q__3.i = -q__4.i;
    c_log(&q__5, &rz);
    q__2.r = q__3.r * q__5.r - q__3.i * q__5.i, q__2.i = q__3.r * q__5.i + 
	    q__3.i * q__5.r;
    q__1.r = q__2.r + pt.r, q__1.i = q__2.i + pt.i;
    p1.r = q__1.r, p1.i = q__1.i;
    r__1 = fnf + 1.f;
    ap = gamln_(&r__1, &idum);
    q__2.r = ap, q__2.i = 0.f;
    q__1.r = p1.r - q__2.r, q__1.i = p1.i - q__2.i;
    pt.r = q__1.r, pt.i = q__1.i;
/* ----------------------------------------------------------------------- */
/*     THE DIVISION CEXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW */
/*     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES */
/* ----------------------------------------------------------------------- */
    q__1.r = p2.r + sum.r, q__1.i = p2.i + sum.i;
    p2.r = q__1.r, p2.i = q__1.i;
    ap = c_abs(&p2);
    r__1 = 1.f / ap;
    q__1.r = r__1, q__1.i = 0.f;
    p1.r = q__1.r, p1.i = q__1.i;
    c_exp(&q__2, &pt);
    q__1.r = q__2.r * p1.r - q__2.i * p1.i, q__1.i = q__2.r * p1.i + q__2.i * 
	    p1.r;
    ck.r = q__1.r, ck.i = q__1.i;
    r_cnjg(&q__2, &p2);
    q__1.r = q__2.r * p1.r - q__2.i * p1.i, q__1.i = q__2.r * p1.i + q__2.i * 
	    p1.r;
    pt.r = q__1.r, pt.i = q__1.i;
    q__1.r = ck.r * pt.r - ck.i * pt.i, q__1.i = ck.r * pt.i + ck.i * pt.r;
    cnorm.r = q__1.r, cnorm.i = q__1.i;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	q__1.r = y[i__3].r * cnorm.r - y[i__3].i * cnorm.i, q__1.i = y[i__3]
		.r * cnorm.i + y[i__3].i * cnorm.r;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L100: */
    }
    return 0;
L110:
    *nz = -2;
    return 0;
} /* cmlri_ */

