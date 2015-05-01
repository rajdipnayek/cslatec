/* blktr1.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    integer npp, k;
    real eps, cnv;
    integer nm, ncmplx, ik;
} cblkt_;

#define cblkt_1 cblkt_

/* Table of constant values */

static integer c__2 = 2;
static integer c__0 = 0;

/* DECK BLKTR1 */
/* Subroutine */ int blktr1_(integer *n, real *an, real *bn, real *cn, 
	integer *m, real *am, real *bm, real *cm, integer *idimy, real *y, 
	real *b, real *w1, real *w2, real *w3, real *wd, real *ww, real *wu, 
	S_fp prdct, S_fp cprdct)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, l, i1, i2, i3, i4, if__, nc, na, ll, ip, ir, np, 
	    iz, nz, im1, im2, im3, ip1, ip2, nm1, nm2, nm3, np1, np2, ip3, 
	    np3, ifd, kdo;
    static real dum;
    static integer izr, imi1, imi2, ipi1, ipi2, ipi3, irm1, idxa, idxc;
    extern /* Subroutine */ int indxa_(integer *, integer *, integer *, 
	    integer *), indxb_(integer *, integer *, integer *, integer *), 
	    indxc_(integer *, integer *, integer *, integer *);

/* ***BEGIN PROLOGUE  BLKTR1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BLKTRI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BLKTR1-S, CBLKT1-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* BLKTR1 solves the linear system set up by BLKTRI. */

/* B  contains the roots of all the B polynomials. */
/* W1,W2,W3,WD,WW,WU  are all working arrays. */
/* PRDCT  is either PRODP or PROD depending on whether the boundary */
/* conditions in the M direction are periodic or not. */
/* CPRDCT is either CPRODP or CPROD which are the complex versions */
/* of PRODP and PROD. These are called in the event that some */
/* of the roots of the B sub P polynomial are complex. */

/* ***SEE ALSO  BLKTRI */
/* ***ROUTINES CALLED  INDXA, INDXB, INDXC */
/* ***COMMON BLOCKS    CBLKT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  BLKTR1 */

/* ***FIRST EXECUTABLE STATEMENT  BLKTR1 */
    /* Parameter adjustments */
    --an;
    --bn;
    --cn;
    --am;
    --bm;
    --cm;
    y_dim1 = *idimy;
    y_offset = 1 + y_dim1;
    y -= y_offset;
    --b;
    --w1;
    --w2;
    --w3;
    --wd;
    --ww;
    --wu;

    /* Function Body */
    kdo = cblkt_1.k - 1;
    i__1 = kdo;
    for (l = 1; l <= i__1; ++l) {
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i3 = i2 + i1;
	i4 = i2 + i2;
	irm1 = ir - 1;
	indxb_(&i2, &ir, &im2, &nm2);
	indxb_(&i1, &irm1, &im3, &nm3);
	indxb_(&i3, &irm1, &im1, &nm1);
	(*prdct)(&nm2, &b[im2], &nm3, &b[im3], &nm1, &b[im1], &c__0, &dum, &y[
		i2 * y_dim1 + 1], &w3[1], m, &am[1], &bm[1], &cm[1], &wd[1], &
		ww[1], &wu[1]);
	if__ = pow_ii(&c__2, &cblkt_1.k);
	i__2 = if__;
	i__3 = i4;
	for (i__ = i4; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    if (i__ - cblkt_1.nm <= 0) {
		goto L101;
	    } else {
		goto L108;
	    }
L101:
	    ipi1 = i__ + i1;
	    ipi2 = i__ + i2;
	    ipi3 = i__ + i3;
	    indxc_(&i__, &ir, &idxc, &nc);
	    if (i__ - if__ >= 0) {
		goto L108;
	    } else {
		goto L102;
	    }
L102:
	    indxa_(&i__, &ir, &idxa, &na);
	    i__4 = i__ - i1;
	    indxb_(&i__4, &irm1, &im1, &nm1);
	    indxb_(&ipi2, &ir, &ip2, &np2);
	    indxb_(&ipi1, &irm1, &ip1, &np1);
	    indxb_(&ipi3, &irm1, &ip3, &np3);
	    (*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], 
		    &w3[1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1],
		     &wu[1]);
	    if (ipi2 - cblkt_1.nm <= 0) {
		goto L105;
	    } else {
		goto L103;
	    }
L103:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		w3[j] = 0.f;
		w2[j] = 0.f;
/* L104: */
	    }
	    goto L106;
L105:
	    (*prdct)(&np2, &b[ip2], &np1, &b[ip1], &np3, &b[ip3], &c__0, &dum,
		     &y[ipi2 * y_dim1 + 1], &w3[1], m, &am[1], &bm[1], &cm[1],
		     &wd[1], &ww[1], &wu[1]);
	    (*prdct)(&np1, &b[ip1], &c__0, &dum, &c__0, &dum, &nc, &cn[idxc], 
		    &w3[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1],
		     &wu[1]);
L106:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		y[j + i__ * y_dim1] = w1[j] + w2[j] + y[j + i__ * y_dim1];
/* L107: */
	    }
L108:
	    ;
	}
/* L109: */
    }
    if (cblkt_1.npp != 0) {
	goto L132;
    } else {
	goto L110;
    }

/*     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD */

L110:
    if__ = pow_ii(&c__2, &cblkt_1.k);
    i__ = if__ / 2;
    i1 = i__ / 2;
    i__1 = i__ - i1;
    i__3 = cblkt_1.k - 2;
    indxb_(&i__1, &i__3, &im1, &nm1);
    i__1 = i__ + i1;
    i__3 = cblkt_1.k - 2;
    indxb_(&i__1, &i__3, &ip1, &np1);
    i__1 = cblkt_1.k - 1;
    indxb_(&i__, &i__1, &iz, &nz);
    (*prdct)(&nz, &b[iz], &nm1, &b[im1], &np1, &b[ip1], &c__0, &dum, &y[i__ * 
	    y_dim1 + 1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &
	    wu[1]);
    izr = i__;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	w2[j] = w1[j];
/* L111: */
    }
    i__1 = cblkt_1.k;
    for (ll = 2; ll <= i__1; ++ll) {
	l = cblkt_1.k - ll + 1;
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i__ = i2;
	indxc_(&i__, &ir, &idxc, &nc);
	indxb_(&i__, &ir, &iz, &nz);
	i__3 = i__ - i1;
	i__2 = ir - 1;
	indxb_(&i__3, &i__2, &im1, &nm1);
	i__3 = i__ + i1;
	i__2 = ir - 1;
	indxb_(&i__3, &i__2, &ip1, &np1);
	(*prdct)(&np1, &b[ip1], &c__0, &dum, &c__0, &dum, &nc, &cn[idxc], &w1[
		1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1])
		;
	i__3 = *m;
	for (j = 1; j <= i__3; ++j) {
	    w1[j] = y[j + i__ * y_dim1] + w1[j];
/* L112: */
	}
	(*prdct)(&nz, &b[iz], &nm1, &b[im1], &np1, &b[ip1], &c__0, &dum, &w1[
		1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1])
		;
/* L113: */
    }
    i__1 = cblkt_1.k;
    for (ll = 2; ll <= i__1; ++ll) {
	l = cblkt_1.k - ll + 1;
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i4 = i2 + i2;
	ifd = if__ - i2;
	i__3 = ifd;
	i__2 = i4;
	for (i__ = i2; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__2) {
	    if (i__ - i2 - izr != 0) {
		goto L117;
	    } else {
		goto L114;
	    }
L114:
	    if (i__ - cblkt_1.nm <= 0) {
		goto L115;
	    } else {
		goto L118;
	    }
L115:
	    indxa_(&i__, &ir, &idxa, &na);
	    indxb_(&i__, &ir, &iz, &nz);
	    i__4 = i__ - i1;
	    i__5 = ir - 1;
	    indxb_(&i__4, &i__5, &im1, &nm1);
	    i__4 = i__ + i1;
	    i__5 = ir - 1;
	    indxb_(&i__4, &i__5, &ip1, &np1);
	    (*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], 
		    &w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1],
		     &wu[1]);
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		w2[j] = y[j + i__ * y_dim1] + w2[j];
/* L116: */
	    }
	    (*prdct)(&nz, &b[iz], &nm1, &b[im1], &np1, &b[ip1], &c__0, &dum, &
		    w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], 
		    &wu[1]);
	    izr = i__;
	    if (i__ - cblkt_1.nm != 0) {
		goto L117;
	    } else {
		goto L119;
	    }
L117:
	    ;
	}
L118:
	;
    }
L119:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	y[j + (cblkt_1.nm + 1) * y_dim1] = y[j + (cblkt_1.nm + 1) * y_dim1] - 
		cn[cblkt_1.nm + 1] * w1[j] - an[cblkt_1.nm + 1] * w2[j];
/* L120: */
    }
    i__1 = if__ / 2;
    i__2 = cblkt_1.k - 1;
    indxb_(&i__1, &i__2, &im1, &nm1);
    i__1 = cblkt_1.k - 1;
    indxb_(&if__, &i__1, &ip, &np);
    if (cblkt_1.ncmplx != 0) {
	goto L121;
    } else {
	goto L122;
    }
L121:
    i__1 = cblkt_1.nm + 1;
    (*cprdct)(&i__1, &b[ip], &nm1, &b[im1], &c__0, &dum, &c__0, &dum, &y[(
	    cblkt_1.nm + 1) * y_dim1 + 1], &y[(cblkt_1.nm + 1) * y_dim1 + 1], 
	    m, &am[1], &bm[1], &cm[1], &w1[1], &w3[1], &ww[1]);
    goto L123;
L122:
    i__1 = cblkt_1.nm + 1;
    (*prdct)(&i__1, &b[ip], &nm1, &b[im1], &c__0, &dum, &c__0, &dum, &y[(
	    cblkt_1.nm + 1) * y_dim1 + 1], &y[(cblkt_1.nm + 1) * y_dim1 + 1], 
	    m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1]);
L123:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	w1[j] = an[1] * y[j + (cblkt_1.nm + 1) * y_dim1];
	w2[j] = cn[cblkt_1.nm] * y[j + (cblkt_1.nm + 1) * y_dim1];
	y[j + y_dim1] -= w1[j];
	y[j + cblkt_1.nm * y_dim1] -= w2[j];
/* L124: */
    }
    i__1 = kdo;
    for (l = 1; l <= i__1; ++l) {
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i4 = i2 + i2;
	i1 = i2 / 2;
	i__ = i4;
	indxa_(&i__, &ir, &idxa, &na);
	i__2 = i__ - i2;
	indxb_(&i__2, &ir, &im2, &nm2);
	i__2 = i__ - i2 - i1;
	i__3 = ir - 1;
	indxb_(&i__2, &i__3, &im3, &nm3);
	i__2 = i__ - i1;
	i__3 = ir - 1;
	indxb_(&i__2, &i__3, &im1, &nm1);
	(*prdct)(&nm2, &b[im2], &nm3, &b[im3], &nm1, &b[im1], &c__0, &dum, &
		w1[1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[
		1]);
	(*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], &w1[
		1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1])
		;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    y[j + i__ * y_dim1] -= w1[j];
/* L125: */
	}
/* L126: */
    }

    izr = cblkt_1.nm;
    i__1 = kdo;
    for (l = 1; l <= i__1; ++l) {
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i3 = i2 + i1;
	i4 = i2 + i2;
	irm1 = ir - 1;
	i__2 = if__;
	i__3 = i4;
	for (i__ = i4; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    ipi1 = i__ + i1;
	    ipi2 = i__ + i2;
	    ipi3 = i__ + i3;
	    if (ipi2 - izr != 0) {
		goto L127;
	    } else {
		goto L128;
	    }
L127:
	    if (i__ - izr != 0) {
		goto L130;
	    } else {
		goto L131;
	    }
L128:
	    indxc_(&i__, &ir, &idxc, &nc);
	    indxb_(&ipi2, &ir, &ip2, &np2);
	    indxb_(&ipi1, &irm1, &ip1, &np1);
	    indxb_(&ipi3, &irm1, &ip3, &np3);
	    (*prdct)(&np2, &b[ip2], &np1, &b[ip1], &np3, &b[ip3], &c__0, &dum,
		     &w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1]
		    , &wu[1]);
	    (*prdct)(&np1, &b[ip1], &c__0, &dum, &c__0, &dum, &nc, &cn[idxc], 
		    &w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1],
		     &wu[1]);
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		y[j + i__ * y_dim1] -= w2[j];
/* L129: */
	    }
	    izr = i__;
	    goto L131;
L130:
	    ;
	}
L131:
	;
    }

/* BEGIN BACK SUBSTITUTION PHASE */

L132:
    i__1 = cblkt_1.k;
    for (ll = 1; ll <= i__1; ++ll) {
	l = cblkt_1.k - ll + 1;
	ir = l - 1;
	irm1 = ir - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i4 = i2 + i2;
	ifd = if__ - i2;
	i__3 = ifd;
	i__2 = i4;
	for (i__ = i2; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__2) {
	    if (i__ - cblkt_1.nm <= 0) {
		goto L133;
	    } else {
		goto L143;
	    }
L133:
	    imi1 = i__ - i1;
	    imi2 = i__ - i2;
	    ipi1 = i__ + i1;
	    ipi2 = i__ + i2;
	    indxa_(&i__, &ir, &idxa, &na);
	    indxc_(&i__, &ir, &idxc, &nc);
	    indxb_(&i__, &ir, &iz, &nz);
	    indxb_(&imi1, &irm1, &im1, &nm1);
	    indxb_(&ipi1, &irm1, &ip1, &np1);
	    if (i__ - i2 <= 0) {
		goto L134;
	    } else {
		goto L136;
	    }
L134:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		w1[j] = 0.f;
/* L135: */
	    }
	    goto L137;
L136:
	    (*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], 
		    &y[imi2 * y_dim1 + 1], &w1[1], m, &am[1], &bm[1], &cm[1], 
		    &wd[1], &ww[1], &wu[1]);
L137:
	    if (ipi2 - cblkt_1.nm <= 0) {
		goto L140;
	    } else {
		goto L138;
	    }
L138:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		w2[j] = 0.f;
/* L139: */
	    }
	    goto L141;
L140:
	    (*prdct)(&np1, &b[ip1], &c__0, &dum, &c__0, &dum, &nc, &cn[idxc], 
		    &y[ipi2 * y_dim1 + 1], &w2[1], m, &am[1], &bm[1], &cm[1], 
		    &wd[1], &ww[1], &wu[1]);
L141:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		w1[j] = y[j + i__ * y_dim1] + w1[j] + w2[j];
/* L142: */
	    }
	    (*prdct)(&nz, &b[iz], &nm1, &b[im1], &np1, &b[ip1], &c__0, &dum, &
		    w1[1], &y[i__ * y_dim1 + 1], m, &am[1], &bm[1], &cm[1], &
		    wd[1], &ww[1], &wu[1]);
L143:
	    ;
	}
/* L144: */
    }
    return 0;
} /* blktr1_ */

