/* cblkt1.f -- translated by f2c (version 12.02.01).
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
} ccblk_;

#define ccblk_1 ccblk_

/* Table of constant values */

static integer c__2 = 2;
static integer c__0 = 0;

/* DECK CBLKT1 */
/* Subroutine */ int cblkt1_(integer *n, real *an, real *bn, real *cn, 
	integer *m, complex *am, complex *bm, complex *cm, integer *idimy, 
	complex *y, real *b, complex *w1, complex *w2, complex *w3, complex *
	wd, complex *ww, complex *wu, S_fp prdct, S_fp cprdct)
{
    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, j, l, i1, i2, i3, i4, if__, nc, na, ll, ip, ir, np, 
	    iz, nz, im1, im2, im3, ip1, ip2, nm1, nm2, nm3, np1, np2, ip3, 
	    np3, ifd, kdo;
    static real dum;
    static integer izr, imi1, imi2, ipi1, ipi2, ipi3, irm1, idxa, idxc;
    extern /* Subroutine */ int inxca_(integer *, integer *, integer *, 
	    integer *), inxcb_(integer *, integer *, integer *, integer *), 
	    inxcc_(integer *, integer *, integer *, integer *);

/* ***BEGIN PROLOGUE  CBLKT1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (BLKTR1-S, CBLKT1-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* CBLKT1 solves the linear system of routine CBLKTR. */

/* B  contains the roots of all the B polynomials. */
/* W1,W2,W3,WD,WW,WU  are all working arrays. */
/* PRDCT is either PROCP or PROC depending on whether the boundary */
/* conditions in the M direction are periodic or not. */
/* CPRDCT is either CPROCP or CPROC which are called if some of the zeros */
/* of the B polynomials are complex. */

/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  INXCA, INXCB, INXCC */
/* ***COMMON BLOCKS    CCBLK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CBLKT1 */

/* ***FIRST EXECUTABLE STATEMENT  CBLKT1 */
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
    kdo = ccblk_1.k - 1;
    i__1 = kdo;
    for (l = 1; l <= i__1; ++l) {
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i3 = i2 + i1;
	i4 = i2 + i2;
	irm1 = ir - 1;
	inxcb_(&i2, &ir, &im2, &nm2);
	inxcb_(&i1, &irm1, &im3, &nm3);
	inxcb_(&i3, &irm1, &im1, &nm1);
	(*prdct)(&nm2, &b[im2], &nm3, &b[im3], &nm1, &b[im1], &c__0, &dum, &y[
		i2 * y_dim1 + 1], &w3[1], m, &am[1], &bm[1], &cm[1], &wd[1], &
		ww[1], &wu[1]);
	if__ = pow_ii(&c__2, &ccblk_1.k);
	i__2 = if__;
	i__3 = i4;
	for (i__ = i4; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {
	    if (i__ - ccblk_1.nm <= 0) {
		goto L101;
	    } else {
		goto L108;
	    }
L101:
	    ipi1 = i__ + i1;
	    ipi2 = i__ + i2;
	    ipi3 = i__ + i3;
	    inxcc_(&i__, &ir, &idxc, &nc);
	    if (i__ - if__ >= 0) {
		goto L108;
	    } else {
		goto L102;
	    }
L102:
	    inxca_(&i__, &ir, &idxa, &na);
	    i__4 = i__ - i1;
	    inxcb_(&i__4, &irm1, &im1, &nm1);
	    inxcb_(&ipi2, &ir, &ip2, &np2);
	    inxcb_(&ipi1, &irm1, &ip1, &np1);
	    inxcb_(&ipi3, &irm1, &ip3, &np3);
	    (*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], 
		    &w3[1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1],
		     &wu[1]);
	    if (ipi2 - ccblk_1.nm <= 0) {
		goto L105;
	    } else {
		goto L103;
	    }
L103:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		i__5 = j;
		w3[i__5].r = 0.f, w3[i__5].i = 0.f;
		i__5 = j;
		w2[i__5].r = 0.f, w2[i__5].i = 0.f;
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
		i__5 = j + i__ * y_dim1;
		i__6 = j;
		i__7 = j;
		q__2.r = w1[i__6].r + w2[i__7].r, q__2.i = w1[i__6].i + w2[
			i__7].i;
		i__8 = j + i__ * y_dim1;
		q__1.r = q__2.r + y[i__8].r, q__1.i = q__2.i + y[i__8].i;
		y[i__5].r = q__1.r, y[i__5].i = q__1.i;
/* L107: */
	    }
L108:
	    ;
	}
/* L109: */
    }
    if (ccblk_1.npp != 0) {
	goto L132;
    } else {
	goto L110;
    }

/*     THE PERIODIC CASE IS TREATED USING THE CAPACITANCE MATRIX METHOD */

L110:
    if__ = pow_ii(&c__2, &ccblk_1.k);
    i__ = if__ / 2;
    i1 = i__ / 2;
    i__1 = i__ - i1;
    i__3 = ccblk_1.k - 2;
    inxcb_(&i__1, &i__3, &im1, &nm1);
    i__1 = i__ + i1;
    i__3 = ccblk_1.k - 2;
    inxcb_(&i__1, &i__3, &ip1, &np1);
    i__1 = ccblk_1.k - 1;
    inxcb_(&i__, &i__1, &iz, &nz);
    (*prdct)(&nz, &b[iz], &nm1, &b[im1], &np1, &b[ip1], &c__0, &dum, &y[i__ * 
	    y_dim1 + 1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &
	    wu[1]);
    izr = i__;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__3 = j;
	i__2 = j;
	w2[i__3].r = w1[i__2].r, w2[i__3].i = w1[i__2].i;
/* L111: */
    }
    i__1 = ccblk_1.k;
    for (ll = 2; ll <= i__1; ++ll) {
	l = ccblk_1.k - ll + 1;
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i__ = i2;
	inxcc_(&i__, &ir, &idxc, &nc);
	inxcb_(&i__, &ir, &iz, &nz);
	i__3 = i__ - i1;
	i__2 = ir - 1;
	inxcb_(&i__3, &i__2, &im1, &nm1);
	i__3 = i__ + i1;
	i__2 = ir - 1;
	inxcb_(&i__3, &i__2, &ip1, &np1);
	(*prdct)(&np1, &b[ip1], &c__0, &dum, &c__0, &dum, &nc, &cn[idxc], &w1[
		1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1])
		;
	i__3 = *m;
	for (j = 1; j <= i__3; ++j) {
	    i__2 = j;
	    i__4 = j + i__ * y_dim1;
	    i__5 = j;
	    q__1.r = y[i__4].r + w1[i__5].r, q__1.i = y[i__4].i + w1[i__5].i;
	    w1[i__2].r = q__1.r, w1[i__2].i = q__1.i;
/* L112: */
	}
	(*prdct)(&nz, &b[iz], &nm1, &b[im1], &np1, &b[ip1], &c__0, &dum, &w1[
		1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1])
		;
/* L113: */
    }
    i__1 = ccblk_1.k;
    for (ll = 2; ll <= i__1; ++ll) {
	l = ccblk_1.k - ll + 1;
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
	    if (i__ - ccblk_1.nm <= 0) {
		goto L115;
	    } else {
		goto L118;
	    }
L115:
	    inxca_(&i__, &ir, &idxa, &na);
	    inxcb_(&i__, &ir, &iz, &nz);
	    i__4 = i__ - i1;
	    i__5 = ir - 1;
	    inxcb_(&i__4, &i__5, &im1, &nm1);
	    i__4 = i__ + i1;
	    i__5 = ir - 1;
	    inxcb_(&i__4, &i__5, &ip1, &np1);
	    (*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], 
		    &w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1],
		     &wu[1]);
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		i__5 = j;
		i__6 = j + i__ * y_dim1;
		i__7 = j;
		q__1.r = y[i__6].r + w2[i__7].r, q__1.i = y[i__6].i + w2[i__7]
			.i;
		w2[i__5].r = q__1.r, w2[i__5].i = q__1.i;
/* L116: */
	    }
	    (*prdct)(&nz, &b[iz], &nm1, &b[im1], &np1, &b[ip1], &c__0, &dum, &
		    w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], 
		    &wu[1]);
	    izr = i__;
	    if (i__ - ccblk_1.nm != 0) {
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
	i__2 = j + (ccblk_1.nm + 1) * y_dim1;
	i__3 = j + (ccblk_1.nm + 1) * y_dim1;
	i__4 = ccblk_1.nm + 1;
	i__5 = j;
	q__3.r = cn[i__4] * w1[i__5].r, q__3.i = cn[i__4] * w1[i__5].i;
	q__2.r = y[i__3].r - q__3.r, q__2.i = y[i__3].i - q__3.i;
	i__6 = ccblk_1.nm + 1;
	i__7 = j;
	q__4.r = an[i__6] * w2[i__7].r, q__4.i = an[i__6] * w2[i__7].i;
	q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L120: */
    }
    i__1 = if__ / 2;
    i__2 = ccblk_1.k - 1;
    inxcb_(&i__1, &i__2, &im1, &nm1);
    i__1 = ccblk_1.k - 1;
    inxcb_(&if__, &i__1, &ip, &np);
    if (ccblk_1.ncmplx != 0) {
	goto L121;
    } else {
	goto L122;
    }
L121:
    i__1 = ccblk_1.nm + 1;
    (*cprdct)(&i__1, &b[ip], &nm1, &b[im1], &c__0, &dum, &c__0, &dum, &y[(
	    ccblk_1.nm + 1) * y_dim1 + 1], &y[(ccblk_1.nm + 1) * y_dim1 + 1], 
	    m, &am[1], &bm[1], &cm[1], &w1[1], &w3[1], &ww[1]);
    goto L123;
L122:
    i__1 = ccblk_1.nm + 1;
    (*prdct)(&i__1, &b[ip], &nm1, &b[im1], &c__0, &dum, &c__0, &dum, &y[(
	    ccblk_1.nm + 1) * y_dim1 + 1], &y[(ccblk_1.nm + 1) * y_dim1 + 1], 
	    m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1]);
L123:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j + (ccblk_1.nm + 1) * y_dim1;
	q__1.r = an[1] * y[i__3].r, q__1.i = an[1] * y[i__3].i;
	w1[i__2].r = q__1.r, w1[i__2].i = q__1.i;
	i__2 = j;
	i__3 = ccblk_1.nm;
	i__4 = j + (ccblk_1.nm + 1) * y_dim1;
	q__1.r = cn[i__3] * y[i__4].r, q__1.i = cn[i__3] * y[i__4].i;
	w2[i__2].r = q__1.r, w2[i__2].i = q__1.i;
	i__2 = j + y_dim1;
	i__3 = j + y_dim1;
	i__4 = j;
	q__1.r = y[i__3].r - w1[i__4].r, q__1.i = y[i__3].i - w1[i__4].i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
	i__2 = j + ccblk_1.nm * y_dim1;
	i__3 = j + ccblk_1.nm * y_dim1;
	i__4 = j;
	q__1.r = y[i__3].r - w2[i__4].r, q__1.i = y[i__3].i - w2[i__4].i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L124: */
    }
    i__1 = kdo;
    for (l = 1; l <= i__1; ++l) {
	ir = l - 1;
	i2 = pow_ii(&c__2, &ir);
	i4 = i2 + i2;
	i1 = i2 / 2;
	i__ = i4;
	inxca_(&i__, &ir, &idxa, &na);
	i__2 = i__ - i2;
	inxcb_(&i__2, &ir, &im2, &nm2);
	i__2 = i__ - i2 - i1;
	i__3 = ir - 1;
	inxcb_(&i__2, &i__3, &im3, &nm3);
	i__2 = i__ - i1;
	i__3 = ir - 1;
	inxcb_(&i__2, &i__3, &im1, &nm1);
	(*prdct)(&nm2, &b[im2], &nm3, &b[im3], &nm1, &b[im1], &c__0, &dum, &
		w1[1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[
		1]);
	(*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], &w1[
		1], &w1[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1], &wu[1])
		;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = j + i__ * y_dim1;
	    i__4 = j + i__ * y_dim1;
	    i__5 = j;
	    q__1.r = y[i__4].r - w1[i__5].r, q__1.i = y[i__4].i - w1[i__5].i;
	    y[i__3].r = q__1.r, y[i__3].i = q__1.i;
/* L125: */
	}
/* L126: */
    }

    izr = ccblk_1.nm;
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
	    inxcc_(&i__, &ir, &idxc, &nc);
	    inxcb_(&ipi2, &ir, &ip2, &np2);
	    inxcb_(&ipi1, &irm1, &ip1, &np1);
	    inxcb_(&ipi3, &irm1, &ip3, &np3);
	    (*prdct)(&np2, &b[ip2], &np1, &b[ip1], &np3, &b[ip3], &c__0, &dum,
		     &w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1]
		    , &wu[1]);
	    (*prdct)(&np1, &b[ip1], &c__0, &dum, &c__0, &dum, &nc, &cn[idxc], 
		    &w2[1], &w2[1], m, &am[1], &bm[1], &cm[1], &wd[1], &ww[1],
		     &wu[1]);
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		i__5 = j + i__ * y_dim1;
		i__6 = j + i__ * y_dim1;
		i__7 = j;
		q__1.r = y[i__6].r - w2[i__7].r, q__1.i = y[i__6].i - w2[i__7]
			.i;
		y[i__5].r = q__1.r, y[i__5].i = q__1.i;
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
    i__1 = ccblk_1.k;
    for (ll = 1; ll <= i__1; ++ll) {
	l = ccblk_1.k - ll + 1;
	ir = l - 1;
	irm1 = ir - 1;
	i2 = pow_ii(&c__2, &ir);
	i1 = i2 / 2;
	i4 = i2 + i2;
	ifd = if__ - i2;
	i__3 = ifd;
	i__2 = i4;
	for (i__ = i2; i__2 < 0 ? i__ >= i__3 : i__ <= i__3; i__ += i__2) {
	    if (i__ - ccblk_1.nm <= 0) {
		goto L133;
	    } else {
		goto L143;
	    }
L133:
	    imi1 = i__ - i1;
	    imi2 = i__ - i2;
	    ipi1 = i__ + i1;
	    ipi2 = i__ + i2;
	    inxca_(&i__, &ir, &idxa, &na);
	    inxcc_(&i__, &ir, &idxc, &nc);
	    inxcb_(&i__, &ir, &iz, &nz);
	    inxcb_(&imi1, &irm1, &im1, &nm1);
	    inxcb_(&ipi1, &irm1, &ip1, &np1);
	    if (i__ - i2 <= 0) {
		goto L134;
	    } else {
		goto L136;
	    }
L134:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		i__5 = j;
		w1[i__5].r = 0.f, w1[i__5].i = 0.f;
/* L135: */
	    }
	    goto L137;
L136:
	    (*prdct)(&nm1, &b[im1], &c__0, &dum, &c__0, &dum, &na, &an[idxa], 
		    &y[imi2 * y_dim1 + 1], &w1[1], m, &am[1], &bm[1], &cm[1], 
		    &wd[1], &ww[1], &wu[1]);
L137:
	    if (ipi2 - ccblk_1.nm <= 0) {
		goto L140;
	    } else {
		goto L138;
	    }
L138:
	    i__4 = *m;
	    for (j = 1; j <= i__4; ++j) {
		i__5 = j;
		w2[i__5].r = 0.f, w2[i__5].i = 0.f;
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
		i__5 = j;
		i__6 = j + i__ * y_dim1;
		i__7 = j;
		q__2.r = y[i__6].r + w1[i__7].r, q__2.i = y[i__6].i + w1[i__7]
			.i;
		i__8 = j;
		q__1.r = q__2.r + w2[i__8].r, q__1.i = q__2.i + w2[i__8].i;
		w1[i__5].r = q__1.r, w1[i__5].i = q__1.i;
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
} /* cblkt1_ */

