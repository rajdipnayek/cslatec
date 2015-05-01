/* cpadd.f -- translated by f2c (version 12.02.01).
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

static complex c_b24 = {1.f,0.f};
static integer c__2 = 2;

/* DECK CPADD */
/* Subroutine */ int cpadd_(integer *n, integer *ierror, real *a, real *c__, 
	complex *cbp, real *bp, real *bh)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Local variables */
    static complex f;
    static integer j, i3;
    static complex r1, r2, r3;
    static real db;
    static complex dd;
    static integer if__, ig;
    static complex fp, cx;
    static integer is, it, nt, iz;
    static real xl, xm, xr;
    static complex fsg, hsg;
    static integer icv;
    static complex fpp;
    static real sgn, psg;
    extern doublereal bcrh_(real *, real *, integer *, real *, real *, real *,
	     E_fp, real *);
    static complex cdis;
    extern doublereal pgsf_(real *, integer *, real *, real *, real *);
    static real scnv;
    static integer nhalf;
    extern doublereal ppgsf_(real *, integer *, real *, real *, real *);
    static integer modiz;
    extern /* Subroutine */ int pppsf_();

/* ***BEGIN PROLOGUE  CPADD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CPADD-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*   CPADD computes the eigenvalues of the periodic tridiagonal matrix */
/*   with coefficients AN,BN,CN. */

/*   N    is the order of the BH and BP polynomials. */
/*   BP   contains the eigenvalues on output. */
/*   CBP  is the same as BP except type complex. */
/*   BH   is used to temporarily store the roots of the B HAT polynomial */
/*        which enters through BP. */

/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  BCRH, PGSF, PPGSF, PPPSF */
/* ***COMMON BLOCKS    CCBLK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CPADD */

/* ***FIRST EXECUTABLE STATEMENT  CPADD */
    /* Parameter adjustments */
    --bh;
    --bp;
    --cbp;
    --c__;
    --a;

    /* Function Body */
    scnv = sqrt(ccblk_1.cnv);
    iz = *n;
    if ((r__1 = bp[*n] - bp[1]) < 0.f) {
	goto L101;
    } else if (r__1 == 0) {
	goto L142;
    } else {
	goto L103;
    }
L101:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	nt = *n - j;
	bh[j] = bp[nt + 1];
/* L102: */
    }
    goto L105;
L103:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	bh[j] = bp[j];
/* L104: */
    }
L105:
    ccblk_1.ncmplx = 0;
    modiz = iz % 2;
    is = 1;
    if (modiz != 0) {
	goto L106;
    } else {
	goto L107;
    }
L106:
    if (a[1] < 0.f) {
	goto L110;
    } else if (a[1] == 0) {
	goto L142;
    } else {
	goto L107;
    }
L107:
    xl = bh[1];
    db = bh[3] - bh[1];
L108:
    xl -= db;
    if (pgsf_(&xl, &iz, &c__[1], &a[1], &bh[1]) <= 0.f) {
	goto L108;
    } else {
	goto L109;
    }
L109:
    sgn = -1.f;
    r__1 = bcrh_(&xl, &bh[1], &iz, &c__[1], &a[1], &bh[1], (E_fp)pgsf_, &sgn);
    q__1.r = r__1, q__1.i = 0.f;
    cbp[1].r = q__1.r, cbp[1].i = q__1.i;
    is = 2;
L110:
    if__ = iz - 1;
    if (modiz != 0) {
	goto L111;
    } else {
	goto L112;
    }
L111:
    if (a[1] < 0.f) {
	goto L112;
    } else if (a[1] == 0) {
	goto L142;
    } else {
	goto L115;
    }
L112:
    xr = bh[iz];
    db = bh[iz] - bh[iz - 2];
L113:
    xr += db;
    if (pgsf_(&xr, &iz, &c__[1], &a[1], &bh[1]) >= 0.f) {
	goto L114;
    } else {
	goto L113;
    }
L114:
    sgn = 1.f;
    i__1 = iz;
    r__1 = bcrh_(&bh[iz], &xr, &iz, &c__[1], &a[1], &bh[1], (E_fp)pgsf_, &sgn)
	    ;
    q__1.r = r__1, q__1.i = 0.f;
    cbp[i__1].r = q__1.r, cbp[i__1].i = q__1.i;
    if__ = iz - 2;
L115:
    i__1 = if__;
    for (ig = is; ig <= i__1; ig += 2) {
	xl = bh[ig];
	xr = bh[ig + 1];
	sgn = -1.f;
	xm = bcrh_(&xl, &xr, &iz, &c__[1], &a[1], &bh[1], (E_fp)pppsf_, &sgn);
	psg = pgsf_(&xm, &iz, &c__[1], &a[1], &bh[1]);
	if (dabs(psg) - ccblk_1.eps <= 0.f) {
	    goto L118;
	} else {
	    goto L116;
	}
L116:
	if ((r__1 = psg * ppgsf_(&xm, &iz, &c__[1], &a[1], &bh[1])) < 0.f) {
	    goto L117;
	} else if (r__1 == 0) {
	    goto L118;
	} else {
	    goto L119;
	}

/*     CASE OF A REAL ZERO */

L117:
	sgn = 1.f;
	i__2 = ig;
	r__1 = bcrh_(&bh[ig], &xm, &iz, &c__[1], &a[1], &bh[1], (E_fp)pgsf_, &
		sgn);
	q__1.r = r__1, q__1.i = 0.f;
	cbp[i__2].r = q__1.r, cbp[i__2].i = q__1.i;
	sgn = -1.f;
	i__2 = ig + 1;
	r__1 = bcrh_(&xm, &bh[ig + 1], &iz, &c__[1], &a[1], &bh[1], (E_fp)
		pgsf_, &sgn);
	q__1.r = r__1, q__1.i = 0.f;
	cbp[i__2].r = q__1.r, cbp[i__2].i = q__1.i;
	goto L136;

/*     CASE OF A MULTIPLE ZERO */

L118:
	i__2 = ig;
	q__1.r = xm, q__1.i = 0.f;
	cbp[i__2].r = q__1.r, cbp[i__2].i = q__1.i;
	i__2 = ig + 1;
	q__1.r = xm, q__1.i = 0.f;
	cbp[i__2].r = q__1.r, cbp[i__2].i = q__1.i;
	goto L136;

/*     CASE OF A COMPLEX ZERO */

L119:
	it = 0;
	icv = 0;
	q__1.r = xm, q__1.i = 0.f;
	cx.r = q__1.r, cx.i = q__1.i;
L120:
	fsg.r = 1.f, fsg.i = 0.f;
	hsg.r = 1.f, hsg.i = 0.f;
	fp.r = 0.f, fp.i = 0.f;
	fpp.r = 0.f, fpp.i = 0.f;
	i__2 = iz;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = j;
	    q__2.r = cx.r - bh[i__3], q__2.i = cx.i;
	    c_div(&q__1, &c_b24, &q__2);
	    dd.r = q__1.r, dd.i = q__1.i;
	    i__3 = j;
	    q__2.r = a[i__3] * fsg.r, q__2.i = a[i__3] * fsg.i;
	    q__1.r = q__2.r * dd.r - q__2.i * dd.i, q__1.i = q__2.r * dd.i + 
		    q__2.i * dd.r;
	    fsg.r = q__1.r, fsg.i = q__1.i;
	    i__3 = j;
	    q__2.r = c__[i__3] * hsg.r, q__2.i = c__[i__3] * hsg.i;
	    q__1.r = q__2.r * dd.r - q__2.i * dd.i, q__1.i = q__2.r * dd.i + 
		    q__2.i * dd.r;
	    hsg.r = q__1.r, hsg.i = q__1.i;
	    q__1.r = fp.r + dd.r, q__1.i = fp.i + dd.i;
	    fp.r = q__1.r, fp.i = q__1.i;
	    q__2.r = dd.r * dd.r - dd.i * dd.i, q__2.i = dd.r * dd.i + dd.i * 
		    dd.r;
	    q__1.r = fpp.r - q__2.r, q__1.i = fpp.i - q__2.i;
	    fpp.r = q__1.r, fpp.i = q__1.i;
/* L121: */
	}
	if (modiz != 0) {
	    goto L123;
	} else {
	    goto L122;
	}
L122:
	q__2.r = 1.f - fsg.r, q__2.i = 0.f - fsg.i;
	q__1.r = q__2.r - hsg.r, q__1.i = q__2.i - hsg.i;
	f.r = q__1.r, f.i = q__1.i;
	goto L124;
L123:
	q__2.r = fsg.r + 1.f, q__2.i = fsg.i + 0.f;
	q__1.r = q__2.r + hsg.r, q__1.i = q__2.i + hsg.i;
	f.r = q__1.r, f.i = q__1.i;
L124:
	i3 = 0;
	if (c_abs(&fp) <= 0.f) {
	    goto L126;
	} else {
	    goto L125;
	}
L125:
	i3 = 1;
	q__2.r = -f.r, q__2.i = -f.i;
	c_div(&q__1, &q__2, &fp);
	r3.r = q__1.r, r3.i = q__1.i;
L126:
	if (c_abs(&fpp) <= 0.f) {
	    goto L132;
	} else {
	    goto L127;
	}
L127:
	pow_ci(&q__3, &fp, &c__2);
	q__5.r = f.r * 2.f, q__5.i = f.i * 2.f;
	q__4.r = q__5.r * fpp.r - q__5.i * fpp.i, q__4.i = q__5.r * fpp.i + 
		q__5.i * fpp.r;
	q__2.r = q__3.r - q__4.r, q__2.i = q__3.i - q__4.i;
	c_sqrt(&q__1, &q__2);
	cdis.r = q__1.r, cdis.i = q__1.i;
	q__1.r = cdis.r - fp.r, q__1.i = cdis.i - fp.i;
	r1.r = q__1.r, r1.i = q__1.i;
	q__2.r = -fp.r, q__2.i = -fp.i;
	q__1.r = q__2.r - cdis.r, q__1.i = q__2.i - cdis.i;
	r2.r = q__1.r, r2.i = q__1.i;
	if (c_abs(&r1) - c_abs(&r2) <= 0.f) {
	    goto L129;
	} else {
	    goto L128;
	}
L128:
	c_div(&q__1, &r1, &fpp);
	r1.r = q__1.r, r1.i = q__1.i;
	goto L130;
L129:
	c_div(&q__1, &r2, &fpp);
	r1.r = q__1.r, r1.i = q__1.i;
L130:
	q__3.r = f.r * 2.f, q__3.i = f.i * 2.f;
	c_div(&q__2, &q__3, &fpp);
	c_div(&q__1, &q__2, &r1);
	r2.r = q__1.r, r2.i = q__1.i;
	if (c_abs(&r2) < c_abs(&r1)) {
	    r1.r = r2.r, r1.i = r2.i;
	}
	if (i3 <= 0) {
	    goto L133;
	} else {
	    goto L131;
	}
L131:
	if (c_abs(&r3) < c_abs(&r1)) {
	    r1.r = r3.r, r1.i = r3.i;
	}
	goto L133;
L132:
	r1.r = r3.r, r1.i = r3.i;
L133:
	q__1.r = cx.r + r1.r, q__1.i = cx.i + r1.i;
	cx.r = q__1.r, cx.i = q__1.i;
	++it;
	if (it > 50) {
	    goto L142;
	}
	if (c_abs(&r1) > scnv) {
	    goto L120;
	}
	if (icv <= 0) {
	    goto L134;
	} else {
	    goto L135;
	}
L134:
	icv = 1;
	goto L120;
L135:
	i__2 = ig;
	cbp[i__2].r = cx.r, cbp[i__2].i = cx.i;
	i__2 = ig + 1;
	r_cnjg(&q__1, &cx);
	cbp[i__2].r = q__1.r, cbp[i__2].i = q__1.i;
L136:
	;
    }
    if ((r__1 = c_abs(&cbp[*n]) - c_abs(&cbp[1])) < 0.f) {
	goto L137;
    } else if (r__1 == 0) {
	goto L142;
    } else {
	goto L139;
    }
L137:
    nhalf = *n / 2;
    i__1 = nhalf;
    for (j = 1; j <= i__1; ++j) {
	nt = *n - j;
	i__2 = j;
	cx.r = cbp[i__2].r, cx.i = cbp[i__2].i;
	i__2 = j;
	i__3 = nt + 1;
	cbp[i__2].r = cbp[i__3].r, cbp[i__2].i = cbp[i__3].i;
	i__2 = nt + 1;
	cbp[i__2].r = cx.r, cbp[i__2].i = cx.i;
/* L138: */
    }
L139:
    ccblk_1.ncmplx = 1;
    i__1 = iz;
    for (j = 2; j <= i__1; ++j) {
	if (r_imag(&cbp[j]) != 0.f) {
	    goto L143;
	} else {
	    goto L140;
	}
L140:
	;
    }
    ccblk_1.ncmplx = 0;
    i__1 = iz;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	bp[j] = cbp[i__2].r;
/* L141: */
    }
    goto L143;
L142:
    *ierror = 4;
L143:
    return 0;
} /* cpadd_ */

