/* cprodp.f -- translated by f2c (version 12.02.01).
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

/* DECK CPRODP */
/* Subroutine */ int cprodp_(integer *nd, complex *bd, integer *nm1, real *
	bm1, integer *nm2, real *bm2, integer *na, real *aa, real *x, real *
	yy, integer *m, real *a, real *b, real *c__, complex *d__, complex *u,
	 complex *y)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Local variables */
    static integer j, k;
    static complex v;
    static integer m1, m2;
    static complex y1, y2, bh;
    static integer ia, id;
    static complex am;
    static integer mm;
    static complex yh, ym;
    static real rt;
    static integer mm2;
    static complex den, crt;
    static integer iflg;

/* ***BEGIN PROLOGUE  CPRODP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BLKTRI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CPRODP-S, CPROCP-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* PRODP applies a sequence of matrix operations to the vector X and */
/* stores the result in YY. (Periodic boundary conditions and COMPLEX */
/* case) */

/* BD,BM1,BM2     are arrays containing roots of certain B polynomials. */
/* ND,NM1,NM2     are the lengths of the arrays BD,BM1,BM2 respectively. */
/* AA             Array containing scalar multipliers of the vector X. */
/* NA             is the length of the array AA. */
/* X,YY      The matrix operations are applied to X and the result is YY. */
/* A,B,C          are arrays which contain the tridiagonal matrix. */
/* M              is the order of the matrix. */
/* D,U,Y          are working arrays. */
/* ISGN           determines whether or not a change in sign is made. */

/* ***SEE ALSO  BLKTRI */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CPRODP */

/* ***FIRST EXECUTABLE STATEMENT  CPRODP */
    /* Parameter adjustments */
    --y;
    --u;
    --d__;
    --c__;
    --b;
    --a;
    --yy;
    --x;
    --aa;
    --bm2;
    --bm1;
    --bd;

    /* Function Body */
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	q__1.r = x[i__3], q__1.i = 0.f;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L101: */
    }
    mm = *m - 1;
    mm2 = *m - 2;
    id = *nd;
    m1 = *nm1;
    m2 = *nm2;
    ia = *na;
L102:
    iflg = 0;
    if (id <= 0) {
	goto L111;
    } else {
	goto L103;
    }
L103:
    i__1 = id;
    crt.r = bd[i__1].r, crt.i = bd[i__1].i;
    --id;
    iflg = 1;

/* BEGIN SOLUTION TO SYSTEM */

    i__1 = *m;
    q__1.r = b[i__1] - crt.r, q__1.i = -crt.i;
    bh.r = q__1.r, bh.i = q__1.i;
    i__1 = *m;
    ym.r = y[i__1].r, ym.i = y[i__1].i;
    q__1.r = b[1] - crt.r, q__1.i = -crt.i;
    den.r = q__1.r, den.i = q__1.i;
    q__2.r = c__[1], q__2.i = 0.f;
    c_div(&q__1, &q__2, &den);
    d__[1].r = q__1.r, d__[1].i = q__1.i;
    q__2.r = a[1], q__2.i = 0.f;
    c_div(&q__1, &q__2, &den);
    u[1].r = q__1.r, u[1].i = q__1.i;
    c_div(&q__1, &y[1], &den);
    y[1].r = q__1.r, y[1].i = q__1.i;
    i__1 = *m;
    q__1.r = c__[i__1], q__1.i = 0.f;
    v.r = q__1.r, v.i = q__1.i;
    if (mm2 - 2 >= 0) {
	goto L104;
    } else {
	goto L106;
    }
L104:
    i__1 = mm2;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	q__2.r = b[i__2] - crt.r, q__2.i = -crt.i;
	i__3 = j;
	i__4 = j - 1;
	q__3.r = a[i__3] * d__[i__4].r, q__3.i = a[i__3] * d__[i__4].i;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	den.r = q__1.r, den.i = q__1.i;
	i__2 = j;
	i__3 = j;
	q__2.r = c__[i__3], q__2.i = 0.f;
	c_div(&q__1, &q__2, &den);
	d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	i__2 = j;
	r__1 = -a[j];
	i__3 = j - 1;
	q__2.r = r__1 * u[i__3].r, q__2.i = r__1 * u[i__3].i;
	c_div(&q__1, &q__2, &den);
	u[i__2].r = q__1.r, u[i__2].i = q__1.i;
	i__2 = j;
	i__3 = j;
	i__4 = j;
	i__5 = j - 1;
	q__3.r = a[i__4] * y[i__5].r, q__3.i = a[i__4] * y[i__5].i;
	q__2.r = y[i__3].r - q__3.r, q__2.i = y[i__3].i - q__3.i;
	c_div(&q__1, &q__2, &den);
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
	i__2 = j - 1;
	q__2.r = v.r * u[i__2].r - v.i * u[i__2].i, q__2.i = v.r * u[i__2].i 
		+ v.i * u[i__2].r;
	q__1.r = bh.r - q__2.r, q__1.i = bh.i - q__2.i;
	bh.r = q__1.r, bh.i = q__1.i;
	i__2 = j - 1;
	q__2.r = v.r * y[i__2].r - v.i * y[i__2].i, q__2.i = v.r * y[i__2].i 
		+ v.i * y[i__2].r;
	q__1.r = ym.r - q__2.r, q__1.i = ym.i - q__2.i;
	ym.r = q__1.r, ym.i = q__1.i;
	q__2.r = -v.r, q__2.i = -v.i;
	i__2 = j - 1;
	q__1.r = q__2.r * d__[i__2].r - q__2.i * d__[i__2].i, q__1.i = q__2.r 
		* d__[i__2].i + q__2.i * d__[i__2].r;
	v.r = q__1.r, v.i = q__1.i;
/* L105: */
    }
L106:
    i__1 = *m - 1;
    q__2.r = b[i__1] - crt.r, q__2.i = -crt.i;
    i__2 = *m - 1;
    i__3 = *m - 2;
    q__3.r = a[i__2] * d__[i__3].r, q__3.i = a[i__2] * d__[i__3].i;
    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
    den.r = q__1.r, den.i = q__1.i;
    i__1 = *m - 1;
    i__2 = *m - 1;
    i__3 = *m - 1;
    i__4 = *m - 2;
    q__3.r = a[i__3] * u[i__4].r, q__3.i = a[i__3] * u[i__4].i;
    q__2.r = c__[i__2] - q__3.r, q__2.i = -q__3.i;
    c_div(&q__1, &q__2, &den);
    d__[i__1].r = q__1.r, d__[i__1].i = q__1.i;
    i__1 = *m - 1;
    i__2 = *m - 1;
    i__3 = *m - 1;
    i__4 = *m - 2;
    q__3.r = a[i__3] * y[i__4].r, q__3.i = a[i__3] * y[i__4].i;
    q__2.r = y[i__2].r - q__3.r, q__2.i = y[i__2].i - q__3.i;
    c_div(&q__1, &q__2, &den);
    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
    i__1 = *m;
    i__2 = *m - 2;
    q__2.r = v.r * d__[i__2].r - v.i * d__[i__2].i, q__2.i = v.r * d__[i__2]
	    .i + v.i * d__[i__2].r;
    q__1.r = a[i__1] - q__2.r, q__1.i = -q__2.i;
    am.r = q__1.r, am.i = q__1.i;
    i__1 = *m - 2;
    q__2.r = v.r * u[i__1].r - v.i * u[i__1].i, q__2.i = v.r * u[i__1].i + 
	    v.i * u[i__1].r;
    q__1.r = bh.r - q__2.r, q__1.i = bh.i - q__2.i;
    bh.r = q__1.r, bh.i = q__1.i;
    i__1 = *m - 2;
    q__2.r = v.r * y[i__1].r - v.i * y[i__1].i, q__2.i = v.r * y[i__1].i + 
	    v.i * y[i__1].r;
    q__1.r = ym.r - q__2.r, q__1.i = ym.i - q__2.i;
    ym.r = q__1.r, ym.i = q__1.i;
    i__1 = *m - 1;
    q__2.r = am.r * d__[i__1].r - am.i * d__[i__1].i, q__2.i = am.r * d__[
	    i__1].i + am.i * d__[i__1].r;
    q__1.r = bh.r - q__2.r, q__1.i = bh.i - q__2.i;
    den.r = q__1.r, den.i = q__1.i;
    if (c_abs(&den) != 0.f) {
	goto L107;
    } else {
	goto L108;
    }
L107:
    i__1 = *m;
    i__2 = *m - 1;
    q__3.r = am.r * y[i__2].r - am.i * y[i__2].i, q__3.i = am.r * y[i__2].i + 
	    am.i * y[i__2].r;
    q__2.r = ym.r - q__3.r, q__2.i = ym.i - q__3.i;
    c_div(&q__1, &q__2, &den);
    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
    goto L109;
L108:
    i__1 = *m;
    y[i__1].r = 1.f, y[i__1].i = 0.f;
L109:
    i__1 = *m - 1;
    i__2 = *m - 1;
    i__3 = *m - 1;
    i__4 = *m;
    q__2.r = d__[i__3].r * y[i__4].r - d__[i__3].i * y[i__4].i, q__2.i = d__[
	    i__3].r * y[i__4].i + d__[i__3].i * y[i__4].r;
    q__1.r = y[i__2].r - q__2.r, q__1.i = y[i__2].i - q__2.i;
    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	k = *m - j;
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k + 1;
	q__3.r = d__[i__4].r * y[i__5].r - d__[i__4].i * y[i__5].i, q__3.i = 
		d__[i__4].r * y[i__5].i + d__[i__4].i * y[i__5].r;
	q__2.r = y[i__3].r - q__3.r, q__2.i = y[i__3].i - q__3.i;
	i__6 = k;
	i__7 = *m;
	q__4.r = u[i__6].r * y[i__7].r - u[i__6].i * y[i__7].i, q__4.i = u[
		i__6].r * y[i__7].i + u[i__6].i * y[i__7].r;
	q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L110: */
    }
L111:
    if (m1 <= 0) {
	goto L112;
    } else {
	goto L114;
    }
L112:
    if (m2 <= 0) {
	goto L123;
    } else {
	goto L113;
    }
L113:
    rt = bm2[m2];
    --m2;
    goto L119;
L114:
    if (m2 <= 0) {
	goto L115;
    } else {
	goto L116;
    }
L115:
    rt = bm1[m1];
    --m1;
    goto L119;
L116:
    if ((r__1 = bm1[m1], dabs(r__1)) - (r__2 = bm2[m2], dabs(r__2)) <= 0.f) {
	goto L118;
    } else {
	goto L117;
    }
L117:
    rt = bm1[m1];
    --m1;
    goto L119;
L118:
    rt = bm2[m2];
    --m2;

/* MATRIX MULTIPLICATION */

L119:
    yh.r = y[1].r, yh.i = y[1].i;
    r__1 = b[1] - rt;
    q__3.r = r__1 * y[1].r, q__3.i = r__1 * y[1].i;
    q__4.r = c__[1] * y[2].r, q__4.i = c__[1] * y[2].i;
    q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
    i__1 = *m;
    q__5.r = a[1] * y[i__1].r, q__5.i = a[1] * y[i__1].i;
    q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
    y1.r = q__1.r, y1.i = q__1.i;
    if (mm - 2 >= 0) {
	goto L120;
    } else {
	goto L122;
    }
L120:
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j - 1;
	q__3.r = a[i__2] * y[i__3].r, q__3.i = a[i__2] * y[i__3].i;
	r__1 = b[j] - rt;
	i__4 = j;
	q__4.r = r__1 * y[i__4].r, q__4.i = r__1 * y[i__4].i;
	q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
	i__5 = j;
	i__6 = j + 1;
	q__5.r = c__[i__5] * y[i__6].r, q__5.i = c__[i__5] * y[i__6].i;
	q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
	y2.r = q__1.r, y2.i = q__1.i;
	i__2 = j - 1;
	y[i__2].r = y1.r, y[i__2].i = y1.i;
	y1.r = y2.r, y1.i = y2.i;
/* L121: */
    }
L122:
    i__1 = *m;
    i__2 = *m;
    i__3 = *m - 1;
    q__3.r = a[i__2] * y[i__3].r, q__3.i = a[i__2] * y[i__3].i;
    r__1 = b[*m] - rt;
    i__4 = *m;
    q__4.r = r__1 * y[i__4].r, q__4.i = r__1 * y[i__4].i;
    q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
    i__5 = *m;
    q__5.r = c__[i__5] * yh.r, q__5.i = c__[i__5] * yh.i;
    q__1.r = q__2.r + q__5.r, q__1.i = q__2.i + q__5.i;
    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
    i__1 = *m - 1;
    y[i__1].r = y1.r, y[i__1].i = y1.i;
    iflg = 1;
    goto L102;
L123:
    if (ia <= 0) {
	goto L126;
    } else {
	goto L124;
    }
L124:
    rt = aa[ia];
    --ia;
    iflg = 1;

/* SCALAR MULTIPLICATION */

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	q__1.r = rt * y[i__3].r, q__1.i = rt * y[i__3].i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L125: */
    }
L126:
    if (iflg <= 0) {
	goto L127;
    } else {
	goto L102;
    }
L127:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	yy[j] = y[i__2].r;
/* L128: */
    }
    return 0;
} /* cprodp_ */

