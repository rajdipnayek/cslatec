/* procp.f -- translated by f2c (version 12.02.01).
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

/* DECK PROCP */
/* Subroutine */ int procp_(integer *nd, real *bd, integer *nm1, real *bm1, 
	integer *nm2, real *bm2, integer *na, real *aa, complex *x, complex *
	y, integer *m, complex *a, complex *b, complex *c__, complex *d__, 
	complex *u, complex *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer j, k;
    static complex v;
    static integer m1, m2;
    static complex bh;
    static integer ia, id;
    static complex am;
    static integer mm;
    static complex ym;
    static real rt;
    static integer mm2;
    static complex den;
    static integer ibr;

/* ***BEGIN PROLOGUE  PROCP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (PRODP-C, PROCP-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* PROCP applies a sequence of matrix operations to the vector X and */
/* stores the result in Y (periodic boundary conditions). */

/* BD,BM1,BM2 are arrays containing roots of certain B polynomials. */
/* ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively. */
/* AA         Array containing scalar multipliers of the vector X. */
/* NA         is the length of the array AA. */
/* X,Y        The matrix operations are applied to X and the result is Y. */
/* A,B,C      are arrays which contain the tridiagonal matrix. */
/* M          is the order of the matrix. */
/* D,U,W      are working arrays. */
/* IS         determines whether or not a change in sign is made. */

/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PROCP */

/* ***FIRST EXECUTABLE STATEMENT  PROCP */
    /* Parameter adjustments */
    --w;
    --u;
    --d__;
    --c__;
    --b;
    --a;
    --y;
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
	y[i__2].r = x[i__3].r, y[i__2].i = x[i__3].i;
	i__2 = j;
	i__3 = j;
	w[i__2].r = y[i__3].r, w[i__2].i = y[i__3].i;
/* L101: */
    }
    mm = *m - 1;
    mm2 = *m - 2;
    id = *nd;
    ibr = 0;
    m1 = *nm1;
    m2 = *nm2;
    ia = *na;
L102:
    if (ia <= 0) {
	goto L105;
    } else {
	goto L103;
    }
L103:
    rt = aa[ia];
    if (*nd == 0) {
	rt = -rt;
    }
    --ia;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	q__1.r = rt * w[i__3].r, q__1.i = rt * w[i__3].i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L104: */
    }
L105:
    if (id <= 0) {
	goto L128;
    } else {
	goto L106;
    }
L106:
    rt = bd[id];
    --id;
    if (id == 0) {
	ibr = 1;
    }

/* BEGIN SOLUTION TO SYSTEM */

    i__1 = *m;
    q__1.r = b[i__1].r - rt, q__1.i = b[i__1].i;
    bh.r = q__1.r, bh.i = q__1.i;
    i__1 = *m;
    ym.r = y[i__1].r, ym.i = y[i__1].i;
    q__1.r = b[1].r - rt, q__1.i = b[1].i;
    den.r = q__1.r, den.i = q__1.i;
    c_div(&q__1, &c__[1], &den);
    d__[1].r = q__1.r, d__[1].i = q__1.i;
    c_div(&q__1, &a[1], &den);
    u[1].r = q__1.r, u[1].i = q__1.i;
    c_div(&q__1, &y[1], &den);
    w[1].r = q__1.r, w[1].i = q__1.i;
    i__1 = *m;
    v.r = c__[i__1].r, v.i = c__[i__1].i;
    if (mm2 - 2 >= 0) {
	goto L107;
    } else {
	goto L109;
    }
L107:
    i__1 = mm2;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	q__2.r = b[i__2].r - rt, q__2.i = b[i__2].i;
	i__3 = j;
	i__4 = j - 1;
	q__3.r = a[i__3].r * d__[i__4].r - a[i__3].i * d__[i__4].i, q__3.i = 
		a[i__3].r * d__[i__4].i + a[i__3].i * d__[i__4].r;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	den.r = q__1.r, den.i = q__1.i;
	i__2 = j;
	c_div(&q__1, &c__[j], &den);
	d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	i__2 = j;
	i__3 = j;
	q__3.r = -a[i__3].r, q__3.i = -a[i__3].i;
	i__4 = j - 1;
	q__2.r = q__3.r * u[i__4].r - q__3.i * u[i__4].i, q__2.i = q__3.r * u[
		i__4].i + q__3.i * u[i__4].r;
	c_div(&q__1, &q__2, &den);
	u[i__2].r = q__1.r, u[i__2].i = q__1.i;
	i__2 = j;
	i__3 = j;
	i__4 = j;
	i__5 = j - 1;
	q__3.r = a[i__4].r * w[i__5].r - a[i__4].i * w[i__5].i, q__3.i = a[
		i__4].r * w[i__5].i + a[i__4].i * w[i__5].r;
	q__2.r = y[i__3].r - q__3.r, q__2.i = y[i__3].i - q__3.i;
	c_div(&q__1, &q__2, &den);
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
	i__2 = j - 1;
	q__2.r = v.r * u[i__2].r - v.i * u[i__2].i, q__2.i = v.r * u[i__2].i 
		+ v.i * u[i__2].r;
	q__1.r = bh.r - q__2.r, q__1.i = bh.i - q__2.i;
	bh.r = q__1.r, bh.i = q__1.i;
	i__2 = j - 1;
	q__2.r = v.r * w[i__2].r - v.i * w[i__2].i, q__2.i = v.r * w[i__2].i 
		+ v.i * w[i__2].r;
	q__1.r = ym.r - q__2.r, q__1.i = ym.i - q__2.i;
	ym.r = q__1.r, ym.i = q__1.i;
	q__2.r = -v.r, q__2.i = -v.i;
	i__2 = j - 1;
	q__1.r = q__2.r * d__[i__2].r - q__2.i * d__[i__2].i, q__1.i = q__2.r 
		* d__[i__2].i + q__2.i * d__[i__2].r;
	v.r = q__1.r, v.i = q__1.i;
/* L108: */
    }
L109:
    i__1 = *m - 1;
    q__2.r = b[i__1].r - rt, q__2.i = b[i__1].i;
    i__2 = *m - 1;
    i__3 = *m - 2;
    q__3.r = a[i__2].r * d__[i__3].r - a[i__2].i * d__[i__3].i, q__3.i = a[
	    i__2].r * d__[i__3].i + a[i__2].i * d__[i__3].r;
    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
    den.r = q__1.r, den.i = q__1.i;
    i__1 = *m - 1;
    i__2 = *m - 1;
    i__3 = *m - 1;
    i__4 = *m - 2;
    q__3.r = a[i__3].r * u[i__4].r - a[i__3].i * u[i__4].i, q__3.i = a[i__3]
	    .r * u[i__4].i + a[i__3].i * u[i__4].r;
    q__2.r = c__[i__2].r - q__3.r, q__2.i = c__[i__2].i - q__3.i;
    c_div(&q__1, &q__2, &den);
    d__[i__1].r = q__1.r, d__[i__1].i = q__1.i;
    i__1 = *m - 1;
    i__2 = *m - 1;
    i__3 = *m - 1;
    i__4 = *m - 2;
    q__3.r = a[i__3].r * w[i__4].r - a[i__3].i * w[i__4].i, q__3.i = a[i__3]
	    .r * w[i__4].i + a[i__3].i * w[i__4].r;
    q__2.r = y[i__2].r - q__3.r, q__2.i = y[i__2].i - q__3.i;
    c_div(&q__1, &q__2, &den);
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    i__1 = *m;
    i__2 = *m - 2;
    q__2.r = v.r * d__[i__2].r - v.i * d__[i__2].i, q__2.i = v.r * d__[i__2]
	    .i + v.i * d__[i__2].r;
    q__1.r = a[i__1].r - q__2.r, q__1.i = a[i__1].i - q__2.i;
    am.r = q__1.r, am.i = q__1.i;
    i__1 = *m - 2;
    q__2.r = v.r * u[i__1].r - v.i * u[i__1].i, q__2.i = v.r * u[i__1].i + 
	    v.i * u[i__1].r;
    q__1.r = bh.r - q__2.r, q__1.i = bh.i - q__2.i;
    bh.r = q__1.r, bh.i = q__1.i;
    i__1 = *m - 2;
    q__2.r = v.r * w[i__1].r - v.i * w[i__1].i, q__2.i = v.r * w[i__1].i + 
	    v.i * w[i__1].r;
    q__1.r = ym.r - q__2.r, q__1.i = ym.i - q__2.i;
    ym.r = q__1.r, ym.i = q__1.i;
    i__1 = *m - 1;
    q__2.r = am.r * d__[i__1].r - am.i * d__[i__1].i, q__2.i = am.r * d__[
	    i__1].i + am.i * d__[i__1].r;
    q__1.r = bh.r - q__2.r, q__1.i = bh.i - q__2.i;
    den.r = q__1.r, den.i = q__1.i;
    if (c_abs(&den) != 0.f) {
	goto L110;
    } else {
	goto L111;
    }
L110:
    i__1 = *m;
    i__2 = *m - 1;
    q__3.r = am.r * w[i__2].r - am.i * w[i__2].i, q__3.i = am.r * w[i__2].i + 
	    am.i * w[i__2].r;
    q__2.r = ym.r - q__3.r, q__2.i = ym.i - q__3.i;
    c_div(&q__1, &q__2, &den);
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    goto L112;
L111:
    i__1 = *m;
    w[i__1].r = 1.f, w[i__1].i = 0.f;
L112:
    i__1 = *m - 1;
    i__2 = *m - 1;
    i__3 = *m - 1;
    i__4 = *m;
    q__2.r = d__[i__3].r * w[i__4].r - d__[i__3].i * w[i__4].i, q__2.i = d__[
	    i__3].r * w[i__4].i + d__[i__3].i * w[i__4].r;
    q__1.r = w[i__2].r - q__2.r, q__1.i = w[i__2].i - q__2.i;
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	k = *m - j;
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k + 1;
	q__3.r = d__[i__4].r * w[i__5].r - d__[i__4].i * w[i__5].i, q__3.i = 
		d__[i__4].r * w[i__5].i + d__[i__4].i * w[i__5].r;
	q__2.r = w[i__3].r - q__3.r, q__2.i = w[i__3].i - q__3.i;
	i__6 = k;
	i__7 = *m;
	q__4.r = u[i__6].r * w[i__7].r - u[i__6].i * w[i__7].i, q__4.i = u[
		i__6].r * w[i__7].i + u[i__6].i * w[i__7].r;
	q__1.r = q__2.r - q__4.r, q__1.i = q__2.i - q__4.i;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
/* L113: */
    }
    if (*na <= 0) {
	goto L116;
    } else {
	goto L102;
    }
L114:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	y[i__2].r = w[i__3].r, y[i__2].i = w[i__3].i;
/* L115: */
    }
    ibr = 1;
    goto L102;
L116:
    if (m1 <= 0) {
	goto L117;
    } else {
	goto L118;
    }
L117:
    if (m2 <= 0) {
	goto L114;
    } else {
	goto L123;
    }
L118:
    if (m2 <= 0) {
	goto L120;
    } else {
	goto L119;
    }
L119:
    if ((r__1 = bm1[m1], dabs(r__1)) - (r__2 = bm2[m2], dabs(r__2)) <= 0.f) {
	goto L123;
    } else {
	goto L120;
    }
L120:
    if (ibr <= 0) {
	goto L121;
    } else {
	goto L122;
    }
L121:
    if ((r__1 = bm1[m1] - bd[id], dabs(r__1)) - (r__2 = bm1[m1] - rt, dabs(
	    r__2)) >= 0.f) {
	goto L122;
    } else {
	goto L114;
    }
L122:
    rt -= bm1[m1];
    --m1;
    goto L126;
L123:
    if (ibr <= 0) {
	goto L124;
    } else {
	goto L125;
    }
L124:
    if ((r__1 = bm2[m2] - bd[id], dabs(r__1)) - (r__2 = bm2[m2] - rt, dabs(
	    r__2)) >= 0.f) {
	goto L125;
    } else {
	goto L114;
    }
L125:
    rt -= bm2[m2];
    --m2;
L126:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	i__4 = j;
	q__2.r = rt * w[i__4].r, q__2.i = rt * w[i__4].i;
	q__1.r = y[i__3].r + q__2.r, q__1.i = y[i__3].i + q__2.i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L127: */
    }
    goto L102;
L128:
    return 0;
} /* procp_ */

