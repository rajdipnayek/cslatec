/* cproc.f -- translated by f2c (version 12.02.01).
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

/* DECK CPROC */
/* Subroutine */ int cproc_(integer *nd, complex *bd, integer *nm1, real *bm1,
	 integer *nm2, real *bm2, integer *na, real *aa, complex *x, complex *
	y, integer *m, complex *a, complex *b, complex *c__, complex *d__, 
	complex *w, real *yy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4, q__5, q__6;

    /* Local variables */
    static integer j, k, m1, m2;
    static complex y1, y2;
    static integer ia, id, mm;
    static real rt;
    static complex den, crt;
    static integer iflg;

/* ***BEGIN PROLOGUE  CPROC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (CPROD-S, CPROC-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* PROC applies a sequence of matrix operations to the vector X and */
/* stores the result in Y. */
/* AA     Array containing scalar multipliers of the vector X. */
/* ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively. */
/* BD,BM1,BM2 are arrays containing roots of certain B polynomials. */
/* NA     is the length of the array AA. */
/* X,Y    The matrix operations are applied to X and the result is Y. */
/* A,B,C  are arrays which contain the tridiagonal matrix. */
/* M      is the order of the matrix. */
/* D,W    are work arrays. */
/* ISGN   determines whether or not a change in sign is made. */

/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CPROC */

/* ***FIRST EXECUTABLE STATEMENT  CPROC */
    /* Parameter adjustments */
    --yy;
    --w;
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
/* L101: */
    }
    mm = *m - 1;
    id = *nd;
    m1 = *nm1;
    m2 = *nm2;
    ia = *na;
L102:
    iflg = 0;
    if (id <= 0) {
	goto L109;
    } else {
	goto L103;
    }
L103:
    i__1 = id;
    crt.r = bd[i__1].r, crt.i = bd[i__1].i;
    --id;

/* BEGIN SOLUTION TO SYSTEM */

    i__1 = *m;
    i__2 = *m;
    q__2.r = b[i__2].r - crt.r, q__2.i = b[i__2].i - crt.i;
    c_div(&q__1, &a[*m], &q__2);
    d__[i__1].r = q__1.r, d__[i__1].i = q__1.i;
    i__1 = *m;
    i__2 = *m;
    q__2.r = b[i__2].r - crt.r, q__2.i = b[i__2].i - crt.i;
    c_div(&q__1, &y[*m], &q__2);
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	k = *m - j;
	i__2 = k + 1;
	q__2.r = b[i__2].r - crt.r, q__2.i = b[i__2].i - crt.i;
	i__3 = k + 1;
	i__4 = k + 2;
	q__3.r = c__[i__3].r * d__[i__4].r - c__[i__3].i * d__[i__4].i, 
		q__3.i = c__[i__3].r * d__[i__4].i + c__[i__3].i * d__[i__4]
		.r;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	den.r = q__1.r, den.i = q__1.i;
	i__2 = k + 1;
	c_div(&q__1, &a[k + 1], &den);
	d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	i__2 = k + 1;
	i__3 = k + 1;
	i__4 = k + 1;
	i__5 = k + 2;
	q__3.r = c__[i__4].r * w[i__5].r - c__[i__4].i * w[i__5].i, q__3.i = 
		c__[i__4].r * w[i__5].i + c__[i__4].i * w[i__5].r;
	q__2.r = y[i__3].r - q__3.r, q__2.i = y[i__3].i - q__3.i;
	c_div(&q__1, &q__2, &den);
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
/* L104: */
    }
    q__2.r = b[1].r - crt.r, q__2.i = b[1].i - crt.i;
    q__3.r = c__[1].r * d__[2].r - c__[1].i * d__[2].i, q__3.i = c__[1].r * 
	    d__[2].i + c__[1].i * d__[2].r;
    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
    den.r = q__1.r, den.i = q__1.i;
    if (c_abs(&den) != 0.f) {
	goto L105;
    } else {
	goto L106;
    }
L105:
    q__3.r = c__[1].r * w[2].r - c__[1].i * w[2].i, q__3.i = c__[1].r * w[2]
	    .i + c__[1].i * w[2].r;
    q__2.r = y[1].r - q__3.r, q__2.i = y[1].i - q__3.i;
    c_div(&q__1, &q__2, &den);
    y[1].r = q__1.r, y[1].i = q__1.i;
    goto L107;
L106:
    y[1].r = 1.f, y[1].i = 0.f;
L107:
    i__1 = *m;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	i__4 = j;
	i__5 = j - 1;
	q__2.r = d__[i__4].r * y[i__5].r - d__[i__4].i * y[i__5].i, q__2.i = 
		d__[i__4].r * y[i__5].i + d__[i__4].i * y[i__5].r;
	q__1.r = w[i__3].r - q__2.r, q__1.i = w[i__3].i - q__2.i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L108: */
    }
L109:
    if (m1 <= 0) {
	goto L110;
    } else {
	goto L112;
    }
L110:
    if (m2 <= 0) {
	goto L121;
    } else {
	goto L111;
    }
L111:
    rt = bm2[m2];
    --m2;
    goto L117;
L112:
    if (m2 <= 0) {
	goto L113;
    } else {
	goto L114;
    }
L113:
    rt = bm1[m1];
    --m1;
    goto L117;
L114:
    if ((r__1 = bm1[m1], dabs(r__1)) - (r__2 = bm2[m2], dabs(r__2)) <= 0.f) {
	goto L116;
    } else {
	goto L115;
    }
L115:
    rt = bm1[m1];
    --m1;
    goto L117;
L116:
    rt = bm2[m2];
    --m2;
L117:
    q__3.r = b[1].r - rt, q__3.i = b[1].i;
    q__2.r = q__3.r * y[1].r - q__3.i * y[1].i, q__2.i = q__3.r * y[1].i + 
	    q__3.i * y[1].r;
    q__4.r = c__[1].r * y[2].r - c__[1].i * y[2].i, q__4.i = c__[1].r * y[2]
	    .i + c__[1].i * y[2].r;
    q__1.r = q__2.r + q__4.r, q__1.i = q__2.i + q__4.i;
    y1.r = q__1.r, y1.i = q__1.i;
    if (mm - 2 >= 0) {
	goto L118;
    } else {
	goto L120;
    }

/* MATRIX MULTIPLICATION */

L118:
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j - 1;
	q__3.r = a[i__2].r * y[i__3].r - a[i__2].i * y[i__3].i, q__3.i = a[
		i__2].r * y[i__3].i + a[i__2].i * y[i__3].r;
	i__4 = j;
	q__5.r = b[i__4].r - rt, q__5.i = b[i__4].i;
	i__5 = j;
	q__4.r = q__5.r * y[i__5].r - q__5.i * y[i__5].i, q__4.i = q__5.r * y[
		i__5].i + q__5.i * y[i__5].r;
	q__2.r = q__3.r + q__4.r, q__2.i = q__3.i + q__4.i;
	i__6 = j;
	i__7 = j + 1;
	q__6.r = c__[i__6].r * y[i__7].r - c__[i__6].i * y[i__7].i, q__6.i = 
		c__[i__6].r * y[i__7].i + c__[i__6].i * y[i__7].r;
	q__1.r = q__2.r + q__6.r, q__1.i = q__2.i + q__6.i;
	y2.r = q__1.r, y2.i = q__1.i;
	i__2 = j - 1;
	y[i__2].r = y1.r, y[i__2].i = y1.i;
	y1.r = y2.r, y1.i = y2.i;
/* L119: */
    }
L120:
    i__1 = *m;
    i__2 = *m;
    i__3 = *m - 1;
    q__2.r = a[i__2].r * y[i__3].r - a[i__2].i * y[i__3].i, q__2.i = a[i__2]
	    .r * y[i__3].i + a[i__2].i * y[i__3].r;
    i__4 = *m;
    q__4.r = b[i__4].r - rt, q__4.i = b[i__4].i;
    i__5 = *m;
    q__3.r = q__4.r * y[i__5].r - q__4.i * y[i__5].i, q__3.i = q__4.r * y[
	    i__5].i + q__4.i * y[i__5].r;
    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
    i__1 = *m - 1;
    y[i__1].r = y1.r, y[i__1].i = y1.i;
    iflg = 1;
    goto L102;
L121:
    if (ia <= 0) {
	goto L124;
    } else {
	goto L122;
    }
L122:
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
/* L123: */
    }
L124:
    if (iflg <= 0) {
	goto L125;
    } else {
	goto L102;
    }
L125:
    return 0;
} /* cproc_ */

