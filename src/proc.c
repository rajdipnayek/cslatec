/* proc.f -- translated by f2c (version 12.02.01).
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

/* DECK PROC */
/* Subroutine */ int proc_(integer *nd, real *bd, integer *nm1, real *bm1, 
	integer *nm2, real *bm2, integer *na, real *aa, complex *x, complex *
	y, integer *m, complex *a, complex *b, complex *c__, complex *d__, 
	complex *w, complex *u)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1, q__2, q__3;

    /* Local variables */
    static integer j, k, m1, m2, ia, id, mm;
    static real rt;
    static complex den;
    static integer ibr;

/* ***BEGIN PROLOGUE  PROC */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (PROD-S, PROC-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* PROC applies a sequence of matrix operations to the vector X and */
/*  stores the result in Y. */
/* BD,BM1,BM2 are arrays containing roots of certain B polynomials. */
/* ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively. */
/* AA         Array containing scalar multipliers of the vector X. */
/* NA         is the length of the array AA. */
/* X,Y        The matrix operations are applied to X and the result is Y. */
/* A,B,C      are arrays which contain the tridiagonal matrix. */
/* M          is the order of the matrix. */
/* D,W,U      are working arrays. */
/* IS         determines whether or not a change in sign is made. */

/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PROC */

/* ***FIRST EXECUTABLE STATEMENT  PROC */
    /* Parameter adjustments */
    --u;
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
	w[i__2].r = x[i__3].r, w[i__2].i = x[i__3].i;
	i__2 = j;
	i__3 = j;
	y[i__2].r = w[i__3].r, y[i__2].i = w[i__3].i;
/* L101: */
    }
    mm = *m - 1;
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

/* SCALAR MULTIPLICATION */

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
	goto L125;
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
    i__2 = *m;
    q__2.r = b[i__2].r - rt, q__2.i = b[i__2].i;
    c_div(&q__1, &a[*m], &q__2);
    d__[i__1].r = q__1.r, d__[i__1].i = q__1.i;
    i__1 = *m;
    i__2 = *m;
    q__2.r = b[i__2].r - rt, q__2.i = b[i__2].i;
    c_div(&q__1, &y[*m], &q__2);
    w[i__1].r = q__1.r, w[i__1].i = q__1.i;
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	k = *m - j;
	i__2 = k + 1;
	q__2.r = b[i__2].r - rt, q__2.i = b[i__2].i;
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
/* L107: */
    }
    q__2.r = b[1].r - rt, q__2.i = b[1].i;
    q__3.r = c__[1].r * d__[2].r - c__[1].i * d__[2].i, q__3.i = c__[1].r * 
	    d__[2].i + c__[1].i * d__[2].r;
    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
    den.r = q__1.r, den.i = q__1.i;
    w[1].r = 1.f, w[1].i = 0.f;
    if (c_abs(&den) != 0.f) {
	goto L108;
    } else {
	goto L109;
    }
L108:
    q__3.r = c__[1].r * w[2].r - c__[1].i * w[2].i, q__3.i = c__[1].r * w[2]
	    .i + c__[1].i * w[2].r;
    q__2.r = y[1].r - q__3.r, q__2.i = y[1].i - q__3.i;
    c_div(&q__1, &q__2, &den);
    w[1].r = q__1.r, w[1].i = q__1.i;
L109:
    i__1 = *m;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	i__4 = j;
	i__5 = j - 1;
	q__2.r = d__[i__4].r * w[i__5].r - d__[i__4].i * w[i__5].i, q__2.i = 
		d__[i__4].r * w[i__5].i + d__[i__4].i * w[i__5].r;
	q__1.r = w[i__3].r - q__2.r, q__1.i = w[i__3].i - q__2.i;
	w[i__2].r = q__1.r, w[i__2].i = q__1.i;
/* L110: */
    }
    if (*na <= 0) {
	goto L113;
    } else {
	goto L102;
    }
L111:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	y[i__2].r = w[i__3].r, y[i__2].i = w[i__3].i;
/* L112: */
    }
    ibr = 1;
    goto L102;
L113:
    if (m1 <= 0) {
	goto L114;
    } else {
	goto L115;
    }
L114:
    if (m2 <= 0) {
	goto L111;
    } else {
	goto L120;
    }
L115:
    if (m2 <= 0) {
	goto L117;
    } else {
	goto L116;
    }
L116:
    if ((r__1 = bm1[m1], dabs(r__1)) - (r__2 = bm2[m2], dabs(r__2)) <= 0.f) {
	goto L120;
    } else {
	goto L117;
    }
L117:
    if (ibr <= 0) {
	goto L118;
    } else {
	goto L119;
    }
L118:
    if ((r__1 = bm1[m1] - bd[id], dabs(r__1)) - (r__2 = bm1[m1] - rt, dabs(
	    r__2)) >= 0.f) {
	goto L119;
    } else {
	goto L111;
    }
L119:
    rt -= bm1[m1];
    --m1;
    goto L123;
L120:
    if (ibr <= 0) {
	goto L121;
    } else {
	goto L122;
    }
L121:
    if ((r__1 = bm2[m2] - bd[id], dabs(r__1)) - (r__2 = bm2[m2] - rt, dabs(
	    r__2)) >= 0.f) {
	goto L122;
    } else {
	goto L111;
    }
L122:
    rt -= bm2[m2];
    --m2;
L123:
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	i__4 = j;
	q__2.r = rt * w[i__4].r, q__2.i = rt * w[i__4].i;
	q__1.r = y[i__3].r + q__2.r, q__1.i = y[i__3].i + q__2.i;
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
/* L124: */
    }
    goto L102;
L125:
    return 0;
} /* proc_ */

