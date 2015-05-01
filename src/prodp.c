/* prodp.f -- translated by f2c (version 12.02.01).
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

/* DECK PRODP */
/* Subroutine */ int prodp_(integer *nd, real *bd, integer *nm1, real *bm1, 
	integer *nm2, real *bm2, integer *na, real *aa, real *x, real *y, 
	integer *m, real *a, real *b, real *c__, real *d__, real *u, real *w)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static integer j, k;
    static real v;
    static integer m1, m2, ia;
    static real bh;
    static integer id;
    static real am;
    static integer mm;
    static real rt, ym;
    static integer mm2;
    static real den;
    static integer ibr;

/* ***BEGIN PROLOGUE  PRODP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BLKTRI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PRODP-S, PROCP-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* PRODP applies a sequence of matrix operations to the vector X and */
/* stores the result in Y (periodic boundary conditions). */

/* BD,BM1,BM2 are arrays containing roots of certain B polynomials. */
/* ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively. */
/* AA         Array containing scalar multipliers of the vector X. */
/* NA         is the length of the array AA. */
/* X,Y        The matrix operations are applied to X and the result is Y. */
/* A,B,C      are arrays which contain the tridiagonal matrix. */
/* M          is the order of the matrix. */
/* D,W,U      are working arrays. */
/* IS         determines whether or not a change in sign is made. */

/* ***SEE ALSO  BLKTRI */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PRODP */

/* ***FIRST EXECUTABLE STATEMENT  PRODP */
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
	y[j] = x[j];
	w[j] = y[j];
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
	y[j] = rt * w[j];
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

    bh = b[*m] - rt;
    ym = y[*m];
    den = b[1] - rt;
    d__[1] = c__[1] / den;
    u[1] = a[1] / den;
    w[1] = y[1] / den;
    v = c__[*m];
    if (mm2 - 2 >= 0) {
	goto L107;
    } else {
	goto L109;
    }
L107:
    i__1 = mm2;
    for (j = 2; j <= i__1; ++j) {
	den = b[j] - rt - a[j] * d__[j - 1];
	d__[j] = c__[j] / den;
	u[j] = -a[j] * u[j - 1] / den;
	w[j] = (y[j] - a[j] * w[j - 1]) / den;
	bh -= v * u[j - 1];
	ym -= v * w[j - 1];
	v = -v * d__[j - 1];
/* L108: */
    }
L109:
    den = b[*m - 1] - rt - a[*m - 1] * d__[*m - 2];
    d__[*m - 1] = (c__[*m - 1] - a[*m - 1] * u[*m - 2]) / den;
    w[*m - 1] = (y[*m - 1] - a[*m - 1] * w[*m - 2]) / den;
    am = a[*m] - v * d__[*m - 2];
    bh -= v * u[*m - 2];
    ym -= v * w[*m - 2];
    den = bh - am * d__[*m - 1];
    if (den != 0.f) {
	goto L110;
    } else {
	goto L111;
    }
L110:
    w[*m] = (ym - am * w[*m - 1]) / den;
    goto L112;
L111:
    w[*m] = 1.f;
L112:
    w[*m - 1] -= d__[*m - 1] * w[*m];
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	k = *m - j;
	w[k] = w[k] - d__[k] * w[k + 1] - u[k] * w[*m];
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
	y[j] = w[j];
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
	y[j] += rt * w[j];
/* L127: */
    }
    goto L102;
L128:
    return 0;
} /* prodp_ */

