/* prod.f -- translated by f2c (version 12.02.01).
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

/* DECK PROD */
/* Subroutine */ int prod_(integer *nd, real *bd, integer *nm1, real *bm1, 
	integer *nm2, real *bm2, integer *na, real *aa, real *x, real *y, 
	integer *m, real *a, real *b, real *c__, real *d__, real *w, real *u)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static integer j, k, m1, m2, ia, id, mm;
    static real rt, den;
    static integer ibr;

/* ***BEGIN PROLOGUE  PROD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to BLKTRI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PROD-S, PROC-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/* PROD applies a sequence of matrix operations to the vector X and */
/* stores the result in Y. */

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
/* ***END PROLOGUE  PROD */

/* ***FIRST EXECUTABLE STATEMENT  PROD */
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
	w[j] = x[j];
	y[j] = w[j];
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
	y[j] = rt * w[j];
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

    d__[*m] = a[*m] / (b[*m] - rt);
    w[*m] = y[*m] / (b[*m] - rt);
    i__1 = mm;
    for (j = 2; j <= i__1; ++j) {
	k = *m - j;
	den = b[k + 1] - rt - c__[k + 1] * d__[k + 2];
	d__[k + 1] = a[k + 1] / den;
	w[k + 1] = (y[k + 1] - c__[k + 1] * w[k + 2]) / den;
/* L107: */
    }
    den = b[1] - rt - c__[1] * d__[2];
    w[1] = 1.f;
    if (den != 0.f) {
	goto L108;
    } else {
	goto L109;
    }
L108:
    w[1] = (y[1] - c__[1] * w[2]) / den;
L109:
    i__1 = *m;
    for (j = 2; j <= i__1; ++j) {
	w[j] -= d__[j] * w[j - 1];
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
	y[j] = w[j];
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
	y[j] += rt * w[j];
/* L124: */
    }
    goto L102;
L125:
    return 0;
} /* prod_ */

