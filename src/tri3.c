/* tri3.f -- translated by f2c (version 12.02.01).
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

/* DECK TRI3 */
/* Subroutine */ int tri3_(integer *m, real *a, real *b, real *c__, integer *
	k, real *y1, real *y2, real *y3, real *tcos, real *d__, real *w1, 
	real *w2, real *w3)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, n;
    static real x, z__;
    static integer k1, k2, k3, k4, l1, l2, l3, ip;
    static real xx;
    static integer mm1, k1p1, k2p1, k3p1, k4p1, k2k3k4, kint1, lint1, lint2, 
	    lint3, kint2, kint3;

/* ***BEGIN PROLOGUE  TRI3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (TRI3-S, CMPTR3-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve three linear systems whose common coefficient */
/*     matrix is a rational function in the matrix given by */

/*                  TRIDIAGONAL (...,A(I),B(I),C(I),...) */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  TRI3 */

/* ***FIRST EXECUTABLE STATEMENT  TRI3 */
    /* Parameter adjustments */
    --w3;
    --w2;
    --w1;
    --d__;
    --tcos;
    --y3;
    --y2;
    --y1;
    --k;
    --c__;
    --b;
    --a;

    /* Function Body */
    mm1 = *m - 1;
    k1 = k[1];
    k2 = k[2];
    k3 = k[3];
    k4 = k[4];
    k1p1 = k1 + 1;
    k2p1 = k2 + 1;
    k3p1 = k3 + 1;
    k4p1 = k4 + 1;
    k2k3k4 = k2 + k3 + k4;
    if (k2k3k4 == 0) {
	goto L101;
    }
    l1 = (k1 + 1) / (k2 + 1);
    l2 = (k1 + 1) / (k3 + 1);
    l3 = (k1 + 1) / (k4 + 1);
    lint1 = 1;
    lint2 = 1;
    lint3 = 1;
    kint1 = k1;
    kint2 = kint1 + k2;
    kint3 = kint2 + k3;
L101:
    i__1 = k1;
    for (n = 1; n <= i__1; ++n) {
	x = tcos[n];
	if (k2k3k4 == 0) {
	    goto L107;
	}
	if (n != l1) {
	    goto L103;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w1[i__] = y1[i__];
/* L102: */
	}
L103:
	if (n != l2) {
	    goto L105;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w2[i__] = y2[i__];
/* L104: */
	}
L105:
	if (n != l3) {
	    goto L107;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w3[i__] = y3[i__];
/* L106: */
	}
L107:
	z__ = 1.f / (b[1] - x);
	d__[1] = c__[1] * z__;
	y1[1] *= z__;
	y2[1] *= z__;
	y3[1] *= z__;
	i__2 = *m;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    z__ = 1.f / (b[i__] - x - a[i__] * d__[i__ - 1]);
	    d__[i__] = c__[i__] * z__;
	    y1[i__] = (y1[i__] - a[i__] * y1[i__ - 1]) * z__;
	    y2[i__] = (y2[i__] - a[i__] * y2[i__ - 1]) * z__;
	    y3[i__] = (y3[i__] - a[i__] * y3[i__ - 1]) * z__;
/* L108: */
	}
	i__2 = mm1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__ = *m - ip;
	    y1[i__] -= d__[i__] * y1[i__ + 1];
	    y2[i__] -= d__[i__] * y2[i__ + 1];
	    y3[i__] -= d__[i__] * y3[i__ + 1];
/* L109: */
	}
	if (k2k3k4 == 0) {
	    goto L115;
	}
	if (n != l1) {
	    goto L111;
	}
	i__ = lint1 + kint1;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y1[i__] = xx * y1[i__] + w1[i__];
/* L110: */
	}
	++lint1;
	l1 = lint1 * k1p1 / k2p1;
L111:
	if (n != l2) {
	    goto L113;
	}
	i__ = lint2 + kint2;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y2[i__] = xx * y2[i__] + w2[i__];
/* L112: */
	}
	++lint2;
	l2 = lint2 * k1p1 / k3p1;
L113:
	if (n != l3) {
	    goto L115;
	}
	i__ = lint3 + kint3;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y3[i__] = xx * y3[i__] + w3[i__];
/* L114: */
	}
	++lint3;
	l3 = lint3 * k1p1 / k4p1;
L115:
	;
    }
    return 0;
} /* tri3_ */

