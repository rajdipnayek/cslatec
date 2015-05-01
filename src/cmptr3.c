/* cmptr3.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static complex c_b10 = {1.f,0.f};

/* DECK CMPTR3 */
/* Subroutine */ int cmptr3_(integer *m, complex *a, complex *b, complex *c__,
	 integer *k, complex *y1, complex *y2, complex *y3, complex *tcos, 
	complex *d__, complex *w1, complex *w2, complex *w3)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, n;
    static complex x, z__;
    static integer k1, k2, k3, k4, l1, l2, l3, ip;
    static complex xx;
    static integer mm1, k1p1, k2p1, k3p1, k4p1, k2k3k4, kint1, lint1, lint2, 
	    lint3, kint2, kint3;

/* ***BEGIN PROLOGUE  CMPTR3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CMGNBN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (TRI3-S, CMPTR3-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve tridiagonal systems. */

/* ***SEE ALSO  CMGNBN */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890206  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CMPTR3 */

/* ***FIRST EXECUTABLE STATEMENT  CMPTR3 */
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
    l1 = k1p1 / k2p1;
    l2 = k1p1 / k3p1;
    l3 = k1p1 / k4p1;
    lint1 = 1;
    lint2 = 1;
    lint3 = 1;
    kint1 = k1;
    kint2 = kint1 + k2;
    kint3 = kint2 + k3;
L101:
    i__1 = k1;
    for (n = 1; n <= i__1; ++n) {
	i__2 = n;
	x.r = tcos[i__2].r, x.i = tcos[i__2].i;
	if (k2k3k4 == 0) {
	    goto L107;
	}
	if (n != l1) {
	    goto L103;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    w1[i__3].r = y1[i__4].r, w1[i__3].i = y1[i__4].i;
/* L102: */
	}
L103:
	if (n != l2) {
	    goto L105;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    w2[i__3].r = y2[i__4].r, w2[i__3].i = y2[i__4].i;
/* L104: */
	}
L105:
	if (n != l3) {
	    goto L107;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    w3[i__3].r = y3[i__4].r, w3[i__3].i = y3[i__4].i;
/* L106: */
	}
L107:
	q__2.r = b[1].r - x.r, q__2.i = b[1].i - x.i;
	c_div(&q__1, &c_b10, &q__2);
	z__.r = q__1.r, z__.i = q__1.i;
	q__1.r = c__[1].r * z__.r - c__[1].i * z__.i, q__1.i = c__[1].r * 
		z__.i + c__[1].i * z__.r;
	d__[1].r = q__1.r, d__[1].i = q__1.i;
	q__1.r = y1[1].r * z__.r - y1[1].i * z__.i, q__1.i = y1[1].r * z__.i 
		+ y1[1].i * z__.r;
	y1[1].r = q__1.r, y1[1].i = q__1.i;
	q__1.r = y2[1].r * z__.r - y2[1].i * z__.i, q__1.i = y2[1].r * z__.i 
		+ y2[1].i * z__.r;
	y2[1].r = q__1.r, y2[1].i = q__1.i;
	q__1.r = y3[1].r * z__.r - y3[1].i * z__.i, q__1.i = y3[1].r * z__.i 
		+ y3[1].i * z__.r;
	y3[1].r = q__1.r, y3[1].i = q__1.i;
	i__2 = *m;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    q__3.r = b[i__3].r - x.r, q__3.i = b[i__3].i - x.i;
	    i__4 = i__;
	    i__5 = i__ - 1;
	    q__4.r = a[i__4].r * d__[i__5].r - a[i__4].i * d__[i__5].i, 
		    q__4.i = a[i__4].r * d__[i__5].i + a[i__4].i * d__[i__5]
		    .r;
	    q__2.r = q__3.r - q__4.r, q__2.i = q__3.i - q__4.i;
	    c_div(&q__1, &c_b10, &q__2);
	    z__.r = q__1.r, z__.i = q__1.i;
	    i__3 = i__;
	    i__4 = i__;
	    q__1.r = c__[i__4].r * z__.r - c__[i__4].i * z__.i, q__1.i = c__[
		    i__4].r * z__.i + c__[i__4].i * z__.r;
	    d__[i__3].r = q__1.r, d__[i__3].i = q__1.i;
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    i__6 = i__ - 1;
	    q__3.r = a[i__5].r * y1[i__6].r - a[i__5].i * y1[i__6].i, q__3.i =
		     a[i__5].r * y1[i__6].i + a[i__5].i * y1[i__6].r;
	    q__2.r = y1[i__4].r - q__3.r, q__2.i = y1[i__4].i - q__3.i;
	    q__1.r = q__2.r * z__.r - q__2.i * z__.i, q__1.i = q__2.r * z__.i 
		    + q__2.i * z__.r;
	    y1[i__3].r = q__1.r, y1[i__3].i = q__1.i;
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    i__6 = i__ - 1;
	    q__3.r = a[i__5].r * y2[i__6].r - a[i__5].i * y2[i__6].i, q__3.i =
		     a[i__5].r * y2[i__6].i + a[i__5].i * y2[i__6].r;
	    q__2.r = y2[i__4].r - q__3.r, q__2.i = y2[i__4].i - q__3.i;
	    q__1.r = q__2.r * z__.r - q__2.i * z__.i, q__1.i = q__2.r * z__.i 
		    + q__2.i * z__.r;
	    y2[i__3].r = q__1.r, y2[i__3].i = q__1.i;
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    i__6 = i__ - 1;
	    q__3.r = a[i__5].r * y3[i__6].r - a[i__5].i * y3[i__6].i, q__3.i =
		     a[i__5].r * y3[i__6].i + a[i__5].i * y3[i__6].r;
	    q__2.r = y3[i__4].r - q__3.r, q__2.i = y3[i__4].i - q__3.i;
	    q__1.r = q__2.r * z__.r - q__2.i * z__.i, q__1.i = q__2.r * z__.i 
		    + q__2.i * z__.r;
	    y3[i__3].r = q__1.r, y3[i__3].i = q__1.i;
/* L108: */
	}
	i__2 = mm1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__ = *m - ip;
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    i__6 = i__ + 1;
	    q__2.r = d__[i__5].r * y1[i__6].r - d__[i__5].i * y1[i__6].i, 
		    q__2.i = d__[i__5].r * y1[i__6].i + d__[i__5].i * y1[i__6]
		    .r;
	    q__1.r = y1[i__4].r - q__2.r, q__1.i = y1[i__4].i - q__2.i;
	    y1[i__3].r = q__1.r, y1[i__3].i = q__1.i;
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    i__6 = i__ + 1;
	    q__2.r = d__[i__5].r * y2[i__6].r - d__[i__5].i * y2[i__6].i, 
		    q__2.i = d__[i__5].r * y2[i__6].i + d__[i__5].i * y2[i__6]
		    .r;
	    q__1.r = y2[i__4].r - q__2.r, q__1.i = y2[i__4].i - q__2.i;
	    y2[i__3].r = q__1.r, y2[i__3].i = q__1.i;
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    i__6 = i__ + 1;
	    q__2.r = d__[i__5].r * y3[i__6].r - d__[i__5].i * y3[i__6].i, 
		    q__2.i = d__[i__5].r * y3[i__6].i + d__[i__5].i * y3[i__6]
		    .r;
	    q__1.r = y3[i__4].r - q__2.r, q__1.i = y3[i__4].i - q__2.i;
	    y3[i__3].r = q__1.r, y3[i__3].i = q__1.i;
/* L109: */
	}
	if (k2k3k4 == 0) {
	    goto L115;
	}
	if (n != l1) {
	    goto L111;
	}
	i__ = lint1 + kint1;
	i__2 = i__;
	q__1.r = x.r - tcos[i__2].r, q__1.i = x.i - tcos[i__2].i;
	xx.r = q__1.r, xx.i = q__1.i;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    q__2.r = xx.r * y1[i__4].r - xx.i * y1[i__4].i, q__2.i = xx.r * 
		    y1[i__4].i + xx.i * y1[i__4].r;
	    i__5 = i__;
	    q__1.r = q__2.r + w1[i__5].r, q__1.i = q__2.i + w1[i__5].i;
	    y1[i__3].r = q__1.r, y1[i__3].i = q__1.i;
/* L110: */
	}
	++lint1;
	l1 = lint1 * k1p1 / k2p1;
L111:
	if (n != l2) {
	    goto L113;
	}
	i__ = lint2 + kint2;
	i__2 = i__;
	q__1.r = x.r - tcos[i__2].r, q__1.i = x.i - tcos[i__2].i;
	xx.r = q__1.r, xx.i = q__1.i;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    q__2.r = xx.r * y2[i__4].r - xx.i * y2[i__4].i, q__2.i = xx.r * 
		    y2[i__4].i + xx.i * y2[i__4].r;
	    i__5 = i__;
	    q__1.r = q__2.r + w2[i__5].r, q__1.i = q__2.i + w2[i__5].i;
	    y2[i__3].r = q__1.r, y2[i__3].i = q__1.i;
/* L112: */
	}
	++lint2;
	l2 = lint2 * k1p1 / k3p1;
L113:
	if (n != l3) {
	    goto L115;
	}
	i__ = lint3 + kint3;
	i__2 = i__;
	q__1.r = x.r - tcos[i__2].r, q__1.i = x.i - tcos[i__2].i;
	xx.r = q__1.r, xx.i = q__1.i;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    q__2.r = xx.r * y3[i__4].r - xx.i * y3[i__4].i, q__2.i = xx.r * 
		    y3[i__4].i + xx.i * y3[i__4].r;
	    i__5 = i__;
	    q__1.r = q__2.r + w3[i__5].r, q__1.i = q__2.i + w3[i__5].i;
	    y3[i__3].r = q__1.r, y3[i__3].i = q__1.i;
/* L114: */
	}
	++lint3;
	l3 = lint3 * k1p1 / k4p1;
L115:
	;
    }
    return 0;
} /* cmptr3_ */

