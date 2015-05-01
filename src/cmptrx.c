/* cmptrx.f -- translated by f2c (version 12.02.01).
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

static complex c_b5 = {1.f,0.f};

/* DECK CMPTRX */
/* Subroutine */ int cmptrx_(integer *idegbr, integer *idegcr, integer *m, 
	complex *a, complex *b, complex *c__, complex *y, complex *tcos, 
	complex *d__, complex *w)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, k, l;
    static complex x, z__;
    static integer kb, kc, ip;
    static complex xx;
    static integer mm1, lint;

/* ***BEGIN PROLOGUE  CMPTRX */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CMGNBN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      COMPLEX (TRIX-S, CMPTRX-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve a system of linear equations where the */
/*     coefficient matrix is a rational function in the matrix given by */
/*     tridiagonal  ( . . . , A(I), B(I), C(I), . . . ). */

/* ***SEE ALSO  CMGNBN */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CMPTRX */

/* ***FIRST EXECUTABLE STATEMENT  CMPTRX */
    /* Parameter adjustments */
    --w;
    --d__;
    --tcos;
    --y;
    --c__;
    --b;
    --a;

    /* Function Body */
    mm1 = *m - 1;
    kb = *idegbr + 1;
    kc = *idegcr + 1;
    l = kb / kc;
    lint = 1;
    i__1 = *idegbr;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	x.r = tcos[i__2].r, x.i = tcos[i__2].i;
	if (k != l) {
	    goto L102;
	}
	i__ = *idegbr + lint;
	i__2 = i__;
	q__1.r = x.r - tcos[i__2].r, q__1.i = x.i - tcos[i__2].i;
	xx.r = q__1.r, xx.i = q__1.i;
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    w[i__3].r = y[i__4].r, w[i__3].i = y[i__4].i;
	    i__3 = i__;
	    i__4 = i__;
	    q__1.r = xx.r * y[i__4].r - xx.i * y[i__4].i, q__1.i = xx.r * y[
		    i__4].i + xx.i * y[i__4].r;
	    y[i__3].r = q__1.r, y[i__3].i = q__1.i;
/* L101: */
	}
L102:
	q__2.r = b[1].r - x.r, q__2.i = b[1].i - x.i;
	c_div(&q__1, &c_b5, &q__2);
	z__.r = q__1.r, z__.i = q__1.i;
	q__1.r = c__[1].r * z__.r - c__[1].i * z__.i, q__1.i = c__[1].r * 
		z__.i + c__[1].i * z__.r;
	d__[1].r = q__1.r, d__[1].i = q__1.i;
	q__1.r = y[1].r * z__.r - y[1].i * z__.i, q__1.i = y[1].r * z__.i + y[
		1].i * z__.r;
	y[1].r = q__1.r, y[1].i = q__1.i;
	i__2 = mm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    q__3.r = b[i__3].r - x.r, q__3.i = b[i__3].i - x.i;
	    i__4 = i__;
	    i__5 = i__ - 1;
	    q__4.r = a[i__4].r * d__[i__5].r - a[i__4].i * d__[i__5].i, 
		    q__4.i = a[i__4].r * d__[i__5].i + a[i__4].i * d__[i__5]
		    .r;
	    q__2.r = q__3.r - q__4.r, q__2.i = q__3.i - q__4.i;
	    c_div(&q__1, &c_b5, &q__2);
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
	    q__3.r = a[i__5].r * y[i__6].r - a[i__5].i * y[i__6].i, q__3.i = 
		    a[i__5].r * y[i__6].i + a[i__5].i * y[i__6].r;
	    q__2.r = y[i__4].r - q__3.r, q__2.i = y[i__4].i - q__3.i;
	    q__1.r = q__2.r * z__.r - q__2.i * z__.i, q__1.i = q__2.r * z__.i 
		    + q__2.i * z__.r;
	    y[i__3].r = q__1.r, y[i__3].i = q__1.i;
/* L103: */
	}
	i__2 = *m;
	q__2.r = b[i__2].r - x.r, q__2.i = b[i__2].i - x.i;
	i__3 = *m;
	i__4 = mm1;
	q__3.r = a[i__3].r * d__[i__4].r - a[i__3].i * d__[i__4].i, q__3.i = 
		a[i__3].r * d__[i__4].i + a[i__3].i * d__[i__4].r;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	z__.r = q__1.r, z__.i = q__1.i;
	if (c_abs(&z__) != 0.f) {
	    goto L104;
	}
	i__2 = *m;
	y[i__2].r = 0.f, y[i__2].i = 0.f;
	goto L105;
L104:
	i__2 = *m;
	i__3 = *m;
	i__4 = *m;
	i__5 = mm1;
	q__3.r = a[i__4].r * y[i__5].r - a[i__4].i * y[i__5].i, q__3.i = a[
		i__4].r * y[i__5].i + a[i__4].i * y[i__5].r;
	q__2.r = y[i__3].r - q__3.r, q__2.i = y[i__3].i - q__3.i;
	c_div(&q__1, &q__2, &z__);
	y[i__2].r = q__1.r, y[i__2].i = q__1.i;
L105:
	i__2 = mm1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__ = *m - ip;
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    i__6 = i__ + 1;
	    q__2.r = d__[i__5].r * y[i__6].r - d__[i__5].i * y[i__6].i, 
		    q__2.i = d__[i__5].r * y[i__6].i + d__[i__5].i * y[i__6]
		    .r;
	    q__1.r = y[i__4].r - q__2.r, q__1.i = y[i__4].i - q__2.i;
	    y[i__3].r = q__1.r, y[i__3].i = q__1.i;
/* L106: */
	}
	if (k != l) {
	    goto L108;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__;
	    q__1.r = y[i__4].r + w[i__5].r, q__1.i = y[i__4].i + w[i__5].i;
	    y[i__3].r = q__1.r, y[i__3].i = q__1.i;
/* L107: */
	}
	++lint;
	l = lint * kb / kc;
L108:
	;
    }
    return 0;
} /* cmptrx_ */

