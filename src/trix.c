/* trix.f -- translated by f2c (version 12.02.01).
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

/* DECK TRIX */
/* Subroutine */ int trix_(integer *idegbr, integer *idegcr, integer *m, real 
	*a, real *b, real *c__, real *y, real *tcos, real *d__, real *w)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, l;
    static real x, z__;
    static integer kb, kc, ip;
    static real xx;
    static integer mm1, lint;

/* ***BEGIN PROLOGUE  TRIX */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to GENBUN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (TRIX-S, CMPTRX-C) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     Subroutine to solve a system of linear equations where the */
/*     coefficient matrix is a rational function in the matrix given by */
/*     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ). */

/* ***SEE ALSO  GENBUN */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  TRIX */

/* ***FIRST EXECUTABLE STATEMENT  TRIX */
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
    l = (*idegbr + 1) / (*idegcr + 1);
    lint = 1;
    i__1 = *idegbr;
    for (k = 1; k <= i__1; ++k) {
	x = tcos[k];
	if (k != l) {
	    goto L102;
	}
	i__ = *idegbr + lint;
	xx = x - tcos[i__];
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w[i__] = y[i__];
	    y[i__] = xx * y[i__];
/* L101: */
	}
L102:
	z__ = 1.f / (b[1] - x);
	d__[1] = c__[1] * z__;
	y[1] *= z__;
	i__2 = mm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    z__ = 1.f / (b[i__] - x - a[i__] * d__[i__ - 1]);
	    d__[i__] = c__[i__] * z__;
	    y[i__] = (y[i__] - a[i__] * y[i__ - 1]) * z__;
/* L103: */
	}
	z__ = b[*m] - x - a[*m] * d__[mm1];
	if (z__ != 0.f) {
	    goto L104;
	}
	y[*m] = 0.f;
	goto L105;
L104:
	y[*m] = (y[*m] - a[*m] * y[mm1]) / z__;
L105:
	i__2 = mm1;
	for (ip = 1; ip <= i__2; ++ip) {
	    i__ = *m - ip;
	    y[i__] -= d__[i__] * y[i__ + 1];
/* L106: */
	}
	if (k != l) {
	    goto L108;
	}
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] += w[i__];
/* L107: */
	}
	++lint;
	l = lint * kb / kc;
L108:
	;
    }
    return 0;
} /* trix_ */

