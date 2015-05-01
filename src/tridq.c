/* tridq.f -- translated by f2c (version 12.02.01).
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

/* DECK TRIDQ */
/* Subroutine */ int tridq_(integer *mr, real *a, real *b, real *c__, real *y,
	 real *d__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m;
    static real z__;
    static integer ip, mm1;

/* ***BEGIN PROLOGUE  TRIDQ */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to POIS3D */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (TRIDQ-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  POIS3D */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900308  Renamed routine from TRID to TRIDQ.  (WRB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  TRIDQ */
/* ***FIRST EXECUTABLE STATEMENT  TRIDQ */
    /* Parameter adjustments */
    --d__;
    --y;
    --c__;
    --b;
    --a;

    /* Function Body */
    m = *mr;
    mm1 = m - 1;
    z__ = 1.f / b[1];
    d__[1] = c__[1] * z__;
    y[1] *= z__;
    i__1 = mm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	z__ = 1.f / (b[i__] - a[i__] * d__[i__ - 1]);
	d__[i__] = c__[i__] * z__;
	y[i__] = (y[i__] - a[i__] * y[i__ - 1]) * z__;
/* L101: */
    }
    z__ = b[m] - a[m] * d__[mm1];
    if (z__ != 0.f) {
	goto L102;
    }
    y[m] = 0.f;
    goto L103;
L102:
    y[m] = (y[m] - a[m] * y[mm1]) / z__;
L103:
    i__1 = mm1;
    for (ip = 1; ip <= i__1; ++ip) {
	i__ = m - ip;
	y[i__] -= d__[i__] * y[i__ + 1];
/* L104: */
    }
    return 0;
} /* tridq_ */

