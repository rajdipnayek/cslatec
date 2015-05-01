/* pgsf.f -- translated by f2c (version 12.02.01).
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

/* DECK PGSF */
doublereal pgsf_(real *x, integer *iz, real *c__, real *a, real *bh)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer j;
    static real dd, fsg, hsg;

/* ***BEGIN PROLOGUE  PGSF */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PGSF-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PGSF */
/* ***FIRST EXECUTABLE STATEMENT  PGSF */
    /* Parameter adjustments */
    --bh;
    --a;
    --c__;

    /* Function Body */
    fsg = 1.f;
    hsg = 1.f;
    i__1 = *iz;
    for (j = 1; j <= i__1; ++j) {
	dd = 1.f / (*x - bh[j]);
	fsg = fsg * a[j] * dd;
	hsg = hsg * c__[j] * dd;
/* L101: */
    }
    if (*iz % 2 != 0) {
	goto L103;
    } else {
	goto L102;
    }
L102:
    ret_val = 1.f - fsg - hsg;
    return ret_val;
L103:
    ret_val = fsg + 1.f + hsg;
    return ret_val;
} /* pgsf_ */

