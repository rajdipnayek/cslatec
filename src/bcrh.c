/* bcrh.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    integer npp, k;
    real eps, cnv;
    integer nm, ncmplx, ik;
} ccblk_;

#define ccblk_1 ccblk_

/* DECK BCRH */
doublereal bcrh_(real *xll, real *xrr, integer *iz, real *c__, real *a, real *
	bh, E_fp f, real *sgn)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Local variables */
    static real x, dx, xl, xr;

/* ***BEGIN PROLOGUE  BCRH */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (BCRH-S, BSRH-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  (NONE) */
/* ***COMMON BLOCKS    CCBLK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  BCRH */
/* ***FIRST EXECUTABLE STATEMENT  BCRH */
    /* Parameter adjustments */
    --bh;
    --a;
    --c__;

    /* Function Body */
    xl = *xll;
    xr = *xrr;
    dx = (r__1 = xr - xl, dabs(r__1)) * .5f;
L101:
    x = (xl + xr) * .5f;
    if ((r__1 = *sgn * (*f)(&x, iz, &c__[1], &a[1], &bh[1])) < 0.f) {
	goto L103;
    } else if (r__1 == 0) {
	goto L105;
    } else {
	goto L102;
    }
L102:
    xr = x;
    goto L104;
L103:
    xl = x;
L104:
    dx *= .5f;
    if (dx - ccblk_1.cnv <= 0.f) {
	goto L105;
    } else {
	goto L101;
    }
L105:
    ret_val = (xl + xr) * .5f;
    return ret_val;
} /* bcrh_ */

