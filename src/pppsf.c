/* pppsf.f -- translated by f2c (version 12.02.01).
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

/* DECK PPPSF */
doublereal pppsf_(real *x, integer *iz, real *c__, real *a, real *bh)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer j;
    static real sum;

/* ***BEGIN PROLOGUE  PPPSF */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CBLKTR */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (PPPSF-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  CBLKTR */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  PPPSF */
/* ***FIRST EXECUTABLE STATEMENT  PPPSF */
    /* Parameter adjustments */
    --bh;
    --a;
    --c__;

    /* Function Body */
    sum = 0.f;
    i__1 = *iz;
    for (j = 1; j <= i__1; ++j) {
	sum += 1.f / (*x - bh[j]);
/* L101: */
    }
    ret_val = sum;
    return ret_val;
} /* pppsf_ */

