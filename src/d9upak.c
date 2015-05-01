/* d9upak.f -- translated by f2c (version 12.02.01).
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

/* DECK D9UPAK */
/* Subroutine */ int d9upak_(doublereal *x, doublereal *y, integer *n)
{
    /* Local variables */
    static doublereal absx;

/* ***BEGIN PROLOGUE  D9UPAK */
/* ***PURPOSE  Unpack a floating point number X so that X = Y*2**N. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  A6B */
/* ***TYPE      DOUBLE PRECISION (R9UPAK-S, D9UPAK-D) */
/* ***KEYWORDS  FNLIB, UNPACK */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*   Unpack a floating point number X so that X = Y*2.0**N, where */
/*   0.5 .LE. ABS(Y) .LT. 1.0. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780701  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900820  Corrected code to find Y between 0.5 and 1.0 rather than */
/*           between 0.05 and 1.0.  (WRB) */
/* ***END PROLOGUE  D9UPAK */
/* ***FIRST EXECUTABLE STATEMENT  D9UPAK */
    absx = abs(*x);
    *n = 0;
    if (*x == 0.) {
	goto L30;
    }

L10:
    if (absx >= .5) {
	goto L20;
    }
    --(*n);
    absx *= 2.;
    goto L10;

L20:
    if (absx < 1.) {
	goto L30;
    }
    ++(*n);
    absx *= .5;
    goto L20;

L30:
    *y = d_sign(&absx, x);
    return 0;

} /* d9upak_ */

