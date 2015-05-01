/* r9upak.f -- translated by f2c (version 12.02.01).
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

/* DECK R9UPAK */
/* Subroutine */ int r9upak_(real *x, real *y, integer *n)
{
    /* Local variables */
    static real absx;

/* ***BEGIN PROLOGUE  R9UPAK */
/* ***PURPOSE  Unpack a floating point number X so that X = Y*2**N. */
/* ***LIBRARY   SLATEC (FNLIB) */
/* ***CATEGORY  A6B */
/* ***TYPE      SINGLE PRECISION (R9UPAK-S, D9UPAK-D) */
/* ***KEYWORDS  FNLIB, UNPACK */
/* ***AUTHOR  Fullerton, W., (LANL) */
/* ***DESCRIPTION */

/*   Unpack a floating point number X so that X = Y*2.0**N, where */
/*   0.5 .LE. ABS(Y) .LT. 1.0. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780701  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/* ***END PROLOGUE  R9UPAK */
/* ***FIRST EXECUTABLE STATEMENT  R9UPAK */
    absx = dabs(*x);
    *n = 0;
    if (*x == 0.f) {
	goto L30;
    }

L10:
    if (absx >= .5f) {
	goto L20;
    }
    --(*n);
    absx *= 2.f;
    goto L10;

L20:
    if (absx < 1.f) {
	goto L30;
    }
    ++(*n);
    absx *= .5f;
    goto L20;

L30:
    *y = r_sign(&absx, x);
    return 0;

} /* r9upak_ */

