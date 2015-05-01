/* cpevlr.f -- translated by f2c (version 12.02.01).
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

/* DECK CPEVLR */
/* Subroutine */ int cpevlr_(integer *n, integer *m, real *a, real *x, real *
	c__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j;
    static real ci;
    static integer np1;
    static real cim1;
    static integer mini;

/* ***BEGIN PROLOGUE  CPEVLR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to CPZERO */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (CPEVLR-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***SEE ALSO  CPZERO */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   810223  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  CPEVLR */
/* ***FIRST EXECUTABLE STATEMENT  CPEVLR */
    /* Parameter adjustments */
    --c__;
    --a;

    /* Function Body */
    np1 = *n + 1;
    i__1 = np1;
    for (j = 1; j <= i__1; ++j) {
	ci = 0.f;
	cim1 = a[j];
/* Computing MIN */
	i__2 = *m + 1, i__3 = *n + 2 - j;
	mini = min(i__2,i__3);
	i__2 = mini;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (j != 1) {
		ci = c__[i__];
	    }
	    if (i__ != 1) {
		cim1 = c__[i__ - 1];
	    }
	    c__[i__] = cim1 + *x * ci;
/* L1: */
	}
    }
    return 0;
} /* cpevlr_ */

