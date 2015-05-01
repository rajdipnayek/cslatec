/* dbdiff.f -- translated by f2c (version 12.02.01).
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

/* DECK DBDIFF */
/* Subroutine */ int dbdiff_(integer *l, doublereal *v)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k;

/* ***BEGIN PROLOGUE  DBDIFF */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBSKIN */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (BDIFF-S, DBDIFF-D) */
/* ***AUTHOR  Amos, D. E., (SNLA) */
/* ***DESCRIPTION */

/*     DBDIFF computes the sum of B(L,K)*V(K)*(-1)**K where B(L,K) */
/*     are the binomial coefficients.  Truncated sums are computed by */
/*     setting last part of the V vector to zero. On return, the binomial */
/*     sum is in V(L). */

/* ***SEE ALSO  DBSKIN */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   820601  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  DBDIFF */

/* ***FIRST EXECUTABLE STATEMENT  DBDIFF */
    /* Parameter adjustments */
    --v;

    /* Function Body */
    if (*l == 1) {
	return 0;
    }
    i__1 = *l;
    for (j = 2; j <= i__1; ++j) {
	k = *l;
	i__2 = *l;
	for (i__ = j; i__ <= i__2; ++i__) {
	    v[k] = v[k - 1] - v[k];
	    --k;
/* L10: */
	}
/* L20: */
    }
    return 0;
} /* dbdiff_ */

