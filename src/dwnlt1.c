/* dwnlt1.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;

/* DECK DWNLT1 */
/* Subroutine */ int dwnlt1_(integer *i__, integer *lend, integer *mend, 
	integer *ir, integer *mdw, logical *recalc, integer *imax, doublereal 
	*hbar, doublereal *h__, doublereal *scale, doublereal *w)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer j, k;
    extern integer idamax_(integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DWNLT1 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to WNLIT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     To update the column Sum Of Squares and find the pivot column. */
/*     The column Sum of Squares Vector will be updated at each step. */
/*     When numerically necessary, these values will be recomputed. */

/* ***SEE ALSO  DWNLIT */
/* ***ROUTINES CALLED  IDAMAX */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890620  Code extracted from WNLIT and made a subroutine.  (RWC)) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DWNLT1 */



/* ***FIRST EXECUTABLE STATEMENT  DWNLT1 */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --h__;
    --scale;

    /* Function Body */
    if (*ir != 1 && ! (*recalc)) {

/*        Update column SS=sum of squares. */

	i__1 = *lend;
	for (j = *i__; j <= i__1; ++j) {
/* Computing 2nd power */
	    d__1 = w[*ir - 1 + j * w_dim1];
	    h__[j] -= scale[*ir - 1] * (d__1 * d__1);
/* L10: */
	}

/*        Test for numerical accuracy. */

	i__1 = *lend - *i__ + 1;
	*imax = idamax_(&i__1, &h__[*i__], &c__1) + *i__ - 1;
	*recalc = *hbar + h__[*imax] * .001f == *hbar;
    }

/*     If required, recalculate column SS, using rows IR through MEND. */

    if (*recalc) {
	i__1 = *lend;
	for (j = *i__; j <= i__1; ++j) {
	    h__[j] = 0.;
	    i__2 = *mend;
	    for (k = *ir; k <= i__2; ++k) {
/* Computing 2nd power */
		d__1 = w[k + j * w_dim1];
		h__[j] += scale[k] * (d__1 * d__1);
/* L20: */
	    }
/* L30: */
	}

/*        Find column with largest SS. */

	i__1 = *lend - *i__ + 1;
	*imax = idamax_(&i__1, &h__[*i__], &c__1) + *i__ - 1;
	*hbar = h__[*imax];
    }
    return 0;
} /* dwnlt1_ */

