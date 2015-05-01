/* dwnlt3.f -- translated by f2c (version 12.02.01).
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

/* DECK DWNLT3 */
/* Subroutine */ int dwnlt3_(integer *i__, integer *imax, integer *m, integer 
	*mdw, integer *ipivot, doublereal *h__, doublereal *w)
{
    /* System generated locals */
    integer w_dim1, w_offset;

    /* Local variables */
    static doublereal t;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer itemp;

/* ***BEGIN PROLOGUE  DWNLT3 */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to WNLIT */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D) */
/* ***AUTHOR  Hanson, R. J., (SNLA) */
/*           Haskell, K. H., (SNLA) */
/* ***DESCRIPTION */

/*     Perform column interchange. */
/*     Exchange elements of permuted index vector and perform column */
/*     interchanges. */

/* ***SEE ALSO  DWNLIT */
/* ***ROUTINES CALLED  DSWAP */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790701  DATE WRITTEN */
/*   890620  Code extracted from WNLIT and made a subroutine.  (RWC)) */
/*   900604  DP version created from SP version.  (RWC) */
/* ***END PROLOGUE  DWNLT3 */



/* ***FIRST EXECUTABLE STATEMENT  DWNLT3 */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --ipivot;
    --h__;

    /* Function Body */
    if (*imax != *i__) {
	itemp = ipivot[*i__];
	ipivot[*i__] = ipivot[*imax];
	ipivot[*imax] = itemp;

	dswap_(m, &w[*imax * w_dim1 + 1], &c__1, &w[*i__ * w_dim1 + 1], &c__1)
		;

	t = h__[*imax];
	h__[*imax] = h__[*i__];
	h__[*i__] = t;
    }
    return 0;
} /* dwnlt3_ */

