/* dbksol.f -- translated by f2c (version 12.02.01).
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

/* DECK DBKSOL */
/* Subroutine */ int dbksol_(integer *n, doublereal *a, doublereal *x)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k, m, nm1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);

/* ***BEGIN PROLOGUE  DBKSOL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to DBVSUP */
/* ***LIBRARY   SLATEC */
/* ***TYPE      DOUBLE PRECISION (BKSOL-S, DBKSOL-D) */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/* ********************************************************************** */
/*     Solution of an upper triangular linear system by */
/*     back-substitution */

/*     The matrix A is assumed to be stored in a linear */
/*     array proceeding in a row-wise manner. The */
/*     vector X contains the given constant vector on input */
/*     and contains the solution on return. */
/*     The actual diagonal of A is unity while a diagonal */
/*     scaling matrix is stored there. */
/* ********************************************************************** */

/* ***SEE ALSO  DBVSUP */
/* ***ROUTINES CALLED  DDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   750601  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/*   910722  Updated AUTHOR section.  (ALS) */
/* ***END PROLOGUE  DBKSOL */


/* ***FIRST EXECUTABLE STATEMENT  DBKSOL */
    /* Parameter adjustments */
    --x;
    --a;

    /* Function Body */
    m = *n * (*n + 1) / 2;
    x[*n] *= a[m];
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L20;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	j = *n - k;
	m = m - k - 1;
	x[j] = x[j] * a[m] - ddot_(&k, &a[m + 1], &c__1, &x[j + 1], &c__1);
/* L10: */
    }
L20:

    return 0;
} /* dbksol_ */

