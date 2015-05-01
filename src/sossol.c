/* sossol.f -- translated by f2c (version 12.02.01).
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

/* DECK SOSSOL */
/* Subroutine */ int sossol_(integer *k, integer *n, integer *l, real *x, 
	real *c__, real *b, integer *m)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, kj, lk, km, kn, km1, np1, jkm, kmm1;
    static real xmax;

/* ***BEGIN PROLOGUE  SOSSOL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SOS */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (SOSSOL-S, DSOSSL-D) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     SOSSOL solves an upper triangular type of linear system by back */
/*     substitution. */

/*     The matrix C is upper trapezoidal and stored as a linear array by */
/*     rows. The equations have been normalized so that the diagonal */
/*     entries of C are understood to be unity. The off diagonal entries */
/*     and the elements of the constant right hand side vector B have */
/*     already been stored as the negatives of the corresponding equation */
/*     values. */
/*     with each call to SOSSOL a (K-1) by (K-1) triangular system is */
/*     resolved. For L greater than K, column L of C is included in the */
/*     right hand side vector. */

/* ***SEE ALSO  SOS */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900328  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  SOSSOL */



/* ***FIRST EXECUTABLE STATEMENT  SOSSOL */
    /* Parameter adjustments */
    --b;
    --c__;
    --x;

    /* Function Body */
    np1 = *n + 1;
    km1 = *k - 1;
    lk = km1;
    if (*l == *k) {
	lk = *k;
    }
    kn = *m;


    i__1 = km1;
    for (kj = 1; kj <= i__1; ++kj) {
	kmm1 = *k - kj;
	km = kmm1 + 1;
	xmax = 0.f;
	kn = kn - np1 + kmm1;
	if (km > lk) {
	    goto L20;
	}
	jkm = kn;

	i__2 = lk;
	for (j = km; j <= i__2; ++j) {
	    ++jkm;
	    xmax += c__[jkm] * x[j];
/* L10: */
	}

L20:
	if (*l <= *k) {
	    goto L30;
	}
	jkm = kn + *l - kmm1;
	xmax += c__[jkm] * x[*l];
L30:
	x[kmm1] = xmax + b[kmm1];
/* L40: */
    }

    return 0;
} /* sossol_ */

