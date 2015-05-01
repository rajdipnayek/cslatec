/* dptsl.f -- translated by f2c (version 12.02.01).
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

/* DECK DPTSL */
/* Subroutine */ int dptsl_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *b)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static doublereal t1, t2;
    static integer ke, kf, kp1, nm1, kbm1, nm1d2;

/* ***BEGIN PROLOGUE  DPTSL */
/* ***PURPOSE  Solve a positive definite tridiagonal linear system. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2B2A */
/* ***TYPE      DOUBLE PRECISION (SPTSL-S, DPTSL-D, CPTSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, POSITIVE DEFINITE, SOLVE, */
/*             TRIDIAGONAL */
/* ***AUTHOR  Dongarra, J., (ANL) */
/* ***DESCRIPTION */

/*     DPTSL, given a positive definite symmetric tridiagonal matrix and */
/*     a right hand side, will find the solution. */

/*     On Entry */

/*        N        INTEGER */
/*                 is the order of the tridiagonal matrix. */

/*        D        DOUBLE PRECISION(N) */
/*                 is the diagonal of the tridiagonal matrix. */
/*                 On output D is destroyed. */

/*        E        DOUBLE PRECISION(N) */
/*                 is the offdiagonal of the tridiagonal matrix. */
/*                 E(1) through E(N-1) should contain the */
/*                 offdiagonal. */

/*        B        DOUBLE PRECISION(N) */
/*                 is the right hand side vector. */

/*     On Return */

/*        B        contains the solution. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890505  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DPTSL */


/*     CHECK FOR 1 X 1 CASE */

/* ***FIRST EXECUTABLE STATEMENT  DPTSL */
    /* Parameter adjustments */
    --b;
    --e;
    --d__;

    /* Function Body */
    if (*n != 1) {
	goto L10;
    }
    b[1] /= d__[1];
    goto L70;
L10:
    nm1 = *n - 1;
    nm1d2 = nm1 / 2;
    if (*n == 2) {
	goto L30;
    }
    kbm1 = *n - 1;

/*           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF */
/*           SUPERDIAGONAL */

    i__1 = nm1d2;
    for (k = 1; k <= i__1; ++k) {
	t1 = e[k] / d__[k];
	d__[k + 1] -= t1 * e[k];
	b[k + 1] -= t1 * b[k];
	t2 = e[kbm1] / d__[kbm1 + 1];
	d__[kbm1] -= t2 * e[kbm1];
	b[kbm1] -= t2 * b[kbm1 + 1];
	--kbm1;
/* L20: */
    }
L30:
    kp1 = nm1d2 + 1;

/*        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER */

    if (*n % 2 != 0) {
	goto L40;
    }
    t1 = e[kp1] / d__[kp1];
    d__[kp1 + 1] -= t1 * e[kp1];
    b[kp1 + 1] -= t1 * b[kp1];
    ++kp1;
L40:

/*        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP */
/*        AND BOTTOM */

    b[kp1] /= d__[kp1];
    if (*n == 2) {
	goto L60;
    }
    k = kp1 - 1;
    ke = kp1 + nm1d2 - 1;
    i__1 = ke;
    for (kf = kp1; kf <= i__1; ++kf) {
	b[k] = (b[k] - e[k] * b[k + 1]) / d__[k];
	b[kf + 1] = (b[kf + 1] - e[kf] * b[kf]) / d__[kf + 1];
	--k;
/* L50: */
    }
L60:
    if (*n % 2 == 0) {
	b[1] = (b[1] - e[1] * b[2]) / d__[1];
    }
L70:
    return 0;
} /* dptsl_ */

