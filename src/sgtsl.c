/* sgtsl.f -- translated by f2c (version 12.02.01).
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

/* DECK SGTSL */
/* Subroutine */ int sgtsl_(integer *n, real *c__, real *d__, real *e, real *
	b, integer *info)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Local variables */
    static integer k;
    static real t;
    static integer kb, kp1, nm1, nm2;

/* ***BEGIN PROLOGUE  SGTSL */
/* ***PURPOSE  Solve a tridiagonal linear system. */
/* ***LIBRARY   SLATEC (LINPACK) */
/* ***CATEGORY  D2A2A */
/* ***TYPE      SINGLE PRECISION (SGTSL-S, DGTSL-D, CGTSL-C) */
/* ***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, TRIDIAGONAL */
/* ***AUTHOR  Dongarra, J., (ANL) */
/* ***DESCRIPTION */

/*     SGTSL given a general tridiagonal matrix and a right hand */
/*     side will find the solution. */

/*     On Entry */

/*        N       INTEGER */
/*                is the order of the tridiagonal matrix. */

/*        C       REAL(N) */
/*                is the subdiagonal of the tridiagonal matrix. */
/*                C(2) through C(N) should contain the subdiagonal. */
/*                On output, C is destroyed. */

/*        D       REAL(N) */
/*                is the diagonal of the tridiagonal matrix. */
/*                On output, D is destroyed. */

/*        E       REAL(N) */
/*                is the superdiagonal of the tridiagonal matrix. */
/*                E(1) through E(N-1) should contain the superdiagonal. */
/*                On output, E is destroyed. */

/*        B       REAL(N) */
/*                is the right hand side vector. */

/*     On Return */

/*        B       is the solution vector. */

/*        INFO    INTEGER */
/*                = 0 normal value. */
/*                = K if the K-th element of the diagonal becomes */
/*                    exactly zero.  The subroutine returns when */
/*                    this is detected. */

/* ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W. */
/*                 Stewart, LINPACK Users' Guide, SIAM, 1979. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   780814  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SGTSL */

/* ***FIRST EXECUTABLE STATEMENT  SGTSL */
    /* Parameter adjustments */
    --b;
    --e;
    --d__;
    --c__;

    /* Function Body */
    *info = 0;
    c__[1] = d__[1];
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L40;
    }
    d__[1] = e[1];
    e[1] = 0.f;
    e[*n] = 0.f;

    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*              FIND THE LARGEST OF THE TWO ROWS */

	if ((r__1 = c__[kp1], dabs(r__1)) < (r__2 = c__[k], dabs(r__2))) {
	    goto L10;
	}

/*                 INTERCHANGE ROW */

	t = c__[kp1];
	c__[kp1] = c__[k];
	c__[k] = t;
	t = d__[kp1];
	d__[kp1] = d__[k];
	d__[k] = t;
	t = e[kp1];
	e[kp1] = e[k];
	e[k] = t;
	t = b[kp1];
	b[kp1] = b[k];
	b[k] = t;
L10:

/*              ZERO ELEMENTS */

	if (c__[k] != 0.f) {
	    goto L20;
	}
	*info = k;
	goto L100;
L20:
	t = -c__[kp1] / c__[k];
	c__[kp1] = d__[kp1] + t * d__[k];
	d__[kp1] = e[kp1] + t * e[k];
	e[kp1] = 0.f;
	b[kp1] += t * b[k];
/* L30: */
    }
L40:
    if (c__[*n] != 0.f) {
	goto L50;
    }
    *info = *n;
    goto L90;
L50:

/*           BACK SOLVE */

    nm2 = *n - 2;
    b[*n] /= c__[*n];
    if (*n == 1) {
	goto L80;
    }
    b[nm1] = (b[nm1] - d__[nm1] * b[*n]) / c__[nm1];
    if (nm2 < 1) {
	goto L70;
    }
    i__1 = nm2;
    for (kb = 1; kb <= i__1; ++kb) {
	k = nm2 - kb + 1;
	b[k] = (b[k] - d__[k] * b[k + 1] - e[k] * b[k + 2]) / c__[k];
/* L60: */
    }
L70:
L80:
L90:
L100:

    return 0;
} /* sgtsl_ */

