/* trisp.f -- translated by f2c (version 12.02.01).
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

/* DECK TRISP */
/* Subroutine */ int trisp_(integer *n, real *a, real *b, real *c__, real *
	d__, real *u, real *z__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, k;
    static real v, an, bn;
    static integer nm1, nm2;
    static real den;

/* ***BEGIN PROLOGUE  TRISP */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subsidiary to SEPELI */
/* ***LIBRARY   SLATEC */
/* ***TYPE      SINGLE PRECISION (TRISP-S) */
/* ***AUTHOR  (UNKNOWN) */
/* ***DESCRIPTION */

/*     This subroutine solves for a non-zero eigenvector corresponding */
/*     to the zero eigenvalue of the transpose of the rank */
/*     deficient ONE matrix with subdiagonal A, diagonal B, and */
/*     superdiagonal C , with A(1) in the (1,N) position, with */
/*     C(N) in the (N,1) position, and all other elements zero. */

/* ***SEE ALSO  SEPELI */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900402  Added TYPE section.  (WRB) */
/* ***END PROLOGUE  TRISP */

/* ***FIRST EXECUTABLE STATEMENT  TRISP */
    /* Parameter adjustments */
    --z__;
    --u;
    --d__;
    --c__;
    --b;
    --a;

    /* Function Body */
    bn = b[*n];
    d__[1] = a[2] / b[1];
    v = a[1];
    u[1] = c__[*n] / b[1];
    nm2 = *n - 2;
    i__1 = nm2;
    for (j = 2; j <= i__1; ++j) {
	den = b[j] - c__[j - 1] * d__[j - 1];
	d__[j] = a[j + 1] / den;
	u[j] = -c__[j - 1] * u[j - 1] / den;
	bn -= v * u[j - 1];
	v = -v * d__[j - 1];
/* L10: */
    }
    den = b[*n - 1] - c__[*n - 2] * d__[*n - 2];
    d__[*n - 1] = (a[*n] - c__[*n - 2] * u[*n - 2]) / den;
    an = c__[*n - 1] - v * d__[*n - 2];
    bn -= v * u[*n - 2];
    den = bn - an * d__[*n - 1];

/*     SET LAST COMPONENT EQUAL TO ONE */

    z__[*n] = 1.f;
    z__[*n - 1] = -d__[*n - 1];
    nm1 = *n - 1;
    i__1 = nm1;
    for (j = 2; j <= i__1; ++j) {
	k = *n - j;
	z__[k] = -d__[k] * z__[k + 1] - u[k] * z__[*n];
/* L20: */
    }
    return 0;
} /* trisp_ */

