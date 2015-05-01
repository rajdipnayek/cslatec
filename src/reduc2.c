/* reduc2.f -- translated by f2c (version 12.02.01).
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

/* DECK REDUC2 */
/* Subroutine */ int reduc2_(integer *nm, integer *n, real *a, real *b, real *
	dl, integer *ierr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k;
    static real x, y;
    static integer i1, j1, nn;

/* ***BEGIN PROLOGUE  REDUC2 */
/* ***PURPOSE  Reduce a certain generalized symmetric eigenproblem to a */
/*            standard symmetric eigenproblem using Cholesky */
/*            factorization. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1C */
/* ***TYPE      SINGLE PRECISION (REDUC2-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure REDUC2, */
/*     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971). */

/*     This subroutine reduces the generalized SYMMETRIC eigenproblems */
/*     ABx=(LAMBDA)x OR BAy=(LAMBDA)y, where B is POSITIVE DEFINITE, */
/*     to the standard symmetric eigenproblem using the Cholesky */
/*     factorization of B. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and B, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrices A and B.  If the Cholesky */
/*          factor L of B is already available, N should be prefixed */
/*          with a minus sign.  N is an INTEGER variable. */

/*        A and B contain the real symmetric input matrices.  Only */
/*          the full upper triangles of the matrices need be supplied. */
/*          If N is negative, the strict lower triangle of B contains, */
/*          instead, the strict lower triangle of its Cholesky factor L. */
/*          A and B are two-dimensional REAL arrays, dimensioned A(NM,N) */
/*          and B(NM,N). */

/*       DL contains, if N is negative, the diagonal elements of L. */
/*          DL is a one-dimensional REAL array, dimensioned DL(N). */

/*     On Output */

/*        A contains in its full lower triangle the full lower triangle */
/*          of the symmetric matrix derived from the reduction to the */
/*          standard form.  The strict upper triangle of A is unaltered. */

/*        B contains in its strict lower triangle the strict lower */
/*          triangle of its Cholesky factor L.  The full upper triangle */
/*          of B is unaltered. */

/*        DL contains the diagonal elements of L. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          7*N+1      if B is not positive definite. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  REDUC2 */


/* ***FIRST EXECUTABLE STATEMENT  REDUC2 */
    /* Parameter adjustments */
    b_dim1 = *nm;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --dl;

    /* Function Body */
    *ierr = 0;
    nn = abs(*n);
    if (*n < 0) {
	goto L100;
    }
/*     .......... FORM L IN THE ARRAYS B AND DL .......... */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = i__ - 1;

	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    x = b[i__ + j * b_dim1];
	    if (i__ == 1) {
		goto L40;
	    }

	    i__3 = i1;
	    for (k = 1; k <= i__3; ++k) {
/* L20: */
		x -= b[i__ + k * b_dim1] * b[j + k * b_dim1];
	    }

L40:
	    if (j != i__) {
		goto L60;
	    }
	    if (x <= 0.f) {
		goto L1000;
	    }
	    y = sqrt(x);
	    dl[i__] = y;
	    goto L80;
L60:
	    b[j + i__ * b_dim1] = x / y;
L80:
	    ;
	}
    }
/*     .......... FORM THE LOWER TRIANGLE OF A*L */
/*                IN THE LOWER TRIANGLE OF THE ARRAY A .......... */
L100:
    i__2 = nn;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i1 = i__ + 1;

	i__1 = i__;
	for (j = 1; j <= i__1; ++j) {
	    x = a[j + i__ * a_dim1] * dl[j];
	    if (j == i__) {
		goto L140;
	    }
	    j1 = j + 1;

	    i__3 = i__;
	    for (k = j1; k <= i__3; ++k) {
/* L120: */
		x += a[k + i__ * a_dim1] * b[k + j * b_dim1];
	    }

L140:
	    if (i__ == nn) {
		goto L180;
	    }

	    i__3 = nn;
	    for (k = i1; k <= i__3; ++k) {
/* L160: */
		x += a[i__ + k * a_dim1] * b[k + j * b_dim1];
	    }

L180:
	    a[i__ + j * a_dim1] = x;
/* L200: */
	}
    }
/*     .......... PRE-MULTIPLY BY TRANSPOSE(L) AND OVERWRITE .......... */
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = i__ + 1;
	y = dl[i__];

	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    x = y * a[i__ + j * a_dim1];
	    if (i__ == nn) {
		goto L280;
	    }

	    i__3 = nn;
	    for (k = i1; k <= i__3; ++k) {
/* L260: */
		x += a[k + j * a_dim1] * b[k + i__ * b_dim1];
	    }

L280:
	    a[i__ + j * a_dim1] = x;
/* L300: */
	}
    }

    goto L1001;
/*     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE .......... */
L1000:
    *ierr = *n * 7 + 1;
L1001:
    return 0;
} /* reduc2_ */

