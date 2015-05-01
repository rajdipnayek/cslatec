/* trbak1.f -- translated by f2c (version 12.02.01).
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

/* DECK TRBAK1 */
/* Subroutine */ int trbak1_(integer *nm, integer *n, real *a, real *e, 
	integer *m, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l;
    static real s;

/* ***BEGIN PROLOGUE  TRBAK1 */
/* ***PURPOSE  Form the eigenvectors of real symmetric matrix from */
/*            the eigenvectors of a symmetric tridiagonal matrix formed */
/*            by TRED1. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (TRBAK1-S) */
/* ***KEYWORDS  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure TRBAK1, */
/*     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     This subroutine forms the eigenvectors of a REAL SYMMETRIC */
/*     matrix by back transforming those of the corresponding */
/*     symmetric tridiagonal matrix determined by  TRED1. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        A contains information about the orthogonal transformations */
/*          used in the reduction by  TRED1  in its strict lower */
/*          triangle.  A is a two-dimensional REAL array, dimensioned */
/*          A(NM,N). */

/*        E contains the subdiagonal elements of the tridiagonal matrix */
/*          in its last N-1 positions.  E(1) is arbitrary.  These */
/*          elements provide the remaining information about the */
/*          orthogonal transformations.  E is a one-dimensional REAL */
/*          array, dimensioned E(N). */

/*        M is the number of columns of Z to be back transformed. */
/*          M is an INTEGER variable. */

/*        Z contains the eigenvectors to be back transformed in its */
/*          first M columns.  Z is a two-dimensional REAL array, */
/*          dimensioned Z(NM,M). */

/*     On Output */

/*        Z contains the transformed eigenvectors in its first M columns. */

/*     Note that TRBAK1 preserves vector Euclidean norms. */

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
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  TRBAK1 */


/* ***FIRST EXECUTABLE STATEMENT  TRBAK1 */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;

    /* Function Body */
    if (*m == 0) {
	goto L200;
    }
    if (*n == 1) {
	goto L200;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	l = i__ - 1;
	if (e[i__] == 0.f) {
	    goto L140;
	}

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    s = 0.f;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
/* L110: */
		s += a[i__ + k * a_dim1] * z__[k + j * z_dim1];
	    }
/*     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1. */
/*                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW .......... */
	    s = s / a[i__ + l * a_dim1] / e[i__];

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
/* L120: */
		z__[k + j * z_dim1] += s * a[i__ + k * a_dim1];
	    }

/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* trbak1_ */

