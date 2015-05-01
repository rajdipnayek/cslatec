/* trbak3.f -- translated by f2c (version 12.02.01).
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

/* DECK TRBAK3 */
/* Subroutine */ int trbak3_(integer *nm, integer *n, integer *nv, real *a, 
	integer *m, real *z__)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static real h__;
    static integer i__, j, k, l;
    static real s;
    static integer ik, iz;

/* ***BEGIN PROLOGUE  TRBAK3 */
/* ***PURPOSE  Form the eigenvectors of a real symmetric matrix from the */
/*            eigenvectors of a symmetric tridiagonal matrix formed */
/*            by TRED3. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (TRBAK3-S) */
/* ***KEYWORDS  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of the ALGOL procedure TRBAK3, */
/*     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     This subroutine forms the eigenvectors of a REAL SYMMETRIC */
/*     matrix by back transforming those of the corresponding */
/*     symmetric tridiagonal matrix determined by  TRED3. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, Z, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        NV is an INTEGER variable set equal to the dimension of the */
/*          array A as specified in the calling program.  NV must not */
/*          be less than  N*(N+1)/2. */

/*        A contains information about the orthogonal transformations */
/*          used in the reduction by  TRED3  in its first N*(N+1)/2 */
/*          positions.  A is a one-dimensional REAL array, dimensioned */
/*          A(NV). */

/*        M is the number of columns of Z to be back transformed. */
/*          M is an INTEGER variable. */

/*        Z contains the eigenvectors to be back transformed in its */
/*          first M columns.  Z is a two-dimensional REAL array, */
/*          dimensioned Z(NM,M). */

/*     On Output */

/*        Z contains the transformed eigenvectors in its first M columns. */

/*     Note that TRBAK3 preserves vector Euclidean norms. */

/*     Questions and comments should be directed to b. s. Garbow, */
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
/* ***END PROLOGUE  TRBAK3 */


/* ***FIRST EXECUTABLE STATEMENT  TRBAK3 */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --a;

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
	iz = i__ * l / 2;
	ik = iz + i__;
	h__ = a[ik];
	if (h__ == 0.f) {
	    goto L140;
	}

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    s = 0.f;
	    ik = iz;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		++ik;
		s += a[ik] * z__[k + j * z_dim1];
/* L110: */
	    }
/*     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW .......... */
	    s = s / h__ / h__;
	    ik = iz;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		++ik;
		z__[k + j * z_dim1] -= s * a[ik];
/* L120: */
	    }

/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* trbak3_ */

