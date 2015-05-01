/* figi.f -- translated by f2c (version 12.02.01).
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

/* DECK FIGI */
/* Subroutine */ int figi_(integer *nm, integer *n, real *t, real *d__, real *
	e, real *e2, integer *ierr)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1;
    real r__1;

    /* Local variables */
    static integer i__;

/* ***BEGIN PROLOGUE  FIGI */
/* ***PURPOSE  Transforms certain real non-symmetric tridiagonal matrix */
/*            to symmetric tridiagonal matrix. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1C */
/* ***TYPE      SINGLE PRECISION (FIGI-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     Given a NONSYMMETRIC TRIDIAGONAL matrix such that the products */
/*     of corresponding pairs of off-diagonal elements are all */
/*     non-negative, this subroutine reduces it to a symmetric */
/*     tridiagonal matrix with the same eigenvalues.  If, further, */
/*     a zero product only occurs when both factors are zero, */
/*     the reduced matrix is similar to the original matrix. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, T, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix T.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        T contains the nonsymmetric matrix.  Its subdiagonal is */
/*          stored in the last N-1 positions of the first column, */
/*          its diagonal in the N positions of the second column, */
/*          and its superdiagonal in the first N-1 positions of */
/*          the third column.  T(1,1) and T(N,3) are arbitrary. */
/*          T is a two-dimensional REAL array, dimensioned T(NM,3). */

/*     On OUTPUT */

/*        T is unaltered. */

/*        D contains the diagonal elements of the tridiagonal symmetric */
/*          matrix.  D is a one-dimensional REAL array, dimensioned D(N). */

/*        E contains the subdiagonal elements of the tridiagonal */
/*          symmetric matrix in its last N-1 positions.  E(1) is not set. */
/*          E is a one-dimensional REAL array, dimensioned E(N). */

/*        E2 contains the squares of the corresponding elements of E. */
/*          E2 may coincide with E if the squares are not needed. */
/*          E2 is a one-dimensional REAL array, dimensioned E2(N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          N+I        if T(I,1)*T(I-1,3) is negative and a symmetric */
/*                     matrix cannot be produced with FIGI, */
/*          -(3*N+I)   if T(I,1)*T(I-1,3) is zero with one factor */
/*                     non-zero.  In this case, the eigenvectors of */
/*                     the symmetric matrix are not simply related */
/*                     to those of  T  and should not be sought. */

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
/* ***END PROLOGUE  FIGI */


/* ***FIRST EXECUTABLE STATEMENT  FIGI */
    /* Parameter adjustments */
    t_dim1 = *nm;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --d__;
    --e;
    --e2;

    /* Function Body */
    *ierr = 0;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ == 1) {
	    goto L90;
	}
	e2[i__] = t[i__ + t_dim1] * t[i__ - 1 + t_dim1 * 3];
	if ((r__1 = e2[i__]) < 0.f) {
	    goto L1000;
	} else if (r__1 == 0) {
	    goto L60;
	} else {
	    goto L80;
	}
L60:
	if (t[i__ + t_dim1] == 0.f && t[i__ - 1 + t_dim1 * 3] == 0.f) {
	    goto L80;
	}
/*     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL */
/*                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO .......... */
	*ierr = -(*n * 3 + i__);
L80:
	e[i__] = sqrt(e2[i__]);
L90:
	d__[i__] = t[i__ + (t_dim1 << 1)];
/* L100: */
    }

    goto L1001;
/*     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL */
/*                ELEMENTS IS NEGATIVE .......... */
L1000:
    *ierr = *n + i__;
L1001:
    return 0;
} /* figi_ */

