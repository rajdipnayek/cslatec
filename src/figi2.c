/* figi2.f -- translated by f2c (version 12.02.01).
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

/* DECK FIGI2 */
/* Subroutine */ int figi2_(integer *nm, integer *n, real *t, real *d__, real 
	*e, real *z__, integer *ierr)
{
    /* System generated locals */
    integer t_dim1, t_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static real h__;
    static integer i__, j;

/* ***BEGIN PROLOGUE  FIGI2 */
/* ***PURPOSE  Transforms certain real non-symmetric tridiagonal matrix */
/*            to symmetric tridiagonal matrix. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C1C */
/* ***TYPE      SINGLE PRECISION (FIGI2-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     Given a NONSYMMETRIC TRIDIAGONAL matrix such that the products */
/*     of corresponding pairs of off-diagonal elements are all */
/*     non-negative, and zero only when both factors are zero, this */
/*     subroutine reduces it to a SYMMETRIC TRIDIAGONAL matrix */
/*     using and accumulating diagonal similarity transformations. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, T and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

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

/*        Z contains the diagonal transformation matrix produced in the */
/*          symmetrization.  Z is a two-dimensional REAL array, */
/*          dimensioned Z(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          N+I        if T(I,1)*T(I-1,3) is negative, */
/*          2*N+I      if T(I,1)*T(I-1,3) is zero with one factor */
/*                     non-zero.  In these cases, there does not exist */
/*                     a symmetrizing similarity transformation which */
/*                     is essential for the validity of the later */
/*                     eigenvector computation. */

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
/* ***END PROLOGUE  FIGI2 */


/* ***FIRST EXECUTABLE STATEMENT  FIGI2 */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    t_dim1 = *nm;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --d__;
    --e;

    /* Function Body */
    *ierr = 0;

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L50: */
	    z__[i__ + j * z_dim1] = 0.f;
	}

	if (i__ == 1) {
	    goto L70;
	}
	h__ = t[i__ + t_dim1] * t[i__ - 1 + t_dim1 * 3];
	if (h__ < 0.f) {
	    goto L900;
	} else if (h__ == 0) {
	    goto L60;
	} else {
	    goto L80;
	}
L60:
	if (t[i__ + t_dim1] != 0.f || t[i__ - 1 + t_dim1 * 3] != 0.f) {
	    goto L1000;
	}
	e[i__] = 0.f;
L70:
	z__[i__ + i__ * z_dim1] = 1.f;
	goto L90;
L80:
	e[i__] = sqrt(h__);
	z__[i__ + i__ * z_dim1] = z__[i__ - 1 + (i__ - 1) * z_dim1] * e[i__] /
		 t[i__ - 1 + t_dim1 * 3];
L90:
	d__[i__] = t[i__ + (t_dim1 << 1)];
/* L100: */
    }

    goto L1001;
/*     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL */
/*                ELEMENTS IS NEGATIVE .......... */
L900:
    *ierr = *n + i__;
    goto L1001;
/*     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL */
/*                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO .......... */
L1000:
    *ierr = (*n << 1) + i__;
L1001:
    return 0;
} /* figi2_ */

