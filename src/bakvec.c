/* bakvec.f -- translated by f2c (version 12.02.01).
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

/* DECK BAKVEC */
/* Subroutine */ int bakvec_(integer *nm, integer *n, real *t, real *e, 
	integer *m, real *z__, integer *ierr)
{
    /* System generated locals */
    integer t_dim1, t_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;

/* ***BEGIN PROLOGUE  BAKVEC */
/* ***PURPOSE  Form the eigenvectors of a certain real non-symmetric */
/*            tridiagonal matrix from a symmetric tridiagonal matrix */
/*            output from FIGI. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (BAKVEC-S) */
/* ***KEYWORDS  EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine forms the eigenvectors of a NONSYMMETRIC */
/*     TRIDIAGONAL matrix by back transforming those of the */
/*     corresponding symmetric matrix determined by  FIGI. */

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

/*        E contains the subdiagonal elements of the symmetric */
/*          matrix in its last N-1 positions.  E(1) is arbitrary. */
/*          E is a one-dimensional REAL array, dimensioned E(N). */

/*        M is the number of eigenvectors to be back transformed. */
/*          M is an INTEGER variable. */

/*        Z contains the eigenvectors to be back transformed */
/*          in its first M columns.  Z is a two-dimensional REAL */
/*          array, dimensioned Z(NM,M). */

/*     On OUTPUT */

/*        T is unaltered. */

/*        E is destroyed. */

/*        Z contains the transformed eigenvectors in its first M columns. */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          2*N+I      if E(I) is zero with T(I,1) or T(I-1,3) non-zero. */
/*                     In this case, the symmetric matrix is not similar */
/*                     to the original matrix, and the eigenvectors */
/*                     cannot be found by this program. */

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
/* ***END PROLOGUE  BAKVEC */


/* ***FIRST EXECUTABLE STATEMENT  BAKVEC */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    t_dim1 = *nm;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --e;

    /* Function Body */
    *ierr = 0;
    if (*m == 0) {
	goto L1001;
    }
    e[1] = 1.f;
    if (*n == 1) {
	goto L1001;
    }

    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (e[i__] != 0.f) {
	    goto L80;
	}
	if (t[i__ + t_dim1] != 0.f || t[i__ - 1 + t_dim1 * 3] != 0.f) {
	    goto L1000;
	}
	e[i__] = 1.f;
	goto L100;
L80:
	e[i__] = e[i__ - 1] * e[i__] / t[i__ - 1 + t_dim1 * 3];
L100:
	;
    }

    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {

	i__2 = *n;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    z__[i__ + j * z_dim1] *= e[i__];
/* L120: */
	}
    }

    goto L1001;
/*     .......... SET ERROR -- EIGENVECTORS CANNOT BE */
/*                FOUND BY THIS PROGRAM .......... */
L1000:
    *ierr = (*n << 1) + i__;
L1001:
    return 0;
} /* bakvec_ */

