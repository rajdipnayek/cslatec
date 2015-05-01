/* htrib3.f -- translated by f2c (version 12.02.01).
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

/* DECK HTRIB3 */
/* Subroutine */ int htrib3_(integer *nm, integer *n, real *a, real *tau, 
	integer *m, real *zr, real *zi)
{
    /* System generated locals */
    integer a_dim1, a_offset, zr_dim1, zr_offset, zi_dim1, zi_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    static real h__;
    static integer i__, j, k, l;
    static real s, si;

/* ***BEGIN PROLOGUE  HTRIB3 */
/* ***PURPOSE  Compute the eigenvectors of a complex Hermitian matrix from */
/*            the eigenvectors of a real symmetric tridiagonal matrix */
/*            output from HTRID3. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4C4 */
/* ***TYPE      SINGLE PRECISION (HTRIB3-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine is a translation of a complex analogue of */
/*     the ALGOL procedure TRBAK3, NUM. MATH. 11, 181-195(1968) */
/*     by Martin, Reinsch, and Wilkinson. */
/*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). */

/*     This subroutine forms the eigenvectors of a COMPLEX HERMITIAN */
/*     matrix by back transforming those of the corresponding */
/*     real symmetric tridiagonal matrix determined by  HTRID3. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A, ZR, and ZI, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        A contains some information about the unitary transformations */
/*          used in the reduction by  HTRID3.  A is a two-dimensional */
/*          REAL array, dimensioned A(NM,N). */

/*        TAU contains further information about the transformations. */
/*          TAU is a one-dimensional REAL array, dimensioned TAU(2,N). */

/*        M is the number of eigenvectors to be back transformed. */
/*          M is an INTEGER variable. */

/*        ZR contains the eigenvectors to be back transformed in its */
/*          first M columns.  The contents of ZI are immaterial.  ZR and */
/*          ZI are two-dimensional REAL arrays, dimensioned ZR(NM,M) and */
/*          ZI(NM,M). */

/*     On OUTPUT */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the transformed eigenvectors in their first M columns. */

/*     NOTE that the last component of each returned vector */
/*     is real and that vector Euclidean norms are preserved. */

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
/* ***END PROLOGUE  HTRIB3 */


/* ***FIRST EXECUTABLE STATEMENT  HTRIB3 */
    /* Parameter adjustments */
    zi_dim1 = *nm;
    zi_offset = 1 + zi_dim1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = 1 + zr_dim1;
    zr -= zr_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    tau -= 3;

    /* Function Body */
    if (*m == 0) {
	goto L200;
    }
/*     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC */
/*                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN */
/*                TRIDIAGONAL MATRIX. .......... */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {

	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    zi[k + j * zi_dim1] = -zr[k + j * zr_dim1] * tau[(k << 1) + 2];
	    zr[k + j * zr_dim1] *= tau[(k << 1) + 1];
/* L50: */
	}
    }

    if (*n == 1) {
	goto L200;
    }
/*     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES .......... */
    i__2 = *n;
    for (i__ = 2; i__ <= i__2; ++i__) {
	l = i__ - 1;
	h__ = a[i__ + i__ * a_dim1];
	if (h__ == 0.f) {
	    goto L140;
	}

	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    s = 0.f;
	    si = 0.f;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		s = s + a[i__ + k * a_dim1] * zr[k + j * zr_dim1] - a[k + i__ 
			* a_dim1] * zi[k + j * zi_dim1];
		si = si + a[i__ + k * a_dim1] * zi[k + j * zi_dim1] + a[k + 
			i__ * a_dim1] * zr[k + j * zr_dim1];
/* L110: */
	    }
/*     .......... DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW .......... */
	    s = s / h__ / h__;
	    si = si / h__ / h__;

	    i__3 = l;
	    for (k = 1; k <= i__3; ++k) {
		zr[k + j * zr_dim1] = zr[k + j * zr_dim1] - s * a[i__ + k * 
			a_dim1] - si * a[k + i__ * a_dim1];
		zi[k + j * zi_dim1] = zi[k + j * zi_dim1] - si * a[i__ + k * 
			a_dim1] + s * a[k + i__ * a_dim1];
/* L120: */
	    }

/* L130: */
	}

L140:
	;
    }

L200:
    return 0;
} /* htrib3_ */

