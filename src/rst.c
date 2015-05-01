/* rst.f -- translated by f2c (version 12.02.01).
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

/* DECK RST */
/* Subroutine */ int rst_(integer *nm, integer *n, real *w, real *e, integer *
	matz, real *z__, integer *ierr)
{
    /* System generated locals */
    integer z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int imtql1_(integer *, real *, real *, integer *),
	     imtql2_(integer *, integer *, real *, real *, real *, integer *);

/* ***BEGIN PROLOGUE  RST */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a real symmetric tridiagonal matrix. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5 */
/* ***TYPE      SINGLE PRECISION (RST-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine calls the recommended sequence of */
/*     subroutines from the eigensystem subroutine package (EISPACK) */
/*     to find the eigenvalues and eigenvectors (if desired) */
/*     of a REAL SYMMETRIC TRIDIAGONAL matrix. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, Z, as declared in the calling program */
/*          dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        W contains the diagonal elements of the real symmetric */
/*          tridiagonal matrix.  W is a one-dimensional REAL array, */
/*          dimensioned W(N). */

/*        E contains the subdiagonal elements of the matrix in its last */
/*          N-1 positions.  E(1) is arbitrary.  E is a one-dimensional */
/*          REAL array, dimensioned E(N). */

/*        MATZ is an INTEGER variable set equal to zero if only */
/*          eigenvalues are desired.  Otherwise, it is set to any */
/*          non-zero integer for both eigenvalues and eigenvectors. */

/*     On Output */

/*        W contains the eigenvalues in ascending order. */

/*        Z contains the eigenvectors if MATZ is not zero.  The eigen- */
/*          vectors are orthonormal.  Z is a two-dimensional REAL array, */
/*          dimensioned Z(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          10*N       if N is greater than NM, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after 30 iterations. */
/*                     The eigenvalues and eigenvectors in the W and Z */
/*                     arrays should be correct for indices */
/*                     1, 2, ..., IERR-1. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  IMTQL1, IMTQL2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RST */


/* ***FIRST EXECUTABLE STATEMENT  RST */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --w;
    --e;

    /* Function Body */
    if (*n <= *nm) {
	goto L10;
    }
    *ierr = *n * 10;
    goto L50;

L10:
    if (*matz != 0) {
	goto L20;
    }
/*     .......... FIND EIGENVALUES ONLY .......... */
    imtql1_(n, &w[1], &e[1], ierr);
    goto L50;
/*     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS .......... */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    z__[j + i__ * z_dim1] = 0.f;
/* L30: */
	}

	z__[i__ + i__ * z_dim1] = 1.f;
/* L40: */
    }

    imtql2_(nm, n, &w[1], &e[1], &z__[z_offset], ierr);
L50:
    return 0;
} /* rst_ */

