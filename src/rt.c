/* rt.f -- translated by f2c (version 12.02.01).
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

/* DECK RT */
/* Subroutine */ int rt_(integer *nm, integer *n, real *a, real *w, integer *
	matz, real *z__, real *fv1, integer *ierr)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset;

    /* Local variables */
    extern /* Subroutine */ int figi_(integer *, integer *, real *, real *, 
	    real *, real *, integer *), figi2_(integer *, integer *, real *, 
	    real *, real *, real *, integer *), imtql1_(integer *, real *, 
	    real *, integer *), imtql2_(integer *, integer *, real *, real *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  RT */
/* ***PURPOSE  Compute the eigenvalues and eigenvectors of a special real */
/*            tridiagonal matrix. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A5 */
/* ***TYPE      SINGLE PRECISION (RT-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine calls the recommended sequence of subroutines */
/*     from the eigensystem subroutine package (EISPACK) to find the */
/*     eigenvalues and eigenvectors (if desired) of a special REAL */
/*     TRIDIAGONAL matrix.  The property of the matrix required for use */
/*     of this subroutine is that the products of pairs of corresponding */
/*     off-diagonal elements be all non-negative.  If eigenvectors are */
/*     desired, no product can be zero unless both factors are zero. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameter, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        A contains the special real tridiagonal matrix in its first */
/*          three columns.  The subdiagonal elements are stored in the */
/*          last N-1 positions of the first column, the diagonal elements */
/*          in the second column, and the superdiagonal elements in the */
/*          first N-1 positions of the third column.  Elements A(1,1) and */
/*          A(N,3) are arbitrary.  A is a two-dimensional REAL array, */
/*          dimensioned A(NM,3). */

/*        MATZ is an INTEGER variable set equal to zero if only */
/*          eigenvalues are desired.  Otherwise, it is set to any */
/*          non-zero integer for both eigenvalues and eigenvectors. */

/*     On Output */

/*        W contains the eigenvalues in ascending order.  W is a */
/*          one-dimensional REAL array, dimensioned W(N). */

/*        Z contains the eigenvectors if MATZ is not zero.  The eigen- */
/*          vectors are not normalized.  Z is a two-dimensional REAL */
/*          array, dimensioned Z(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          10*N       if N is greater than NM, */
/*          N+J        if A(J,1)*A(J-1,3) is negative, */
/*          2*N+J      if the product is zero with one factor non-zero, */
/*                     and MATZ is non-zero; */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after 30 iterations. */
/*                     The eigenvalues and eigenvectors in the W and Z */
/*                     arrays should be correct for indices */
/*                     1, 2, ..., IERR-1. */

/*        FV1 is a one-dimensional REAL array used for temporary storage, */
/*          dimensioned FV1(N). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  FIGI, FIGI2, IMTQL1, IMTQL2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RT */


/* ***FIRST EXECUTABLE STATEMENT  RT */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --fv1;

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
    figi_(nm, n, &a[a_offset], &w[1], &fv1[1], &fv1[1], ierr);
    if (*ierr > 0) {
	goto L50;
    }
    imtql1_(n, &w[1], &fv1[1], ierr);
    goto L50;
/*     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS .......... */
L20:
    figi2_(nm, n, &a[a_offset], &w[1], &fv1[1], &z__[z_offset], ierr);
    if (*ierr != 0) {
	goto L50;
    }
    imtql2_(nm, n, &w[1], &fv1[1], &z__[z_offset], ierr);
L50:
    return 0;
} /* rt_ */

