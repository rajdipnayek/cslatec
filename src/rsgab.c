/* rsgab.f -- translated by f2c (version 12.02.01).
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

/* DECK RSGAB */
/* Subroutine */ int rsgab_(integer *nm, integer *n, real *a, real *b, real *
	w, integer *matz, real *z__, real *fv1, real *fv2, integer *ierr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset;

    /* Local variables */
    extern /* Subroutine */ int tql2_(integer *, integer *, real *, real *, 
	    real *, integer *), tred1_(integer *, integer *, real *, real *, 
	    real *, real *), tred2_(integer *, integer *, real *, real *, 
	    real *, real *), rebak_(integer *, integer *, real *, real *, 
	    integer *, real *), reduc2_(integer *, integer *, real *, real *, 
	    real *, integer *), tqlrat_(integer *, real *, real *, integer *);

/* ***BEGIN PROLOGUE  RSGAB */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a symmetric generalized eigenproblem. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4B1 */
/* ***TYPE      SINGLE PRECISION (RSGAB-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine calls the recommended sequence of */
/*     subroutines from the eigensystem subroutine package (EISPACK) */
/*     to find the eigenvalues and eigenvectors (if desired) */
/*     for the REAL SYMMETRIC generalized eigenproblem  ABx = (LAMBDA)x. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A, B, and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrices A and B.  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        A contains a real symmetric matrix.  A is a two-dimensional */
/*          REAL array, dimensioned A(NM,N). */

/*        B contains a positive definite real symmetric matrix.  B is a */
/*          two-dimensional REAL array, dimensioned B(NM,N). */

/*        MATZ is an INTEGER variable set equal to zero if only */
/*          eigenvalues are desired.  Otherwise, it is set to any */
/*          non-zero integer for both eigenvalues and eigenvectors. */

/*     On Output */

/*        W contains the eigenvalues in ascending order.  W is a */
/*          one-dimensional REAL array, dimensioned W(N). */

/*        Z contains the eigenvectors if MATZ is not zero.  Z is a */
/*          two-dimensional REAL array, dimensioned Z(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          10*N       if N is greater than NM, */
/*          7*N+1      if B is not positive definite, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after 30 iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     1, 2, ..., IERR-1, but no eigenvectors are */
/*                     computed. */

/*        FV1 and FV2 are one-dimensional REAL arrays used for temporary */
/*          storage, dimensioned FV1(N) and FV2(N). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  REBAK, REDUC2, TQL2, TQLRAT, TRED1, TRED2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RSGAB */


/* ***FIRST EXECUTABLE STATEMENT  RSGAB */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    b_dim1 = *nm;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --w;
    --fv1;
    --fv2;

    /* Function Body */
    if (*n <= *nm) {
	goto L10;
    }
    *ierr = *n * 10;
    goto L50;

L10:
    reduc2_(nm, n, &a[a_offset], &b[b_offset], &fv2[1], ierr);
    if (*ierr != 0) {
	goto L50;
    }
    if (*matz != 0) {
	goto L20;
    }
/*     .......... FIND EIGENVALUES ONLY .......... */
    tred1_(nm, n, &a[a_offset], &w[1], &fv1[1], &fv2[1]);
    tqlrat_(n, &w[1], &fv2[1], ierr);
    goto L50;
/*     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS .......... */
L20:
    tred2_(nm, n, &a[a_offset], &w[1], &fv1[1], &z__[z_offset]);
    tql2_(nm, n, &w[1], &fv1[1], &z__[z_offset], ierr);
    if (*ierr != 0) {
	goto L50;
    }
    rebak_(nm, n, &b[b_offset], &fv2[1], n, &z__[z_offset]);
L50:
    return 0;
} /* rsgab_ */

