/* rgg.f -- translated by f2c (version 12.02.01).
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

/* Table of constant values */

static real c_b5 = 0.f;

/* DECK RGG */
/* Subroutine */ int rgg_(integer *nm, integer *n, real *a, real *b, real *
	alfr, real *alfi, real *beta, integer *matz, real *z__, integer *ierr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, z_dim1, z_offset;

    /* Local variables */
    static logical tf;
    extern /* Subroutine */ int qzit_(integer *, integer *, real *, real *, 
	    real *, logical *, real *, integer *), qzvec_(integer *, integer *
	    , real *, real *, real *, real *, real *, real *), qzhes_(integer 
	    *, integer *, real *, real *, logical *, real *), qzval_(integer *
	    , integer *, real *, real *, real *, real *, real *, logical *, 
	    real *);

/* ***BEGIN PROLOGUE  RGG */
/* ***PURPOSE  Compute the eigenvalues and eigenvectors for a real */
/*            generalized eigenproblem. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4B2 */
/* ***TYPE      SINGLE PRECISION (RGG-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine calls the recommended sequence of */
/*     subroutines from the eigensystem subroutine package (EISPACK) */
/*     to find the eigenvalues and eigenvectors (if desired) */
/*     for the REAL GENERAL GENERALIZED eigenproblem  Ax = (LAMBDA)Bx. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A, B, and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrices A and B.  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        A contains a real general matrix.  A is a two-dimensional */
/*          REAL array, dimensioned A(NM,N). */

/*        B contains a real general matrix.  B is a two-dimensional */
/*          REAL array, dimensioned B(NM,N). */

/*        MATZ is an INTEGER variable set equal to zero if only */
/*          eigenvalues are desired.  Otherwise, it is set to any */
/*          non-zero integer for both eigenvalues and eigenvectors. */

/*     On Output */

/*        A and B have been destroyed. */

/*        ALFR and ALFI contain the real and imaginary parts, */
/*          respectively, of the numerators of the eigenvalues. */
/*          ALFR and ALFI are one-dimensional REAL arrays, */
/*          dimensioned ALFR(N) and ALFI(N). */

/*        BETA contains the denominators of the eigenvalues, */
/*          which are thus given by the ratios  (ALFR+I*ALFI)/BETA. */
/*          Complex conjugate pairs of eigenvalues appear consecutively */
/*          with the eigenvalue having the positive imaginary part first. */
/*          BETA is a one-dimensional REAL array, dimensioned BETA(N). */

/*        Z contains the real and imaginary parts of the eigenvectors */
/*          if MATZ is not zero.  If the J-th eigenvalue is real, the */
/*          J-th column of  Z  contains its eigenvector.  If the J-th */
/*          eigenvalue is complex with positive imaginary part, the */
/*          J-th and (J+1)-th columns of  Z  contain the real and */
/*          imaginary parts of its eigenvector.  The conjugate of this */
/*          vector is the eigenvector for the conjugate eigenvalue. */
/*          Z is a two-dimensional REAL array, dimensioned Z(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          10*N       if N is greater than NM, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30*N iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     IERR+1, IERR+2, ..., N, but no eigenvectors are */
/*                     computed. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  QZHES, QZIT, QZVAL, QZVEC */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  RGG */


/* ***FIRST EXECUTABLE STATEMENT  RGG */
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
    --alfr;
    --alfi;
    --beta;

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
    tf = FALSE_;
    qzhes_(nm, n, &a[a_offset], &b[b_offset], &tf, &z__[z_offset]);
    qzit_(nm, n, &a[a_offset], &b[b_offset], &c_b5, &tf, &z__[z_offset], ierr)
	    ;
    qzval_(nm, n, &a[a_offset], &b[b_offset], &alfr[1], &alfi[1], &beta[1], &
	    tf, &z__[z_offset]);
    goto L50;
/*     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS .......... */
L20:
    tf = TRUE_;
    qzhes_(nm, n, &a[a_offset], &b[b_offset], &tf, &z__[z_offset]);
    qzit_(nm, n, &a[a_offset], &b[b_offset], &c_b5, &tf, &z__[z_offset], ierr)
	    ;
    qzval_(nm, n, &a[a_offset], &b[b_offset], &alfr[1], &alfi[1], &beta[1], &
	    tf, &z__[z_offset]);
    if (*ierr != 0) {
	goto L50;
    }
    qzvec_(nm, n, &a[a_offset], &b[b_offset], &alfr[1], &alfi[1], &beta[1], &
	    z__[z_offset]);
L50:
    return 0;
} /* rgg_ */

