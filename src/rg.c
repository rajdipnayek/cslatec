/* rg.f -- translated by f2c (version 12.02.01).
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

/* DECK RG */
/* Subroutine */ int rg_(integer *nm, integer *n, real *a, real *wr, real *wi,
	 integer *matz, real *z__, integer *iv1, real *fv1, integer *ierr)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset;

    /* Local variables */
    static integer is1, is2;
    extern /* Subroutine */ int hqr_(integer *, integer *, integer *, integer 
	    *, real *, real *, real *, integer *), hqr2_(integer *, integer *,
	     integer *, integer *, real *, real *, real *, real *, integer *),
	     balbak_(integer *, integer *, integer *, integer *, real *, 
	    integer *, real *), balanc_(integer *, integer *, real *, integer 
	    *, integer *, real *), elmhes_(integer *, integer *, integer *, 
	    integer *, real *, integer *), eltran_(integer *, integer *, 
	    integer *, integer *, real *, integer *, real *);

/* ***BEGIN PROLOGUE  RG */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a real general matrix. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A2 */
/* ***TYPE      SINGLE PRECISION (RG-S, CG-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine calls the recommended sequence of */
/*     subroutines from the eigensystem subroutine package (EISPACK) */
/*     To find the eigenvalues and eigenvectors (if desired) */
/*     of a REAL GENERAL matrix. */

/*     On Input */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, A and Z, as declared in the calling */
/*          program dimension statement.  NM is an INTEGER variable. */

/*        N is the order of the matrix A.  N is an INTEGER variable. */
/*          N must be less than or equal to NM. */

/*        A contains the real general matrix.  A is a two-dimensional */
/*          REAL array, dimensioned A(NM,N). */

/*        MATZ is an INTEGER variable set equal to zero if only */
/*          eigenvalues are desired.  Otherwise, it is set to any */
/*          non-zero integer for both eigenvalues and eigenvectors. */

/*     On Output */

/*        A has been destroyed. */

/*        WR and WI contain the real and imaginary parts, respectively, */
/*          of the eigenvalues.  The eigenvalues are unordered except */
/*          that complex conjugate pairs of eigenvalues appear consecu- */
/*          tively with the eigenvalue having the positive imaginary part */
/*          first.  If an error exit is made, the eigenvalues should be */
/*          correct for indices IERR+1, IERR+2, ..., N.  WR and WI are */
/*          one-dimensional REAL arrays, dimensioned WR(N) and WI(N). */

/*        Z contains the real and imaginary parts of the eigenvectors */
/*          if MATZ is not zero.  If the J-th eigenvalue is real, the */
/*          J-th column of Z contains its eigenvector.  If the J-th */
/*          eigenvalue is complex with positive imaginary part, the */
/*          J-th and (J+1)-th columns of Z contain the real and */
/*          imaginary parts of its eigenvector.  The conjugate of this */
/*          vector is the eigenvector for the conjugate eigenvalue. */
/*          Z is a two-dimensional REAL array, dimensioned Z(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          10*N       if N is greater than NM, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30 iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     IERR+1, IERR+2, ..., N, but no eigenvectors are */
/*                     computed. */

/*        IV1 and FV1 are one-dimensional temporary storage arrays of */
/*          dimension N.  IV1 is of type INTEGER and FV1 of type REAL. */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  BALANC, BALBAK, ELMHES, ELTRAN, HQR, HQR2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/*   921103  Corrected description of IV1.  (DWL, FNF and WRB) */
/* ***END PROLOGUE  RG */


/* ***FIRST EXECUTABLE STATEMENT  RG */
    /* Parameter adjustments */
    z_dim1 = *nm;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    a_dim1 = *nm;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --wr;
    --wi;
    --iv1;
    --fv1;

    /* Function Body */
    if (*n <= *nm) {
	goto L10;
    }
    *ierr = *n * 10;
    goto L50;

L10:
    balanc_(nm, n, &a[a_offset], &is1, &is2, &fv1[1]);
    elmhes_(nm, n, &is1, &is2, &a[a_offset], &iv1[1]);
    if (*matz != 0) {
	goto L20;
    }
/*     .......... FIND EIGENVALUES ONLY .......... */
    hqr_(nm, n, &is1, &is2, &a[a_offset], &wr[1], &wi[1], ierr);
    goto L50;
/*     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS .......... */
L20:
    eltran_(nm, n, &is1, &is2, &a[a_offset], &iv1[1], &z__[z_offset]);
    hqr2_(nm, n, &is1, &is2, &a[a_offset], &wr[1], &wi[1], &z__[z_offset], 
	    ierr);
    if (*ierr != 0) {
	goto L50;
    }
    balbak_(nm, n, &is1, &is2, &fv1[1], n, &z__[z_offset]);
L50:
    return 0;
} /* rg_ */

