/* ch.f -- translated by f2c (version 12.02.01).
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

/* DECK CH */
/* Subroutine */ int ch_(integer *nm, integer *n, real *ar, real *ai, real *w,
	 integer *matz, real *zr, real *zi, real *fv1, real *fv2, real *fm1, 
	integer *ierr)
{
    /* System generated locals */
    integer ar_dim1, ar_offset, ai_dim1, ai_offset, zr_dim1, zr_offset, 
	    zi_dim1, zi_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int tql2_(integer *, integer *, real *, real *, 
	    real *, integer *), htridi_(integer *, integer *, real *, real *, 
	    real *, real *, real *, real *), htribk_(integer *, integer *, 
	    real *, real *, real *, integer *, real *, real *), tqlrat_(
	    integer *, real *, real *, integer *);

/* ***BEGIN PROLOGUE  CH */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a complex Hermitian matrix. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A3 */
/* ***TYPE      COMPLEX (RS-S, CH-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK */
/* ***AUTHOR  Smith, B. T., et al. */
/* ***DESCRIPTION */

/*     This subroutine calls the recommended sequence of */
/*     subroutines from the eigensystem subroutine package (EISPACK) */
/*     to find the eigenvalues and eigenvectors (if desired) */
/*     of a COMPLEX HERMITIAN matrix. */

/*     On INPUT */

/*        NM must be set to the row dimension of the two-dimensional */
/*          array parameters, AR, AI, ZR and ZI, as declared in the */
/*          calling program dimension statement.  NM is an INTEGER */
/*          variable. */

/*        N is the order of the matrix A=(AR,AI).  N is an INTEGER */
/*          variable.  N must be less than or equal to NM. */

/*        AR and AI contain the real and imaginary parts, respectively, */
/*          of the complex Hermitian matrix.  AR and AI are */
/*          two-dimensional REAL arrays, dimensioned AR(NM,N) */
/*          and AI(NM,N). */

/*        MATZ is an INTEGER variable set equal to zero if only */
/*          eigenvalues are desired.  Otherwise, it is set to any */
/*          non-zero integer for both eigenvalues and eigenvectors. */

/*     On OUTPUT */

/*        W contains the eigenvalues in ascending order. */
/*          W is a one-dimensional REAL array, dimensioned W(N). */

/*        ZR and ZI contain the real and imaginary parts, respectively, */
/*          of the eigenvectors if MATZ is not zero.  ZR and ZI are */
/*          two-dimensional REAL arrays, dimensioned ZR(NM,N) and */
/*          ZI(NM,N). */

/*        IERR is an INTEGER flag set to */
/*          Zero       for normal return, */
/*          10*N       if N is greater than NM, */
/*          J          if the J-th eigenvalue has not been */
/*                     determined after a total of 30 iterations. */
/*                     The eigenvalues should be correct for indices */
/*                     1, 2, ..., IERR-1, but no eigenvectors are */
/*                     computed. */

/*        FV1 and FV2 are one-dimensional REAL arrays used for */
/*          temporary storage, dimensioned FV1(N) and FV2(N). */

/*        FM1 is a two-dimensional REAL array used for temporary */
/*          storage, dimensioned FM1(2,N). */

/*     Questions and comments should be directed to B. S. Garbow, */
/*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY */
/*     ------------------------------------------------------------------ */

/* ***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow, */
/*                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen- */
/*                 system Routines - EISPACK Guide, Springer-Verlag, */
/*                 1976. */
/* ***ROUTINES CALLED  HTRIBK, HTRIDI, TQL2, TQLRAT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   760101  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  CH */


/* ***FIRST EXECUTABLE STATEMENT  CH */
    /* Parameter adjustments */
    zi_dim1 = *nm;
    zi_offset = 1 + zi_dim1;
    zi -= zi_offset;
    zr_dim1 = *nm;
    zr_offset = 1 + zr_dim1;
    zr -= zr_offset;
    ai_dim1 = *nm;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    ar_dim1 = *nm;
    ar_offset = 1 + ar_dim1;
    ar -= ar_offset;
    --w;
    --fv1;
    --fv2;
    fm1 -= 3;

    /* Function Body */
    if (*n <= *nm) {
	goto L10;
    }
    *ierr = *n * 10;
    goto L50;

L10:
    htridi_(nm, n, &ar[ar_offset], &ai[ai_offset], &w[1], &fv1[1], &fv2[1], &
	    fm1[3]);
    if (*matz != 0) {
	goto L20;
    }
/*     .......... FIND EIGENVALUES ONLY .......... */
    tqlrat_(n, &w[1], &fv2[1], ierr);
    goto L50;
/*     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS .......... */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {

	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    zr[j + i__ * zr_dim1] = 0.f;
/* L30: */
	}

	zr[i__ + i__ * zr_dim1] = 1.f;
/* L40: */
    }

    tql2_(nm, n, &w[1], &fv1[1], &zr[zr_offset], ierr);
    if (*ierr != 0) {
	goto L50;
    }
    htribk_(nm, n, &ar[ar_offset], &ai[ai_offset], &fm1[3], n, &zr[zr_offset],
	     &zi[zi_offset]);
L50:
    return 0;
} /* ch_ */

