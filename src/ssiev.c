/* ssiev.f -- translated by f2c (version 12.02.01).
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

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK SSIEV */
/* Subroutine */ int ssiev_(real *a, integer *lda, integer *n, real *e, real *
	work, integer *job, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int tred1_(integer *, integer *, real *, real *, 
	    real *, real *), tred2_(integer *, integer *, real *, real *, 
	    real *, real *), imtql2_(integer *, integer *, real *, real *, 
	    real *, integer *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), tqlrat_(integer *, real *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SSIEV */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a real symmetric matrix. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D4A1 */
/* ***TYPE      SINGLE PRECISION (SSIEV-S, CHIEV-C) */
/* ***KEYWORDS  COMPLEX HERMITIAN, EIGENVALUES, EIGENVECTORS, MATRIX, */
/*             SYMMETRIC */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/*           Moler, C. B., (U. of New Mexico) */
/*           Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     Abstract */
/*      SSIEV computes the eigenvalues and, optionally, the eigenvectors */
/*      of a real symmetric matrix. */

/*     Call Sequence Parameters- */
/*       (The values of parameters marked with * (star) will be  changed */
/*         by SSIEV.) */

/*       A*      REAL (LDA,N) */
/*               real symmetric input matrix. */
/*               Only the diagonal and upper triangle of A must be input, */
/*               as SSIEV copies the upper triangle to the lower. */
/*               That is, the user must define A(I,J), I=1,..N, and J=I,. */
/*               ..,N. */
/*               On return from SSIEV, if the user has set JOB */
/*               = 0        the lower triangle of A has been altered. */
/*               = nonzero  the N eigenvectors of A are stored in its */
/*               first N columns.  See also INFO below. */

/*       LDA     INTEGER */
/*               set by the user to */
/*               the leading dimension of the array A. */

/*       N       INTEGER */
/*               set by the user to */
/*               the order of the matrix A and */
/*               the number of elements in E. */

/*       E*      REAL (N) */
/*               on return from SSIEV, E contains the N */
/*               eigenvalues of A.  See also INFO below. */

/*       WORK*   REAL (2*N) */
/*               temporary storage vector.  Contents changed by SSIEV. */

/*       JOB     INTEGER */
/*               set by user on input */
/*               = 0         only calculate eigenvalues of A. */
/*               = nonzero   calculate eigenvalues and eigenvectors of A. */

/*       INFO*   INTEGER */
/*               on return from SSIEV, the value of INFO is */
/*               = 0 for normal return. */
/*               = K if the eigenvalue iteration fails to converge. */
/*                   eigenvalues and vectors 1 through K-1 are correct. */


/*     Error Messages- */
/*          No. 1   recoverable  N is greater than LDA */
/*          No. 2   recoverable  N is less than one */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  IMTQL2, TQLRAT, TRED1, TRED2, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800808  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  SSIEV */
/* ***FIRST EXECUTABLE STATEMENT  SSIEV */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --e;
    --work;

    /* Function Body */
    if (*n > *lda) {
	xermsg_("SLATEC", "SSIEV", "N .GT. LDA.", &c__1, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)11);
    }
    if (*n > *lda) {
	return 0;
    }
    if (*n < 1) {
	xermsg_("SLATEC", "SSIEV", "N .LT. 1", &c__2, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)8);
    }
    if (*n < 1) {
	return 0;
    }

/*       CHECK N=1 CASE */

    e[1] = a[a_dim1 + 1];
    *info = 0;
    if (*n == 1) {
	return 0;
    }

/*     COPY UPPER TRIANGLE TO LOWER */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    a[j + i__ * a_dim1] = a[i__ + j * a_dim1];
/* L10: */
	}
    }

    if (*job != 0) {
	goto L20;
    }

/*     EIGENVALUES ONLY */

    tred1_(lda, n, &a[a_offset], &e[1], &work[1], &work[*n + 1]);
    tqlrat_(n, &e[1], &work[*n + 1], info);
    return 0;

/*     EIGENVALUES AND EIGENVECTORS */

L20:
    tred2_(lda, n, &a[a_offset], &e[1], &work[1], &a[a_offset]);
    imtql2_(lda, n, &e[1], &work[1], &a[a_offset], info);
    return 0;
} /* ssiev_ */

