/* chiev.f -- translated by f2c (version 12.02.01).
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
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__0 = 0;
static integer c__4 = 4;
static integer c__6 = 6;

/* DECK CHIEV */
/* Subroutine */ int chiev_(real *a, integer *lda, integer *n, real *e, real *
	v, integer *ldv, real *work, integer *job, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m, mdim;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), imtql2_(integer *, integer *, real *, real *, real *, 
	    integer *), htridi_(integer *, integer *, real *, real *, real *, 
	    real *, real *, real *), htribk_(integer *, integer *, real *, 
	    real *, real *, integer *, real *, real *), xermsg_(char *, char *
	    , char *, integer *, integer *, ftnlen, ftnlen, ftnlen), tqlrat_(
	    integer *, real *, real *, integer *), scopym_(integer *, real *, 
	    integer *, real *, integer *);

/* ***BEGIN PROLOGUE  CHIEV */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a complex Hermitian matrix. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D4A3 */
/* ***TYPE      COMPLEX (SSIEV-S, CHIEV-C) */
/* ***KEYWORDS  COMPLEX HERMITIAN, EIGENVALUES, EIGENVECTORS, MATRIX, */
/*             SYMMETRIC */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/*           Moler, C. B., (U. of New Mexico) */
/*           Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     David Kahaner, Cleve Moler, G. W. Stewart, */
/*       N.B.S.         U.N.M.      N.B.S./U.MD. */

/*     Abstract */
/*      CHIEV computes the eigenvalues and, optionally, */
/*      the eigenvectors of a complex Hermitian matrix. */

/*     Call Sequence Parameters- */
/*       (the values of parameters marked with * (star) will be changed */
/*         by CHIEV.) */

/*        A*      COMPLEX(LDA,N) */
/*                complex Hermitian input matrix. */
/*                Only the upper triangle of A need be */
/*                filled in.  Elements on diagonal must be real. */

/*        LDA     INTEGER */
/*                set by the user to */
/*                the leading dimension of the complex array A. */

/*        N       INTEGER */
/*                set by the user to */
/*                the order of the matrices A and V, and */
/*                the number of elements in E. */

/*        E*      REAL(N) */
/*                on return from CHIEV E contains the eigenvalues of A. */
/*                See also INFO below. */

/*        V*      COMPLEX(LDV,N) */
/*                on return from CHIEV if the user has set JOB */
/*                = 0        V is not referenced. */
/*                = nonzero  the N eigenvectors of A are stored in the */
/*                first N columns of V.  See also INFO below. */

/*        LDV     INTEGER */
/*                set by the user to */
/*                the leading dimension of the array V if JOB is also */
/*                set nonzero.  In that case N must be .LE. LDV. */
/*                If JOB is set to zero LDV is not referenced. */

/*        WORK*   REAL(4N) */
/*                temporary storage vector.  Contents changed by CHIEV. */

/*        JOB     INTEGER */
/*                set by the user to */
/*                = 0        eigenvalues only to be calculated by CHIEV. */
/*                           Neither V nor LDV are referenced. */
/*                = nonzero  eigenvalues and vectors to be calculated. */
/*                           In this case A and V must be distinct arrays */
/*                           also if LDA .GT. LDV CHIEV changes all the */
/*                           elements of A thru column N.  If LDA < LDV */
/*                           CHIEV changes all the elements of V through */
/*                           column N.  If LDA = LDV only A(I,J) and V(I, */
/*                           J) for I,J = 1,...,N are changed by CHIEV. */

/*        INFO*   INTEGER */
/*                on return from CHIEV the value of INFO is */
/*                = 0  normal return, calculation successful. */
/*                = K  if the eigenvalue iteration fails to converge, */
/*                     eigenvalues (and eigenvectors if requested) */
/*                     1 through K-1 are correct. */

/*      Error Messages */
/*           No. 1  recoverable  N is greater than LDA */
/*           No. 2  recoverable  N is less than one. */
/*           No. 3  recoverable  JOB is nonzero and N is greater than LDV */
/*           No. 4  warning      LDA > LDV,  elements of A other than the */
/*                               N by N input elements have been changed */
/*           No. 5  warning      LDA < LDV,  elements of V other than the */
/*                               N by N output elements have been changed */
/*           No. 6  recoverable  nonreal element on diagonal of A. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  HTRIBK, HTRIDI, IMTQL2, SCOPY, SCOPYM, TQLRAT, */
/*                    XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800808  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/* ***END PROLOGUE  CHIEV */
/* ***FIRST EXECUTABLE STATEMENT  CHIEV */
    /* Parameter adjustments */
    --work;
    --v;
    --e;
    --a;

    /* Function Body */
    if (*n > *lda) {
	xermsg_("SLATEC", "CHIEV", "N .GT. LDA.", &c__1, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)11);
    }
    if (*n > *lda) {
	return 0;
    }
    if (*n < 1) {
	xermsg_("SLATEC", "CHIEV", "N .LT. 1", &c__2, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)8);
    }
    if (*n < 1) {
	return 0;
    }
    if (*n == 1 && *job == 0) {
	goto L35;
    }
    mdim = *lda << 1;
    if (*job == 0) {
	goto L5;
    }
    if (*n > *ldv) {
	xermsg_("SLATEC", "CHIEV", "JOB .NE. 0 AND N .GT. LDV.", &c__3, &c__1,
		 (ftnlen)6, (ftnlen)5, (ftnlen)26);
    }
    if (*n > *ldv) {
	return 0;
    }
    if (*n == 1) {
	goto L35;
    }

/*       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0 */

/* Computing MIN */
    i__1 = mdim, i__2 = *ldv << 1;
    mdim = min(i__1,i__2);
    if (*lda < *ldv) {
	xermsg_("SLATEC", "CHIEV", "LDA.LT.LDV,  ELEMENTS OF V OTHER THAN TH"
		"E N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.", &c__5, &c__0, (
		ftnlen)6, (ftnlen)5, (ftnlen)83);
    }
    if (*lda <= *ldv) {
	goto L5;
    }
    xermsg_("SLATEC", "CHIEV", "LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N B"
	    "Y N INPUT ELEMENTS HAVE BEEN CHANGED.", &c__4, &c__0, (ftnlen)6, (
	    ftnlen)5, (ftnlen)81);
    l = *n - 1;
    i__1 = l;
    for (j = 1; j <= i__1; ++j) {
	m = (j << 1) * *ldv + 1;
	k = (j << 1) * *lda + 1;
	i__2 = *n << 1;
	scopy_(&i__2, &a[k], &c__1, &a[m], &c__1);
/* L4: */
    }
L5:

/*     FILL IN LOWER TRIANGLE OF A, COLUMN BY COLUMN. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = (j - 1) * (mdim + 2) + 1;
	if (a[k + 1] != 0.f) {
	    xermsg_("SLATEC", "CHIEV", "NONREAL ELEMENT ON DIAGONAL OF A", &
		    c__6, &c__1, (ftnlen)6, (ftnlen)5, (ftnlen)32);
	}
	if (a[k + 1] != 0.f) {
	    return 0;
	}
	i__2 = *n - j + 1;
	scopy_(&i__2, &a[k], &mdim, &a[k], &c__2);
	i__2 = *n - j + 1;
	scopym_(&i__2, &a[k + 1], &mdim, &a[k + 1], &c__2);
/* L6: */
    }

/*     SEPARATE REAL AND IMAGINARY PARTS */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = (j - 1) * mdim + 1;
	l = k + *n;
	scopy_(n, &a[k + 1], &c__2, &work[1], &c__1);
	scopy_(n, &a[k], &c__2, &a[k], &c__1);
	scopy_(n, &work[1], &c__1, &a[l], &c__1);
/* L10: */
    }

/*    REDUCE A TO TRIDIAGONAL MATRIX. */

    htridi_(&mdim, n, &a[1], &a[*n + 1], &e[1], &work[1], &work[*n + 1], &
	    work[(*n << 1) + 1]);
    if (*job != 0) {
	goto L15;
    }

/*     EIGENVALUES ONLY. */

    tqlrat_(n, &e[1], &work[*n + 1], info);
    return 0;

/*     EIGENVALUES AND EIGENVECTORS. */

L15:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = (j - 1) * mdim + 1;
	m = k + *n - 1;
	i__2 = m;
	for (i__ = k; i__ <= i__2; ++i__) {
/* L16: */
	    v[i__] = 0.f;
	}
	i__ = k + j - 1;
	v[i__] = 1.f;
/* L17: */
    }
    imtql2_(&mdim, n, &e[1], &work[1], &v[1], info);
    if (*info != 0) {
	return 0;
    }
    htribk_(&mdim, n, &a[1], &a[*n + 1], &work[(*n << 1) + 1], n, &v[1], &v[*
	    n + 1]);

/*    CONVERT EIGENVECTORS TO COMPLEX STORAGE. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = (j - 1) * mdim + 1;
	i__ = (j - 1 << 1) * *ldv + 1;
	l = k + *n;
	scopy_(n, &v[k], &c__1, &work[1], &c__1);
	scopy_(n, &v[l], &c__1, &v[i__ + 1], &c__2);
	scopy_(n, &work[1], &c__1, &v[i__], &c__2);
/* L20: */
    }
    return 0;

/*     TAKE CARE OF N=1 CASE. */

L35:
    if (a[2] != 0.f) {
	xermsg_("SLATEC", "CHIEV", "NONREAL ELEMENT ON DIAGONAL OF A", &c__6, 
		&c__1, (ftnlen)6, (ftnlen)5, (ftnlen)32);
    }
    if (a[2] != 0.f) {
	return 0;
    }
    e[1] = a[1];
    *info = 0;
    if (*job == 0) {
	return 0;
    }
    v[1] = a[1];
    v[2] = 0.f;
    return 0;
} /* chiev_ */

