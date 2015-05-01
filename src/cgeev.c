/* cgeev.f -- translated by f2c (version 12.02.01).
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

/* DECK CGEEV */
/* Subroutine */ int cgeev_(real *a, integer *lda, integer *n, real *e, real *
	v, integer *ldv, real *work, integer *job, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m, ihi, ilo;
    extern /* Subroutine */ int cbal_(integer *, integer *, real *, real *, 
	    integer *, integer *, real *);
    static integer mdim;
    extern /* Subroutine */ int corth_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *), comqr_(integer *, 
	    integer *, integer *, integer *, real *, real *, real *, real *, 
	    integer *), cbabk2_(integer *, integer *, integer *, integer *, 
	    real *, integer *, real *, real *), scopy_(integer *, real *, 
	    integer *, real *, integer *), comqr2_(integer *, integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, integer *), xermsg_(char *, char *, char *
	    , integer *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  CGEEV */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a complex general matrix. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D4A4 */
/* ***TYPE      COMPLEX (SGEEV-S, CGEEV-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/*           Moler, C. B., (U. of New Mexico) */
/*           Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     Abstract */
/*      CGEEV computes the eigenvalues and, optionally, */
/*      the eigenvectors of a general complex matrix. */

/*     Call Sequence Parameters- */
/*       (The values of parameters marked with * (star) will be changed */
/*         by CGEEV.) */

/*        A*      COMPLEX(LDA,N) */
/*                complex nonsymmetric input matrix. */

/*        LDA     INTEGER */
/*                set by the user to */
/*                the leading dimension of the complex array A. */

/*        N       INTEGER */
/*                set by the user to */
/*                the order of the matrices A and V, and */
/*                the number of elements in E. */

/*        E*      COMPLEX(N) */
/*                on return from CGEEV E contains the eigenvalues of A. */
/*                See also INFO below. */

/*        V*      COMPLEX(LDV,N) */
/*                on return from CGEEV if the user has set JOB */
/*                = 0        V is not referenced. */
/*                = nonzero  the N eigenvectors of A are stored in the */
/*                first N columns of V.  See also INFO below. */
/*                (If the input matrix A is nearly degenerate, V */
/*                 will be badly conditioned, i.e. have nearly */
/*                 dependent columns.) */

/*        LDV     INTEGER */
/*                set by the user to */
/*                the leading dimension of the array V if JOB is also */
/*                set nonzero.  In that case N must be .LE. LDV. */
/*                If JOB is set to zero LDV is not referenced. */

/*        WORK*   REAL(3N) */
/*                temporary storage vector.  Contents changed by CGEEV. */

/*        JOB     INTEGER */
/*                set by the user to */
/*                = 0        eigenvalues only to be calculated by CGEEV. */
/*                           neither V nor LDV are referenced. */
/*                = nonzero  eigenvalues and vectors to be calculated. */
/*                           In this case A & V must be distinct arrays. */
/*                           Also,  if LDA > LDV,  CGEEV changes all the */
/*                           elements of A thru column N.  If LDA < LDV, */
/*                           CGEEV changes all the elements of V through */
/*                           column N.  If LDA = LDV only A(I,J) and V(I, */
/*                           J) for I,J = 1,...,N are changed by CGEEV. */

/*        INFO*   INTEGER */
/*                on return from CGEEV the value of INFO is */
/*                = 0  normal return, calculation successful. */
/*                = K  if the eigenvalue iteration fails to converge, */
/*                     eigenvalues K+1 through N are correct, but */
/*                     no eigenvectors were computed even if they were */
/*                     requested (JOB nonzero). */

/*      Error Messages */
/*           No. 1  recoverable  N is greater than LDA */
/*           No. 2  recoverable  N is less than one. */
/*           No. 3  recoverable  JOB is nonzero and N is greater than LDV */
/*           No. 4  warning      LDA > LDV,  elements of A other than the */
/*                               N by N input elements have been changed */
/*           No. 5  warning      LDA < LDV,  elements of V other than the */
/*                               N by N output elements have been changed */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  CBABK2, CBAL, COMQR, COMQR2, CORTH, SCOPY, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800808  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  CGEEV */
/* ***FIRST EXECUTABLE STATEMENT  CGEEV */
    /* Parameter adjustments */
    --work;
    --v;
    --e;
    --a;

    /* Function Body */
    if (*n > *lda) {
	xermsg_("SLATEC", "CGEEV", "N .GT. LDA.", &c__1, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)11);
    }
    if (*n > *lda) {
	return 0;
    }
    if (*n < 1) {
	xermsg_("SLATEC", "CGEEV", "N .LT. 1", &c__2, &c__1, (ftnlen)6, (
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
	xermsg_("SLATEC", "CGEEV", "JOB .NE. 0 AND N .GT. LDV.", &c__3, &c__1,
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
	xermsg_("SLATEC", "CGEEV", "LDA.LT.LDV,  ELEMENTS OF V OTHER THAN TH"
		"E N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.", &c__5, &c__0, (
		ftnlen)6, (ftnlen)5, (ftnlen)83);
    }
    if (*lda <= *ldv) {
	goto L5;
    }
    xermsg_("SLATEC", "CGEEV", "LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N B"
	    "Y N INPUT ELEMENTS HAVE BEEN CHANGED.", &c__4, &c__0, (ftnlen)6, (
	    ftnlen)5, (ftnlen)81);
    l = *n - 1;
    i__1 = l;
    for (j = 1; j <= i__1; ++j) {
	i__ = *n << 1;
	m = (j << 1) * *ldv + 1;
	k = (j << 1) * *lda + 1;
	scopy_(&i__, &a[k], &c__1, &a[m], &c__1);
/* L4: */
    }
L5:

/*     SEPARATE REAL AND IMAGINARY PARTS */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	k = (j - 1) * mdim + 1;
	l = k + *n;
	scopy_(n, &a[k + 1], &c__2, &work[1], &c__1);
	scopy_(n, &a[k], &c__2, &a[k], &c__1);
	scopy_(n, &work[1], &c__1, &a[l], &c__1);
/* L6: */
    }

/*     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG. */

    cbal_(&mdim, n, &a[1], &a[*n + 1], &ilo, &ihi, &work[1]);
    corth_(&mdim, n, &ilo, &ihi, &a[1], &a[*n + 1], &work[*n + 1], &work[(*n 
	    << 1) + 1]);
    if (*job != 0) {
	goto L10;
    }

/*     EIGENVALUES ONLY */

    comqr_(&mdim, n, &ilo, &ihi, &a[1], &a[*n + 1], &e[1], &e[*n + 1], info);
    goto L30;

/*     EIGENVALUES AND EIGENVECTORS. */

L10:
    comqr2_(&mdim, n, &ilo, &ihi, &work[*n + 1], &work[(*n << 1) + 1], &a[1], 
	    &a[*n + 1], &e[1], &e[*n + 1], &v[1], &v[*n + 1], info);
    if (*info != 0) {
	goto L30;
    }
    cbabk2_(&mdim, n, &ilo, &ihi, &work[1], n, &v[1], &v[*n + 1]);

/*     CONVERT EIGENVECTORS TO COMPLEX STORAGE. */

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

/*     CONVERT EIGENVALUES TO COMPLEX STORAGE. */

L30:
    scopy_(n, &e[1], &c__1, &work[1], &c__1);
    scopy_(n, &e[*n + 1], &c__1, &e[2], &c__2);
    scopy_(n, &work[1], &c__1, &e[1], &c__2);
    return 0;

/*     TAKE CARE OF N=1 CASE */

L35:
    e[1] = a[1];
    e[2] = a[2];
    *info = 0;
    if (*job == 0) {
	return 0;
    }
    v[1] = a[1];
    v[2] = a[2];
    return 0;
} /* cgeev_ */

