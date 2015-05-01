/* sgeev.f -- translated by f2c (version 12.02.01).
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
static real c_b39 = 0.f;

/* DECK SGEEV */
/* Subroutine */ int sgeev_(real *a, integer *lda, integer *n, real *e, real *
	v, integer *ldv, real *work, integer *job, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, m, jb, km, kp, ihi, ilo;
    extern /* Subroutine */ int hqr_(integer *, integer *, integer *, integer 
	    *, real *, real *, real *, integer *), hqr2_(integer *, integer *,
	     integer *, integer *, real *, real *, real *, real *, integer *);
    static integer mdim;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), balbak_(integer *, integer *, integer *, integer *, 
	    real *, integer *, real *), balanc_(integer *, integer *, real *, 
	    integer *, integer *, real *), orthes_(integer *, integer *, 
	    integer *, integer *, real *, real *), xermsg_(char *, char *, 
	    char *, integer *, integer *, ftnlen, ftnlen, ftnlen), ortran_(
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    ), scopym_(integer *, real *, integer *, real *, integer *);

/* ***BEGIN PROLOGUE  SGEEV */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a real general matrix. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  D4A2 */
/* ***TYPE      SINGLE PRECISION (SGEEV-S, CGEEV-C) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, GENERAL MATRIX */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/*           Moler, C. B., (U. of New Mexico) */
/*           Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     Abstract */
/*      SGEEV computes the eigenvalues and, optionally, */
/*      the eigenvectors of a general real matrix. */

/*     Call Sequence Parameters- */
/*       (The values of parameters marked with * (star) will be changed */
/*         by SGEEV.) */

/*        A*      REAL(LDA,N) */
/*                real nonsymmetric input matrix. */

/*        LDA     INTEGER */
/*                set by the user to */
/*                the leading dimension of the real array A. */

/*        N       INTEGER */
/*                set by the user to */
/*                the order of the matrices A and V, and */
/*                the number of elements in E. */

/*        E*      COMPLEX(N) */
/*                on return from SGEEV, E contains the eigenvalues of A. */
/*                See also INFO below. */

/*        V*      COMPLEX(LDV,N) */
/*                on return from SGEEV, if the user has set JOB */
/*                = 0        V is not referenced. */
/*                = nonzero  the N eigenvectors of A are stored in the */
/*                first N columns of V.  See also INFO below. */
/*                (Note that if the input matrix A is nearly degenerate, */
/*                 V may be badly conditioned, i.e., may have nearly */
/*                 dependent columns.) */

/*        LDV     INTEGER */
/*                set by the user to */
/*                the leading dimension of the array V if JOB is also */
/*                set nonzero.  In that case, N must be .LE. LDV. */
/*                If JOB is set to zero, LDV is not referenced. */

/*        WORK*   REAL(2N) */
/*                temporary storage vector.  Contents changed by SGEEV. */

/*        JOB     INTEGER */
/*                set by the user to */
/*                = 0        eigenvalues only to be calculated by SGEEV. */
/*                           Neither V nor LDV is referenced. */
/*                = nonzero  eigenvalues and vectors to be calculated. */
/*                           In this case, A & V must be distinct arrays. */
/*                           Also, if LDA .GT. LDV, SGEEV changes all the */
/*                           elements of A thru column N.  If LDA < LDV, */
/*                           SGEEV changes all the elements of V through */
/*                           column N. If LDA = LDV, only A(I,J) and V(I, */
/*                           J) for I,J = 1,...,N are changed by SGEEV. */

/*        INFO*   INTEGER */
/*                on return from SGEEV the value of INFO is */
/*                = 0  normal return, calculation successful. */
/*                = K  if the eigenvalue iteration fails to converge, */
/*                     eigenvalues K+1 through N are correct, but */
/*                     no eigenvectors were computed even if they were */
/*                     requested (JOB nonzero). */

/*      Error Messages */
/*           No. 1  recoverable  N is greater than LDA */
/*           No. 2  recoverable  N is less than one. */
/*           No. 3  recoverable  JOB is nonzero and N is greater than LDV */
/*           No. 4  warning      LDA > LDV, elements of A other than the */
/*                               N by N input elements have been changed. */
/*           No. 5  warning      LDA < LDV, elements of V other than the */
/*                               N x N output elements have been changed. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  BALANC, BALBAK, HQR, HQR2, ORTHES, ORTRAN, SCOPY, */
/*                    SCOPYM, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800808  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890531  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  SGEEV */
/* ***FIRST EXECUTABLE STATEMENT  SGEEV */
    /* Parameter adjustments */
    --work;
    --v;
    --e;
    --a;

    /* Function Body */
    if (*n > *lda) {
	xermsg_("SLATEC", "SGEEV", "N .GT. LDA.", &c__1, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)11);
    }
    if (*n > *lda) {
	return 0;
    }
    if (*n < 1) {
	xermsg_("SLATEC", "SGEEV", "N .LT. 1", &c__2, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)8);
    }
    if (*n < 1) {
	return 0;
    }
    if (*n == 1 && *job == 0) {
	goto L35;
    }
    mdim = *lda;
    if (*job == 0) {
	goto L5;
    }
    if (*n > *ldv) {
	xermsg_("SLATEC", "SGEEV", "JOB .NE. 0 AND N .GT. LDV.", &c__3, &c__1,
		 (ftnlen)6, (ftnlen)5, (ftnlen)26);
    }
    if (*n > *ldv) {
	return 0;
    }
    if (*n == 1) {
	goto L35;
    }

/*       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0 */

    mdim = min(*lda,*ldv);
    if (*lda < *ldv) {
	xermsg_("SLATEC", "SGEEV", "LDA.LT.LDV,  ELEMENTS OF V OTHER THAN TH"
		"E N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.", &c__5, &c__0, (
		ftnlen)6, (ftnlen)5, (ftnlen)83);
    }
    if (*lda <= *ldv) {
	goto L5;
    }
    xermsg_("SLATEC", "SGEEV", "LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N B"
	    "Y N INPUT ELEMENTS HAVE BEEN CHANGED.", &c__4, &c__0, (ftnlen)6, (
	    ftnlen)5, (ftnlen)81);
    l = *n - 1;
    i__1 = l;
    for (j = 1; j <= i__1; ++j) {
	m = j * *ldv + 1;
	k = j * *lda + 1;
	scopy_(n, &a[k], &c__1, &a[m], &c__1);
/* L4: */
    }
L5:

/*     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG. */

    balanc_(&mdim, n, &a[1], &ilo, &ihi, &work[1]);
    orthes_(&mdim, n, &ilo, &ihi, &a[1], &work[*n + 1]);
    if (*job != 0) {
	goto L10;
    }

/*     EIGENVALUES ONLY */

    hqr_(lda, n, &ilo, &ihi, &a[1], &e[1], &e[*n + 1], info);
    goto L30;

/*     EIGENVALUES AND EIGENVECTORS. */

L10:
    ortran_(&mdim, n, &ilo, &ihi, &a[1], &work[*n + 1], &v[1]);
    hqr2_(&mdim, n, &ilo, &ihi, &a[1], &e[1], &e[*n + 1], &v[1], info);
    if (*info != 0) {
	goto L30;
    }
    balbak_(&mdim, n, &ilo, &ihi, &work[1], n, &v[1]);

/*     CONVERT EIGENVECTORS TO COMPLEX STORAGE. */

    i__1 = *n;
    for (jb = 1; jb <= i__1; ++jb) {
	j = *n + 1 - jb;
	i__ = *n + j;
	k = (j - 1) * mdim + 1;
	kp = k + mdim;
	km = k - mdim;
	if (e[i__] >= 0.f) {
	    scopy_(n, &v[k], &c__1, &work[1], &c__2);
	}
	if (e[i__] < 0.f) {
	    scopy_(n, &v[km], &c__1, &work[1], &c__2);
	}
	if (e[i__] == 0.f) {
	    scopy_(n, &c_b39, &c__0, &work[2], &c__2);
	}
	if (e[i__] > 0.f) {
	    scopy_(n, &v[kp], &c__1, &work[2], &c__2);
	}
	if (e[i__] < 0.f) {
	    scopym_(n, &v[k], &c__1, &work[2], &c__2);
	}
	l = (j - 1 << 1) * *ldv + 1;
	i__2 = *n << 1;
	scopy_(&i__2, &work[1], &c__1, &v[l], &c__1);
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
    e[2] = 0.f;
    *info = 0;
    if (*job == 0) {
	return 0;
    }
    v[1] = a[1];
    v[2] = 0.f;
    return 0;
} /* sgeev_ */

