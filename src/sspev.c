/* sspev.f -- translated by f2c (version 12.02.01).
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

/* DECK SSPEV */
/* Subroutine */ int sspev_(real *a, integer *n, real *e, real *v, integer *
	ldv, real *work, integer *job, integer *info)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, m;
    extern /* Subroutine */ int tred3_(integer *, integer *, real *, real *, 
	    real *, real *), trbak3_(integer *, integer *, integer *, real *, 
	    integer *, real *), imtql2_(integer *, integer *, real *, real *, 
	    real *, integer *), xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), tqlrat_(integer *, real *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SSPEV */
/* ***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors */
/*            of a real symmetric matrix stored in packed form. */
/* ***LIBRARY   SLATEC (EISPACK) */
/* ***CATEGORY  D4A1 */
/* ***TYPE      SINGLE PRECISION (SSPEV-S) */
/* ***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK, PACKED, SYMMETRIC */
/* ***AUTHOR  Kahaner, D. K., (NBS) */
/*           Moler, C. B., (U. of New Mexico) */
/*           Stewart, G. W., (U. of Maryland) */
/* ***DESCRIPTION */

/*     Abstract */
/*      SSPEV computes the eigenvalues and, optionally, the eigenvectors */
/*      of a real symmetric matrix stored in packed form. */

/*     Call Sequence Parameters- */
/*       (The values of parameters marked with * (star) will be  changed */
/*         by SSPEV.) */

/*        A*      REAL(N*(N+1)/2) */
/*                real symmetric packed input matrix.  Contains upper */
/*                triangle and diagonal of A, by column (elements */
/*                11, 12, 22, 13, 23, 33, ...). */

/*        N       INTEGER */
/*                set by the user to */
/*                the order of the matrix A. */

/*        E*      REAL(N) */
/*                on return from SSPEV, E contains the eigenvalues of A. */
/*                See also INFO below. */

/*        V*      REAL(LDV,N) */
/*                on return from SSPEV, if the user has set JOB */
/*                = 0        V is not referenced. */
/*                = nonzero  the N eigenvectors of A are stored in the */
/*                first N columns of V.  See also INFO below. */

/*        LDV     INTEGER */
/*                set by the user to */
/*                the leading dimension of the array V if JOB is also */
/*                set nonzero.  In that case, N must be .LE. LDV. */
/*                If JOB is set to zero, LDV is not referenced. */

/*        WORK*   REAL(2N) */
/*                temporary storage vector.  Contents changed by SSPEV. */

/*        JOB     INTEGER */
/*                set by the user to */
/*                = 0        eigenvalues only to be calculated by SSPEV. */
/*                           Neither V nor LDV are referenced. */
/*                = nonzero  eigenvalues and vectors to be calculated. */
/*                           In this case, A & V must be distinct arrays. */
/*                           Also, if LDA .GT. LDV, SSPEV changes all the */
/*                           elements of A thru column N.  If LDA < LDV, */
/*                           SSPEV changes all the elements of V through */
/*                           column N.  If LDA=LDV, only A(I,J) and V(I, */
/*                           J) for I,J = 1,...,N are changed by SSPEV. */

/*       INFO*   INTEGER */
/*               on return from SSPEV, the value of INFO is */
/*               = 0 for normal return. */
/*               = K if the eigenvalue iteration fails to converge. */
/*                   Eigenvalues and vectors 1 through K-1 are correct. */


/*     Error Messages- */
/*          No. 1   recoverable  N is greater than LDV and JOB is nonzero */
/*          No. 2   recoverable  N is less than one */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  IMTQL2, TQLRAT, TRBAK3, TRED3, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800808  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   900326  Removed duplicate information from DESCRIPTION section. */
/*           (WRB) */
/* ***END PROLOGUE  SSPEV */
/* ***FIRST EXECUTABLE STATEMENT  SSPEV */
    /* Parameter adjustments */
    --a;
    --e;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;

    /* Function Body */
    if (*n > *ldv) {
	xermsg_("SLATEC", "SSPEV", "N .GT. LDV.", &c__1, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)11);
    }
    if (*n > *ldv) {
	return 0;
    }
    if (*n < 1) {
	xermsg_("SLATEC", "SSPEV", "N .LT. 1", &c__2, &c__1, (ftnlen)6, (
		ftnlen)5, (ftnlen)8);
    }
    if (*n < 1) {
	return 0;
    }

/*       CHECK N=1 CASE */

    e[1] = a[1];
    *info = 0;
    if (*n == 1) {
	return 0;
    }

    if (*job != 0) {
	goto L20;
    }

/*     EIGENVALUES ONLY */

    tred3_(n, &c__1, &a[1], &e[1], &work[1], &work[*n + 1]);
    tqlrat_(n, &e[1], &work[*n + 1], info);
    return 0;

/*     EIGENVALUES AND EIGENVECTORS */

L20:
    tred3_(n, &c__1, &a[1], &e[1], &work[1], &work[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
/* L25: */
	    v[i__ + j * v_dim1] = 0.f;
	}
/* L30: */
	v[i__ + i__ * v_dim1] = 1.f;
    }
    imtql2_(ldv, n, &e[1], &work[1], &v[v_offset], info);
    m = *n;
    if (*info != 0) {
	m = *info - 1;
    }
    trbak3_(ldv, n, &c__1, &a[1], &m, &v[v_offset]);
    return 0;
} /* sspev_ */

