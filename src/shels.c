/* shels.f -- translated by f2c (version 12.02.01).
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

/* DECK SHELS */
/* Subroutine */ int shels_(real *a, integer *lda, integer *n, real *q, real *
	b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static real c__;
    static integer k;
    static real s, t, t1, t2;
    static integer kb, iq, kp1;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SHELS */
/* ***SUBSIDIARY */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      SINGLE PRECISION (SHELS-S, DHELS-D) */
/* ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION, */
/*             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE */
/* ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@llnl.gov */
/*           Seager, Mark K., (LLNL), seager@llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO Box 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/* ***DESCRIPTION */
/*        This routine is extracted from the LINPACK routine SGESL with */
/*        changes due to the fact that A is an upper Hessenberg matrix. */

/*        SHELS solves the least squares problem: */

/*                   MIN(B-A*X,B-A*X) */

/*        using the factors computed by SHEQR. */

/* *Usage: */
/*      INTEGER LDA, N */
/*      REAL A(LDA,N), Q(2*N), B(N+1) */

/*      CALL SHELS(A, LDA, N, Q, B) */

/* *Arguments: */
/* A       :IN       Real A(LDA,N) */
/*          The output from SHEQR which contains the upper */
/*          triangular factor R in the QR decomposition of A. */
/* LDA     :IN       Integer */
/*          The leading dimension of the array A. */
/* N       :IN       Integer */
/*          A is originally an (N+1) by N matrix. */
/* Q       :IN       Real Q(2*N) */
/*          The coefficients of the N Givens rotations */
/*          used in the QR factorization of A. */
/* B       :INOUT    Real B(N+1) */
/*          On input, B is the right hand side vector. */
/*          On output, B is the solution vector X. */

/* ***SEE ALSO  SGMRES */
/* ***ROUTINES CALLED  SAXPY */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871001  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF) */
/*   910506  Made subsidiary to SGMRES.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/* ***END PROLOGUE  SHELS */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  SHELS */

/*         Minimize(B-A*X,B-A*X).  First form Q*B. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	iq = (k - 1 << 1) + 1;
	c__ = q[iq];
	s = q[iq + 1];
	t1 = b[k];
	t2 = b[kp1];
	b[k] = c__ * t1 - s * t2;
	b[kp1] = s * t1 + c__ * t2;
/* L20: */
    }

/*         Now solve  R*X = Q*B. */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L40: */
    }
    return 0;
/* ------------- LAST LINE OF SHELS FOLLOWS ---------------------------- */
} /* shels_ */

