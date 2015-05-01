/* sheqr.f -- translated by f2c (version 12.02.01).
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

/* DECK SHEQR */
/* Subroutine */ int sheqr_(real *a, integer *lda, integer *n, real *q, 
	integer *info, integer *ijob)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static real c__;
    static integer i__, j, k;
    static real s, t, t1, t2;
    static integer iq, km1, kp1, nm1;

/* ***BEGIN PROLOGUE  SHEQR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      SINGLE PRECISION (SHEQR-S, DHEQR-D) */
/* ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION, */
/*             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE */
/* ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@llnl.gov */
/*           Seager, Mark K., (LLNL), seager@llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO Box 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/* ***DESCRIPTION */
/*        This   routine  performs  a QR   decomposition  of an  upper */
/*        Hessenberg matrix A using Givens  rotations.  There  are two */
/*        options  available: 1)  Performing  a fresh decomposition 2) */
/*        updating the QR factors by adding a row and  a column to the */
/*        matrix A. */

/* *Usage: */
/*      INTEGER LDA, N, INFO, IJOB */
/*      REAL A(LDA,N), Q(2*N) */

/*      CALL SHEQR(A, LDA, N, Q, INFO, IJOB) */

/* *Arguments: */
/* A      :INOUT    Real A(LDA,N) */
/*         On input, the matrix to be decomposed. */
/*         On output, the upper triangular matrix R. */
/*         The factorization can be written Q*A = R, where */
/*         Q is a product of Givens rotations and R is upper */
/*         triangular. */
/* LDA    :IN       Integer */
/*         The leading dimension of the array A. */
/* N      :IN       Integer */
/*         A is an (N+1) by N Hessenberg matrix. */
/* Q      :OUT      Real Q(2*N) */
/*         The factors c and s of each Givens rotation used */
/*         in decomposing A. */
/* INFO   :OUT      Integer */
/*         = 0  normal value. */
/*         = K  if  A(K,K) .eq. 0.0 .  This is not an error */
/*           condition for this subroutine, but it does */
/*           indicate that SHELS will divide by zero */
/*           if called. */
/* IJOB   :IN       Integer */
/*         = 1     means that a fresh decomposition of the */
/*                 matrix A is desired. */
/*         .ge. 2  means that the current decomposition of A */
/*                 will be updated by the addition of a row */
/*                 and a column. */

/* ***SEE ALSO  SGMRES */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871001  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910506  Made subsidiary to SGMRES.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/* ***END PROLOGUE  SHEQR */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  SHEQR */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;

    /* Function Body */
    if (*ijob > 1) {
	goto L70;
    }
/*   ------------------------------------------------------------------- */
/*         A new factorization is desired. */
/*   ------------------------------------------------------------------- */
/*         QR decomposition without pivoting. */

    *info = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	km1 = k - 1;
	kp1 = k + 1;

/*           Compute K-th column of R. */
/*           First, multiply the K-th column of A by the previous */
/*           K-1 Givens rotations. */

	if (km1 < 1) {
	    goto L20;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    i__ = (j - 1 << 1) + 1;
	    t1 = a[j + k * a_dim1];
	    t2 = a[j + 1 + k * a_dim1];
	    c__ = q[i__];
	    s = q[i__ + 1];
	    a[j + k * a_dim1] = c__ * t1 - s * t2;
	    a[j + 1 + k * a_dim1] = s * t1 + c__ * t2;
/* L10: */
	}

/*         Compute Givens components C and S. */

L20:
	iq = (km1 << 1) + 1;
	t1 = a[k + k * a_dim1];
	t2 = a[kp1 + k * a_dim1];
	if (t2 == 0.f) {
	    c__ = 1.f;
	    s = 0.f;
	} else if (dabs(t2) >= dabs(t1)) {
	    t = t1 / t2;
	    s = -1.f / sqrt(t * t + 1.f);
	    c__ = -s * t;
	} else {
	    t = t2 / t1;
	    c__ = 1.f / sqrt(t * t + 1.f);
	    s = -c__ * t;
	}
	q[iq] = c__;
	q[iq + 1] = s;
	a[k + k * a_dim1] = c__ * t1 - s * t2;
	if (a[k + k * a_dim1] == 0.f) {
	    *info = k;
	}
/* L60: */
    }
    return 0;
/*   ------------------------------------------------------------------- */
/*         The old factorization of a will be updated.  A row and a */
/*         column has been added to the matrix A.  N by N-1 is now */
/*         the old size of the matrix. */
/*   ------------------------------------------------------------------- */
L70:
    nm1 = *n - 1;
/*   ------------------------------------------------------------------- */
/*         Multiply the new column by the N previous Givens rotations. */
/*   ------------------------------------------------------------------- */
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	i__ = (k - 1 << 1) + 1;
	t1 = a[k + *n * a_dim1];
	t2 = a[k + 1 + *n * a_dim1];
	c__ = q[i__];
	s = q[i__ + 1];
	a[k + *n * a_dim1] = c__ * t1 - s * t2;
	a[k + 1 + *n * a_dim1] = s * t1 + c__ * t2;
/* L100: */
    }
/*   ------------------------------------------------------------------- */
/*         Complete update of decomposition by forming last Givens */
/*         rotation, and multiplying it times the column */
/*         vector(A(N,N),A(NP1,N)). */
/*   ------------------------------------------------------------------- */
    *info = 0;
    t1 = a[*n + *n * a_dim1];
    t2 = a[*n + 1 + *n * a_dim1];
    if (t2 == 0.f) {
	c__ = 1.f;
	s = 0.f;
    } else if (dabs(t2) >= dabs(t1)) {
	t = t1 / t2;
	s = -1.f / sqrt(t * t + 1.f);
	c__ = -s * t;
    } else {
	t = t2 / t1;
	c__ = 1.f / sqrt(t * t + 1.f);
	s = -c__ * t;
    }
    iq = (*n << 1) - 1;
    q[iq] = c__;
    q[iq + 1] = s;
    a[*n + *n * a_dim1] = c__ * t1 - s * t2;
    if (a[*n + *n * a_dim1] == 0.f) {
	*info = *n;
    }
    return 0;
/* ------------- LAST LINE OF SHEQR FOLLOWS ---------------------------- */
} /* sheqr_ */

