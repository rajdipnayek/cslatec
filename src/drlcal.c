/* drlcal.f -- translated by f2c (version 12.02.01).
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

/* DECK DRLCAL */
/* Subroutine */ int drlcal_(integer *n, integer *kmp, integer *ll, integer *
	maxl, doublereal *v, doublereal *q, doublereal *rl, doublereal *
	snormw, doublereal *prod, doublereal *r0nrm)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer i__, k;
    static doublereal s;
    static integer i2, ip1;
    static doublereal tem;
    static integer llm1, llp1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);

/* ***BEGIN PROLOGUE  DRLCAL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Internal routine for DGMRES. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      DOUBLE PRECISION (SRLCAL-S, DRLCAL-D) */
/* ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION, */
/*             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE */
/* ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@llnl.gov */
/*           Seager, Mark K., (LLNL), seager@llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO Box 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/* ***DESCRIPTION */
/*         This routine calculates the scaled residual RL from the */
/*         V(I)'s. */
/* *Usage: */
/*      INTEGER N, KMP, LL, MAXL */
/*      DOUBLE PRECISION V(N,LL), Q(2*MAXL), RL(N), SNORMW, PROD, R0NORM */

/*      CALL DRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, R0NRM) */

/* *Arguments: */
/* N      :IN       Integer */
/*         The order of the matrix A, and the lengths */
/*         of the vectors SR, SZ, R0 and Z. */
/* KMP    :IN       Integer */
/*         The number of previous V vectors the new vector VNEW */
/*         must be made orthogonal to. (KMP .le. MAXL) */
/* LL     :IN       Integer */
/*         The current dimension of the Krylov subspace. */
/* MAXL   :IN       Integer */
/*         The maximum dimension of the Krylov subspace. */
/* V      :IN       Double Precision V(N,LL) */
/*         The N x LL array containing the orthogonal vectors */
/*         V(*,1) to V(*,LL). */
/* Q      :IN       Double Precision Q(2*MAXL) */
/*         A double precision array of length 2*MAXL containing the */
/*         components of the Givens rotations used in the QR */
/*         decomposition of HES.  It is loaded in DHEQR and used in */
/*         DHELS. */
/* RL     :OUT      Double Precision RL(N) */
/*         The residual vector RL.  This is either SB*(B-A*XL) if */
/*         not preconditioning or preconditioning on the right, */
/*         or SB*(M-inverse)*(B-A*XL) if preconditioning on the */
/*         left. */
/* SNORMW :IN       Double Precision */
/*         Scale factor. */
/* PROD   :IN       Double Precision */
/*         The product s1*s2*...*sl = the product of the sines of the */
/*         Givens rotations used in the QR factorization of */
/*         the Hessenberg matrix HES. */
/* R0NRM  :IN       Double Precision */
/*         The scaled norm of initial residual R0. */

/* ***SEE ALSO  DGMRES */
/* ***ROUTINES CALLED  DCOPY, DSCAL */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890404  DATE WRITTEN */
/*   890404  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910506  Made subsidiary to DGMRES.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/* ***END PROLOGUE  DRLCAL */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  DRLCAL */
    /* Parameter adjustments */
    --rl;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --q;

    /* Function Body */
    if (*kmp == *maxl) {

/*         calculate RL.  Start by copying V(*,1) into RL. */

	dcopy_(n, &v[v_dim1 + 1], &c__1, &rl[1], &c__1);
	llm1 = *ll - 1;
	i__1 = llm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ip1 = i__ + 1;
	    i2 = i__ << 1;
	    s = q[i2];
	    c__ = q[i2 - 1];
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		rl[k] = s * rl[k] + c__ * v[k + ip1 * v_dim1];
/* L10: */
	    }
/* L20: */
	}
	s = q[*ll * 2];
	c__ = q[(*ll << 1) - 1] / *snormw;
	llp1 = *ll + 1;
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    rl[k] = s * rl[k] + c__ * v[k + llp1 * v_dim1];
/* L30: */
	}
    }

/*         When KMP < MAXL, RL vector already partially calculated. */
/*         Scale RL by R0NRM*PROD to obtain the residual RL. */

    tem = *r0nrm * *prod;
    dscal_(n, &tem, &rl[1], &c__1);
    return 0;
/* ------------- LAST LINE OF DRLCAL FOLLOWS ---------------------------- */
} /* drlcal_ */

