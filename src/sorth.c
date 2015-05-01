/* sorth.f -- translated by f2c (version 12.02.01).
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

/* DECK SORTH */
/* Subroutine */ int sorth_(real *vnew, real *v, real *hes, integer *n, 
	integer *ll, integer *ldhes, integer *kmp, real *snormw)
{
    /* System generated locals */
    integer hes_dim1, hes_offset, v_dim1, v_offset, i__1, i__2;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__, i0;
    static real arg, tem;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real vnrm;
    extern doublereal snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);
    static real sumdsq;

/* ***BEGIN PROLOGUE  SORTH */
/* ***SUBSIDIARY */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      SINGLE PRECISION (SORTH-S, DORTH-D) */
/* ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION, */
/*             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE */
/* ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@llnl.gov */
/*           Seager, Mark K., (LLNL), seager@llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO Box 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/* ***DESCRIPTION */
/*        This routine  orthogonalizes  the  vector  VNEW  against the */
/*        previous KMP  vectors in the   V array.  It uses  a modified */
/*        Gram-Schmidt   orthogonalization procedure with  conditional */
/*        reorthogonalization. */

/* *Usage: */
/*      INTEGER N, LL, LDHES, KMP */
/*      REAL VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW */

/*      CALL SORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW) */

/* *Arguments: */
/* VNEW   :INOUT    Real VNEW(N) */
/*         On input, the vector of length N containing a scaled */
/*         product of the Jacobian and the vector V(*,LL). */
/*         On output, the new vector orthogonal to V(*,i0) to V(*,LL), */
/*         where i0 = max(1, LL-KMP+1). */
/* V      :IN       Real V(N,LL) */
/*         The N x LL array containing the previous LL */
/*         orthogonal vectors V(*,1) to V(*,LL). */
/* HES    :INOUT    Real HES(LDHES,LL) */
/*         On input, an LL x LL upper Hessenberg matrix containing, */
/*         in HES(I,K), K.lt.LL, the scaled inner products of */
/*         A*V(*,K) and V(*,i). */
/*         On return, column LL of HES is filled in with */
/*         the scaled inner products of A*V(*,LL) and V(*,i). */
/* N      :IN       Integer */
/*         The order of the matrix A, and the length of VNEW. */
/* LL     :IN       Integer */
/*         The current order of the matrix HES. */
/* LDHES  :IN       Integer */
/*         The leading dimension of the HES array. */
/* KMP    :IN       Integer */
/*         The number of previous vectors the new vector VNEW */
/*         must be made orthogonal to (KMP .le. MAXL). */
/* SNORMW :OUT      REAL */
/*         Scalar containing the l-2 norm of VNEW. */

/* ***SEE ALSO  SGMRES */
/* ***ROUTINES CALLED  SAXPY, SDOT, SNRM2 */
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
/* ***END PROLOGUE  SORTH */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  SORTH */

/*         Get norm of unaltered VNEW for later use. */

    /* Parameter adjustments */
    --vnew;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    hes_dim1 = *ldhes;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;

    /* Function Body */
    vnrm = snrm2_(n, &vnew[1], &c__1);
/*   ------------------------------------------------------------------- */
/*         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL). */
/*         Scaled inner products give new column of HES. */
/*         Projections of earlier vectors are subtracted from VNEW. */
/*   ------------------------------------------------------------------- */
/* Computing MAX */
    i__1 = 1, i__2 = *ll - *kmp + 1;
    i0 = max(i__1,i__2);
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	hes[i__ + *ll * hes_dim1] = sdot_(n, &v[i__ * v_dim1 + 1], &c__1, &
		vnew[1], &c__1);
	tem = -hes[i__ + *ll * hes_dim1];
	saxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* L10: */
    }
/*   ------------------------------------------------------------------- */
/*         Compute SNORMW = norm of VNEW.  If VNEW is small compared */
/*         to its input value (in norm), then reorthogonalize VNEW to */
/*         V(*,1) through V(*,LL).  Correct if relative correction */
/*         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using */
/*         the dot products involved. */
/*   ------------------------------------------------------------------- */
    *snormw = snrm2_(n, &vnew[1], &c__1);
    if (vnrm + *snormw * .001f != vnrm) {
	return 0;
    }
    sumdsq = 0.f;
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	tem = -sdot_(n, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
	if (hes[i__ + *ll * hes_dim1] + tem * .001f == hes[i__ + *ll * 
		hes_dim1]) {
	    goto L30;
	}
	hes[i__ + *ll * hes_dim1] -= tem;
	saxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* Computing 2nd power */
	r__1 = tem;
	sumdsq += r__1 * r__1;
L30:
	;
    }
    if (sumdsq == 0.f) {
	return 0;
    }
/* Computing MAX */
/* Computing 2nd power */
    r__3 = *snormw;
    r__1 = 0.f, r__2 = r__3 * r__3 - sumdsq;
    arg = dmax(r__1,r__2);
    *snormw = sqrt(arg);

    return 0;
/* ------------- LAST LINE OF SORTH FOLLOWS ---------------------------- */
} /* sorth_ */

