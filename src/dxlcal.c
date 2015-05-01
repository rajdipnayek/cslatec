/* dxlcal.f -- translated by f2c (version 12.02.01).
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

/* DECK DXLCAL */
/* Subroutine */ int dxlcal_(integer *n, integer *lgmr, doublereal *x, 
	doublereal *xl, doublereal *zl, doublereal *hes, integer *maxlp1, 
	doublereal *q, doublereal *v, doublereal *r0nrm, doublereal *wk, 
	doublereal *sz, integer *jscal, integer *jpre, S_fp msolve, integer *
	nmsl, doublereal *rpar, integer *ipar, integer *nelt, integer *ia, 
	integer *ja, doublereal *a, integer *isym)
{
    /* System generated locals */
    integer hes_dim1, hes_offset, v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k, ll, llp1;
    extern /* Subroutine */ int dhels_(doublereal *, integer *, integer *, 
	    doublereal *, doublereal *), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), daxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *);

/* ***BEGIN PROLOGUE  DXLCAL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Internal routine for DGMRES. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      DOUBLE PRECISION (SXLCAL-S, DXLCAL-D) */
/* ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION, */
/*             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE */
/* ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@llnl.gov */
/*           Seager, Mark K., (LLNL), seager@llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO Box 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/* ***DESCRIPTION */
/*        This  routine computes the solution  XL,  the current DGMRES */
/*        iterate, given the  V(I)'s and  the  QR factorization of the */
/*        Hessenberg  matrix HES.   This routine  is  only called when */
/*        ITOL=11. */

/* *Usage: */
/*      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, NMSL, IPAR(USER DEFINED) */
/*      INTEGER NELT, IA(NELT), JA(NELT), ISYM */
/*      DOUBLE PRECISION X(N), XL(N), ZL(N), HES(MAXLP1,MAXL), Q(2*MAXL), */
/*     $                 V(N,MAXLP1), R0NRM, WK(N), SZ(N), */
/*     $                 RPAR(USER DEFINED), A(NELT) */
/*      EXTERNAL MSOLVE */

/*      CALL DXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM, */
/*     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, */
/*     $     NELT, IA, JA, A, ISYM) */

/* *Arguments: */
/* N      :IN       Integer */
/*         The order of the matrix A, and the lengths */
/*         of the vectors SR, SZ, R0 and Z. */
/* LGMR   :IN       Integer */
/*         The number of iterations performed and */
/*         the current order of the upper Hessenberg */
/*         matrix HES. */
/* X      :IN       Double Precision X(N) */
/*         The current approximate solution as of the last restart. */
/* XL     :OUT      Double Precision XL(N) */
/*         An array of length N used to hold the approximate */
/*         solution X(L). */
/*         Warning: XL and ZL are the same array in the calling routine. */
/* ZL     :IN       Double Precision ZL(N) */
/*         An array of length N used to hold the approximate */
/*         solution Z(L). */
/* HES    :IN       Double Precision HES(MAXLP1,MAXL) */
/*         The upper triangular factor of the QR decomposition */
/*         of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*         entries are the scaled inner-products of A*V(*,i) and V(*,k). */
/* MAXLP1 :IN       Integer */
/*         MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES. */
/*         MAXL is the maximum allowable order of the matrix HES. */
/* Q      :IN       Double Precision Q(2*MAXL) */
/*         A double precision array of length 2*MAXL containing the */
/*         components of the Givens rotations used in the QR */
/*         decomposition of HES.  It is loaded in DHEQR. */
/* V      :IN       Double Precision V(N,MAXLP1) */
/*         The N by(LGMR+1) array containing the LGMR */
/*         orthogonal vectors V(*,1) to V(*,LGMR). */
/* R0NRM  :IN       Double Precision */
/*         The scaled norm of the initial residual for the */
/*         current call to DPIGMR. */
/* WK     :IN       Double Precision WK(N) */
/*         A double precision work array of length N. */
/* SZ     :IN       Double Precision SZ(N) */
/*         A vector of length N containing the non-zero */
/*         elements of the diagonal scaling matrix for Z. */
/* JSCAL  :IN       Integer */
/*         A flag indicating whether arrays SR and SZ are used. */
/*         JSCAL=0 means SR and SZ are not used and the */
/*                 algorithm will perform as if all */
/*                 SR(i) = 1 and SZ(i) = 1. */
/*         JSCAL=1 means only SZ is used, and the algorithm */
/*                 performs as if all SR(i) = 1. */
/*         JSCAL=2 means only SR is used, and the algorithm */
/*                 performs as if all SZ(i) = 1. */
/*         JSCAL=3 means both SR and SZ are used. */
/* JPRE   :IN       Integer */
/*         The preconditioner type flag. */
/* MSOLVE :EXT      External. */
/*         Name of the routine which solves a linear system Mz = r for */
/*         z given r with the preconditioning matrix M (M is supplied via */
/*         RPAR and IPAR arrays.  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSOLVE is: */
/*             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR) */
/*         Where N is the number of unknowns, R is the right-hand side */
/*         vector and Z is the solution upon return.  NELT, IA, JA, A and */
/*         ISYM are defined as below.  RPAR is a double precision array */
/*         that can be used to pass necessary preconditioning information */
/*         and/or workspace to MSOLVE.  IPAR is an integer work array */
/*         for the same purpose as RPAR. */
/* NMSL   :IN       Integer */
/*         The number of calls to MSOLVE. */
/* RPAR   :IN       Double Precision RPAR(USER DEFINED) */
/*         Double Precision workspace passed directly to the MSOLVE */
/*         routine. */
/* IPAR   :IN       Integer IPAR(USER DEFINED) */
/*         Integer workspace passed directly to the MSOLVE routine. */
/* NELT   :IN       Integer */
/*         The length of arrays IA, JA and A. */
/* IA     :IN       Integer IA(NELT) */
/*         An integer array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* JA     :IN       Integer JA(NELT) */
/*         An integer array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* A      :IN       Double Precision A(NELT) */
/*         A double precision array of length NELT containing matrix */
/*         data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* ISYM   :IN       Integer */
/*         A flag to indicate symmetric matrix storage. */
/*         If ISYM=0, all non-zero entries of the matrix are */
/*         stored.  If ISYM=1, the matrix is symmetric and */
/*         only the upper or lower triangular part is stored. */

/* ***SEE ALSO  DGMRES */
/* ***ROUTINES CALLED  DAXPY, DCOPY, DHELS */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890404  DATE WRITTEN */
/*   890404  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF) */
/*   910506  Made subsidiary to DGMRES.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/* ***END PROLOGUE  DXLCAL */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  DXLCAL */
    /* Parameter adjustments */
    --wk;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --zl;
    --xl;
    --x;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --q;
    --sz;
    --rpar;
    --ipar;
    --a;
    --ja;
    --ia;

    /* Function Body */
    ll = *lgmr;
    llp1 = ll + 1;
    i__1 = llp1;
    for (k = 1; k <= i__1; ++k) {
	wk[k] = 0.;
/* L10: */
    }
    wk[1] = *r0nrm;
    dhels_(&hes[hes_offset], maxlp1, &ll, &q[1], &wk[1]);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	zl[k] = 0.;
/* L20: */
    }
    i__1 = ll;
    for (i__ = 1; i__ <= i__1; ++i__) {
	daxpy_(n, &wk[i__], &v[i__ * v_dim1 + 1], &c__1, &zl[1], &c__1);
/* L30: */
    }
    if (*jscal == 1 || *jscal == 3) {
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    zl[k] /= sz[k];
/* L40: */
	}
    }
    if (*jpre > 0) {
	dcopy_(n, &zl[1], &c__1, &wk[1], &c__1);
	(*msolve)(n, &wk[1], &zl[1], nelt, &ia[1], &ja[1], &a[1], isym, &rpar[
		1], &ipar[1]);
	++(*nmsl);
    }
/*         calculate XL from X and ZL. */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	xl[k] = x[k] + zl[k];
/* L50: */
    }
    return 0;
/* ------------- LAST LINE OF DXLCAL FOLLOWS ---------------------------- */
} /* dxlcal_ */

