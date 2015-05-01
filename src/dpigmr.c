/* dpigmr.f -- translated by f2c (version 12.02.01).
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

/* DECK DPIGMR */
/* Subroutine */ int dpigmr_(integer *n, doublereal *r0, doublereal *sr, 
	doublereal *sz, integer *jscal, integer *maxl, integer *maxlp1, 
	integer *kmp, integer *nrsts, integer *jpre, S_fp matvec, S_fp msolve,
	 integer *nmsl, doublereal *z__, doublereal *v, doublereal *hes, 
	doublereal *q, integer *lgmr, doublereal *rpar, integer *ipar, 
	doublereal *wk, doublereal *dl, doublereal *rhol, integer *nrmax, 
	doublereal *b, doublereal *bnrm, doublereal *x, doublereal *xl, 
	integer *itol, doublereal *tol, integer *nelt, integer *ia, integer *
	ja, doublereal *a, integer *isym, integer *iunit, integer *iflag, 
	doublereal *err)
{
    /* System generated locals */
    integer hes_dim1, hes_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, j, k;
    static doublereal s;
    static integer i2, ll, ip1;
    static doublereal tem, rho;
    static integer llp1, info;
    static doublereal prod;
    static integer iter;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal r0nrm;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dhels_(doublereal *, integer *, integer *, doublereal 
	    *, doublereal *), dheqr_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *);
    static doublereal dlnrm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dorth_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *);
    static integer itmax;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), drlcal_(integer *, integer *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern integer isdgmr_(integer *, doublereal *, doublereal *, doublereal *
	    , integer *, integer *, integer *, doublereal *, integer *, S_fp, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    , doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *);
    static doublereal snormw;

/* ***BEGIN PROLOGUE  DPIGMR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Internal routine for DGMRES. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      DOUBLE PRECISION (SPIGMR-S, DPIGMR-D) */
/* ***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION, */
/*             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE */
/* ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@llnl.gov */
/*           Seager, Mark K., (LLNL), seager@llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO Box 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/* ***DESCRIPTION */
/*         This routine solves the linear system A * Z = R0 using a */
/*         scaled preconditioned version of the generalized minimum */
/*         residual method.  An initial guess of Z = 0 is assumed. */

/* *Usage: */
/*      INTEGER N, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, NMSL, LGMR */
/*      INTEGER IPAR(USER DEFINED), NRMAX, ITOL, NELT, IA(NELT), JA(NELT) */
/*      INTEGER ISYM, IUNIT, IFLAG */
/*      DOUBLE PRECISION R0(N), SR(N), SZ(N), Z(N), V(N,MAXLP1), */
/*     $                 HES(MAXLP1,MAXL), Q(2*MAXL), RPAR(USER DEFINED), */
/*     $                 WK(N), DL(N), RHOL, B(N), BNRM, X(N), XL(N), */
/*     $                 TOL, A(NELT), ERR */
/*      EXTERNAL MATVEC, MSOLVE */

/*      CALL DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, */
/*     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, */
/*     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL, */
/*     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR) */

/* *Arguments: */
/* N      :IN       Integer */
/*         The order of the matrix A, and the lengths */
/*         of the vectors SR, SZ, R0 and Z. */
/* R0     :IN       Double Precision R0(N) */
/*         R0 = the right hand side of the system A*Z = R0. */
/*         R0 is also used as workspace when computing */
/*         the final approximation. */
/*         (R0 is the same as V(*,MAXL+1) in the call to DPIGMR.) */
/* SR     :IN       Double Precision SR(N) */
/*         SR is a vector of length N containing the non-zero */
/*         elements of the diagonal scaling matrix for R0. */
/* SZ     :IN       Double Precision SZ(N) */
/*         SZ is a vector of length N containing the non-zero */
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
/* MAXL   :IN       Integer */
/*         The maximum allowable order of the matrix H. */
/* MAXLP1 :IN       Integer */
/*         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES. */
/* KMP    :IN       Integer */
/*         The number of previous vectors the new vector VNEW */
/*         must be made orthogonal to.  (KMP .le. MAXL) */
/* NRSTS  :IN       Integer */
/*         Counter for the number of restarts on the current */
/*         call to DGMRES.  If NRSTS .gt. 0, then the residual */
/*         R0 is already scaled, and so scaling of it is */
/*         not necessary. */
/* JPRE   :IN       Integer */
/*         Preconditioner type flag. */
/* MATVEC :EXT      External. */
/*         Name of a routine which performs the matrix vector multiply */
/*         Y = A*X given A and X.  The name of the MATVEC routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MATVEC is: */
/*             CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM) */
/*         where N is the number of unknowns, Y is the product A*X */
/*         upon return, X is an input vector, and NELT is the number of */
/*         non-zeros in the SLAP IA, JA, A storage for the matrix A. */
/*         ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
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
/* NMSL   :OUT      Integer */
/*         The number of calls to MSOLVE. */
/* Z      :OUT      Double Precision Z(N) */
/*         The final computed approximation to the solution */
/*         of the system A*Z = R0. */
/* V      :OUT      Double Precision V(N,MAXLP1) */
/*         The N by (LGMR+1) array containing the LGMR */
/*         orthogonal vectors V(*,1) to V(*,LGMR). */
/* HES    :OUT      Double Precision HES(MAXLP1,MAXL) */
/*         The upper triangular factor of the QR decomposition */
/*         of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*         entries are the scaled inner-products of A*V(*,I) */
/*         and V(*,K). */
/* Q      :OUT      Double Precision Q(2*MAXL) */
/*         A double precision array of length 2*MAXL containing the */
/*         components of the Givens rotations used in the QR */
/*         decomposition of HES.  It is loaded in DHEQR and used in */
/*         DHELS. */
/* LGMR   :OUT      Integer */
/*         The number of iterations performed and */
/*         the current order of the upper Hessenberg */
/*         matrix HES. */
/* RPAR   :IN       Double Precision RPAR(USER DEFINED) */
/*         Double Precision workspace passed directly to the MSOLVE */
/*         routine. */
/* IPAR   :IN       Integer IPAR(USER DEFINED) */
/*         Integer workspace passed directly to the MSOLVE routine. */
/* WK     :IN       Double Precision WK(N) */
/*         A double precision work array of length N used by routines */
/*         MATVEC and MSOLVE. */
/* DL     :INOUT    Double Precision DL(N) */
/*         On input, a double precision work array of length N used for */
/*         calculation of the residual norm RHO when the method is */
/*         incomplete (KMP.lt.MAXL), and/or when using restarting. */
/*         On output, the scaled residual vector RL.  It is only loaded */
/*         when performing restarts of the Krylov iteration. */
/* RHOL   :OUT      Double Precision */
/*         A double precision scalar containing the norm of the final */
/*         residual. */
/* NRMAX  :IN       Integer */
/*         The maximum number of restarts of the Krylov iteration. */
/*         NRMAX .gt. 0 means restarting is active, while */
/*         NRMAX = 0 means restarting is not being used. */
/* B      :IN       Double Precision B(N) */
/*         The right hand side of the linear system A*X = b. */
/* BNRM   :IN       Double Precision */
/*         The scaled norm of b. */
/* X      :IN       Double Precision X(N) */
/*         The current approximate solution as of the last */
/*         restart. */
/* XL     :IN       Double Precision XL(N) */
/*         An array of length N used to hold the approximate */
/*         solution X(L) when ITOL=11. */
/* ITOL   :IN       Integer */
/*         A flag to indicate the type of convergence criterion */
/*         used.  See the driver for its description. */
/* TOL    :IN       Double Precision */
/*         The tolerance on residuals R0-A*Z in scaled norm. */
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
/*         data. It is passed directly to the MATVEC and MSOLVE routines. */
/* ISYM   :IN       Integer */
/*         A flag to indicate symmetric matrix storage. */
/*         If ISYM=0, all non-zero entries of the matrix are */
/*         stored.  If ISYM=1, the matrix is symmetric and */
/*         only the upper or lower triangular part is stored. */
/* IUNIT  :IN       Integer */
/*         The i/o unit number for writing intermediate residual */
/*         norm values. */
/* IFLAG  :OUT      Integer */
/*         An integer error flag.. */
/*         0 means convergence in LGMR iterations, LGMR.le.MAXL. */
/*         1 means the convergence test did not pass in MAXL */
/*           iterations, but the residual norm is .lt. norm(R0), */
/*           and so Z is computed. */
/*         2 means the convergence test did not pass in MAXL */
/*           iterations, residual .ge. norm(R0), and Z = 0. */
/* ERR    :OUT      Double Precision. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL. */

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */

/* ***SEE ALSO  DGMRES */
/* ***ROUTINES CALLED  DAXPY, DCOPY, DHELS, DHEQR, DNRM2, DORTH, DRLCAL, */
/*                    DSCAL, ISDGMR */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890404  DATE WRITTEN */
/*   890404  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF) */
/*   910506  Made subsidiary to DGMRES.  (FNF) */
/*   920511  Added complete declaration section.  (WRB) */
/* ***END PROLOGUE  DPIGMR */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  DPIGMR */

/*         Zero out the Z array. */

    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --r0;
    --sr;
    --sz;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --z__;
    --q;
    --rpar;
    --ipar;
    --wk;
    --dl;
    --b;
    --x;
    --xl;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = 0.;
/* L5: */
    }

    *iflag = 0;
    *lgmr = 0;
    *nmsl = 0;
/*         Load ITMAX, the maximum number of iterations. */
    itmax = (*nrmax + 1) * *maxl;
/*   ------------------------------------------------------------------- */
/*         The initial residual is the vector R0. */
/*         Apply left precon. if JPRE < 0 and this is not a restart. */
/*         Apply scaling to R0 if JSCAL = 2 or 3. */
/*   ------------------------------------------------------------------- */
    if (*jpre < 0 && *nrsts == 0) {
	dcopy_(n, &r0[1], &c__1, &wk[1], &c__1);
	(*msolve)(n, &wk[1], &r0[1], nelt, &ia[1], &ja[1], &a[1], isym, &rpar[
		1], &ipar[1]);
	++(*nmsl);
    }
    if ((*jscal == 2 || *jscal == 3) && *nrsts == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__ + v_dim1] = r0[i__] * sr[i__];
/* L10: */
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__ + v_dim1] = r0[i__];
/* L20: */
	}
    }
    r0nrm = dnrm2_(n, &v[v_offset], &c__1);
    iter = *nrsts * *maxl;

/*         Call stopping routine ISDGMR. */

    if (isdgmr_(n, &b[1], &x[1], &xl[1], nelt, &ia[1], &ja[1], &a[1], isym, (
	    S_fp)msolve, nmsl, itol, tol, &itmax, &iter, err, iunit, &v[
	    v_dim1 + 1], &z__[1], &wk[1], &rpar[1], &ipar[1], &r0nrm, bnrm, &
	    sr[1], &sz[1], jscal, kmp, lgmr, maxl, maxlp1, &v[v_offset], &q[1]
	    , &snormw, &prod, &r0nrm, &hes[hes_offset], jpre) != 0) {
	return 0;
    }
    tem = 1. / r0nrm;
    dscal_(n, &tem, &v[v_dim1 + 1], &c__1);

/*         Zero out the HES array. */

    i__1 = *maxl;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *maxlp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    hes[i__ + j * hes_dim1] = 0.;
/* L40: */
	}
/* L50: */
    }
/*   ------------------------------------------------------------------- */
/*         Main loop to compute the vectors V(*,2) to V(*,MAXL). */
/*         The running product PROD is needed for the convergence test. */
/*   ------------------------------------------------------------------- */
    prod = 1.;
    i__1 = *maxl;
    for (ll = 1; ll <= i__1; ++ll) {
	*lgmr = ll;
/*   ------------------------------------------------------------------- */
/*        Unscale  the  current V(LL)  and store  in WK.  Call routine */
/*        MSOLVE    to   compute(M-inverse)*WK,   where    M   is  the */
/*        preconditioner matrix.  Save the answer in Z.   Call routine */
/*        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system */
/*        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call */
/*        routine DORTH  to  orthogonalize the    new vector VNEW   = */
/*        V(*,LL+1).  Call routine DHEQR to update the factors of HES. */
/*   ------------------------------------------------------------------- */
	if (*jscal == 1 || *jscal == 3) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		wk[i__] = v[i__ + ll * v_dim1] / sz[i__];
/* L60: */
	    }
	} else {
	    dcopy_(n, &v[ll * v_dim1 + 1], &c__1, &wk[1], &c__1);
	}
	if (*jpre > 0) {
	    (*msolve)(n, &wk[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		    rpar[1], &ipar[1]);
	    ++(*nmsl);
	    (*matvec)(n, &z__[1], &v[(ll + 1) * v_dim1 + 1], nelt, &ia[1], &
		    ja[1], &a[1], isym);
	} else {
	    (*matvec)(n, &wk[1], &v[(ll + 1) * v_dim1 + 1], nelt, &ia[1], &ja[
		    1], &a[1], isym);
	}
	if (*jpre < 0) {
	    dcopy_(n, &v[(ll + 1) * v_dim1 + 1], &c__1, &wk[1], &c__1);
	    (*msolve)(n, &wk[1], &v[(ll + 1) * v_dim1 + 1], nelt, &ia[1], &ja[
		    1], &a[1], isym, &rpar[1], &ipar[1]);
	    ++(*nmsl);
	}
	if (*jscal == 2 || *jscal == 3) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		v[i__ + (ll + 1) * v_dim1] *= sr[i__];
/* L65: */
	    }
	}
	dorth_(&v[(ll + 1) * v_dim1 + 1], &v[v_offset], &hes[hes_offset], n, &
		ll, maxlp1, kmp, &snormw);
	hes[ll + 1 + ll * hes_dim1] = snormw;
	dheqr_(&hes[hes_offset], maxlp1, &ll, &q[1], &info, &ll);
	if (info == ll) {
	    goto L120;
	}
/*   ------------------------------------------------------------------- */
/*         Update RHO, the estimate of the norm of the residual R0-A*ZL. */
/*         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not */
/*         necessarily orthogonal for LL > KMP.  The vector DL must then */
/*         be computed, and its norm used in the calculation of RHO. */
/*   ------------------------------------------------------------------- */
	prod *= q[ll * 2];
	rho = (d__1 = prod * r0nrm, abs(d__1));
	if (ll > *kmp && *kmp < *maxl) {
	    if (ll == *kmp + 1) {
		dcopy_(n, &v[v_dim1 + 1], &c__1, &dl[1], &c__1);
		i__2 = *kmp;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ip1 = i__ + 1;
		    i2 = i__ << 1;
		    s = q[i2];
		    c__ = q[i2 - 1];
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
			dl[k] = s * dl[k] + c__ * v[k + ip1 * v_dim1];
/* L70: */
		    }
/* L75: */
		}
	    }
	    s = q[ll * 2];
	    c__ = q[(ll << 1) - 1] / snormw;
	    llp1 = ll + 1;
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		dl[k] = s * dl[k] + c__ * v[k + llp1 * v_dim1];
/* L80: */
	    }
	    dlnrm = dnrm2_(n, &dl[1], &c__1);
	    rho *= dlnrm;
	}
	*rhol = rho;
/*   ------------------------------------------------------------------- */
/*         Test for convergence.  If passed, compute approximation ZL. */
/*         If failed and LL < MAXL, then continue iterating. */
/*   ------------------------------------------------------------------- */
	iter = *nrsts * *maxl + *lgmr;
	if (isdgmr_(n, &b[1], &x[1], &xl[1], nelt, &ia[1], &ja[1], &a[1], 
		isym, (S_fp)msolve, nmsl, itol, tol, &itmax, &iter, err, 
		iunit, &dl[1], &z__[1], &wk[1], &rpar[1], &ipar[1], rhol, 
		bnrm, &sr[1], &sz[1], jscal, kmp, lgmr, maxl, maxlp1, &v[
		v_offset], &q[1], &snormw, &prod, &r0nrm, &hes[hes_offset], 
		jpre) != 0) {
	    goto L200;
	}
	if (ll == *maxl) {
	    goto L100;
	}
/*   ------------------------------------------------------------------- */
/*         Rescale so that the norm of V(1,LL+1) is one. */
/*   ------------------------------------------------------------------- */
	tem = 1. / snormw;
	dscal_(n, &tem, &v[(ll + 1) * v_dim1 + 1], &c__1);
/* L90: */
    }
L100:
    if (rho < r0nrm) {
	goto L150;
    }
L120:
    *iflag = 2;

/*         Load approximate solution with zero. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = 0.;
/* L130: */
    }
    return 0;
L150:
    *iflag = 1;

/*         Tolerance not met, but residual norm reduced. */

    if (*nrmax > 0) {

/*        If performing restarting (NRMAX > 0)  calculate the residual */
/*        vector RL and  store it in the DL  array.  If the incomplete */
/*        version is being used (KMP < MAXL) then DL has  already been */
/*        calculated up to a scaling factor.   Use DRLCAL to calculate */
/*        the scaled residual vector. */

	drlcal_(n, kmp, maxl, maxl, &v[v_offset], &q[1], &dl[1], &snormw, &
		prod, &r0nrm);
    }
/*   ------------------------------------------------------------------- */
/*         Compute the approximation ZL to the solution.  Since the */
/*         vector Z was used as workspace, and the initial guess */
/*         of the linear iteration is zero, Z must be reset to zero. */
/*   ------------------------------------------------------------------- */
L200:
    ll = *lgmr;
    llp1 = ll + 1;
    i__1 = llp1;
    for (k = 1; k <= i__1; ++k) {
	r0[k] = 0.;
/* L210: */
    }
    r0[1] = r0nrm;
    dhels_(&hes[hes_offset], maxlp1, &ll, &q[1], &r0[1]);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	z__[k] = 0.;
/* L220: */
    }
    i__1 = ll;
    for (i__ = 1; i__ <= i__1; ++i__) {
	daxpy_(n, &r0[i__], &v[i__ * v_dim1 + 1], &c__1, &z__[1], &c__1);
/* L230: */
    }
    if (*jscal == 1 || *jscal == 3) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__[i__] /= sz[i__];
/* L240: */
	}
    }
    if (*jpre > 0) {
	dcopy_(n, &z__[1], &c__1, &wk[1], &c__1);
	(*msolve)(n, &wk[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		rpar[1], &ipar[1]);
	++(*nmsl);
    }
    return 0;
/* ------------- LAST LINE OF DPIGMR FOLLOWS ---------------------------- */
} /* dpigmr_ */

