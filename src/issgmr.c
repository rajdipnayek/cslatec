/* issgmr.f -- translated by f2c (version 12.02.01).
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

/* Common Block Declarations */

struct {
    real soln[1];
} sslblk_;

#define sslblk_1 sslblk_

/* Table of constant values */

static integer c__1 = 1;

/* DECK ISSGMR */
integer issgmr_(integer *n, real *b, real *x, real *xl, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, S_fp msolve, 
	integer *nmsl, integer *itol, real *tol, integer *itmax, integer *
	iter, real *err, integer *iunit, real *r__, real *z__, real *dz, real 
	*rwork, integer *iwork, real *rnrm, real *bnrm, real *sb, real *sx, 
	integer *jscal, integer *kmp, integer *lgmr, integer *maxl, integer *
	maxlp1, real *v, real *q, real *snormw, real *prod, real *r0nrm, real 
	*hes, integer *jpre)
{
    /* Format strings */
    static char fmt_1020[] = "(1x,\002 ITER = \002,i5,\002 IELMAX = \002,i5"
	    ",\002 |R(IELMAX)/X(IELMAX)| = \002,e12.5)";
    static char fmt_1000[] = "(\002 Generalized Minimum Residual(\002,i3,i3"
	    ",\002) for \002,\002N, ITOL = \002,i5,i5,/\002 ITER\002,\002   N"
	    "atural Err Est\002,\002   Error Estimate\002)";
    static char fmt_1010[] = "(1x,i4,1x,e16.7,1x,e16.7)";

    /* System generated locals */
    integer hes_dim1, hes_offset, v_dim1, v_offset, ret_val, i__1;
    real r__1, r__2, r__3;

    /* Local variables */
    static integer i__;
    static real tem, rat, fuzz;
    extern doublereal snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static real dxnrm;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    extern doublereal r1mach_(integer *);
    static integer ielmax;
    extern /* Subroutine */ int srlcal_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *), 
	    sxlcal_(integer *, integer *, real *, real *, real *, real *, 
	    integer *, real *, real *, real *, real *, real *, integer *, 
	    integer *, S_fp, integer *, real *, integer *, integer *, integer 
	    *, integer *, real *, integer *);
    static real ratmax, solnrm;

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_1010, 0 };


/* ***BEGIN PROLOGUE  ISSGMR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Generalized Minimum Residual Stop Test. */
/*            This routine calculates the stop test for the Generalized */
/*            Minimum RESidual (GMRES) iteration scheme.  It returns a */
/*            non-zero if the error estimate (the type of which is */
/*            determined by ITOL) is less than the user specified */
/*            tolerance TOL. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      SINGLE PRECISION (ISSGMR-S, ISDGMR-D) */
/* ***KEYWORDS  GMRES, LINEAR SYSTEM, SLAP, SPARSE, STOP TEST */
/* ***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@llnl.gov */
/*           Seager, Mark K., (LLNL), seager@llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO Box 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/* ***DESCRIPTION */

/* *Usage: */
/*      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NMSL, ITOL */
/*      INTEGER ITMAX, ITER, IUNIT, IWORK(USER DEFINED), JSCAL */
/*      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE */
/*      REAL B(N), X(N), XL(MAXL), A(NELT), TOL, ERR, R(N), Z(N), */
/*     $     DZ(N), RWORK(USER DEFINED), RNRM, BNRM, SB(N), SX(N), */
/*     $     V(N,MAXLP1), Q(2*MAXL), SNORMW, PROD, R0NRM, */
/*     $     HES(MAXLP1,MAXL) */
/*      EXTERNAL MSOLVE */

/*      IF (ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE, */
/*     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ, */
/*     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL, */
/*     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM, */
/*     $     HES, JPRE) .NE. 0) THEN ITERATION DONE */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand-side vector. */
/* X      :IN       Real X(N). */
/*         Approximate solution vector as of the last restart. */
/* XL     :OUT      Real XL(N) */
/*         An array of length N used to hold the approximate */
/*         solution as of the current iteration.  Only computed by */
/*         this routine when ITOL=11. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays contain the matrix data structure for A. */
/*         It could take any form.  See "Description", in the SGMRES, */
/*         SSLUGM and SSDGMR routines for more details. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* MSOLVE :EXT      External. */
/*         Name of a routine which solves a linear system Mz = r for  z */
/*         given r with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays.  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSOLVE is: */
/*             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) */
/*         Where N is the number of unknowns, R is the right-hand side */
/*         vector and Z is the solution upon return.  NELT, IA, JA, A and */
/*         ISYM are defined as above.  RWORK is a real array that can */
/*         be used to pass necessary preconditioning information and/or */
/*         workspace to MSOLVE.  IWORK is an integer work array for */
/*         the same purpose as RWORK. */
/* NMSL   :INOUT    Integer. */
/*         A counter for the number of calls to MSOLVE. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate the type of convergence criterion used. */
/*         ITOL=0  Means the  iteration stops when the test described */
/*                 below on  the  residual RL  is satisfied.  This is */
/*                 the  "Natural Stopping Criteria" for this routine. */
/*                 Other values  of   ITOL  cause  extra,   otherwise */
/*                 unnecessary, computation per iteration and     are */
/*                 therefore much less efficient. */
/*         ITOL=1  Means   the  iteration stops   when the first test */
/*                 described below on  the residual RL  is satisfied, */
/*                 and there  is either right  or  no preconditioning */
/*                 being used. */
/*         ITOL=2  Implies     that   the  user    is   using    left */
/*                 preconditioning, and the second stopping criterion */
/*                 below is used. */
/*         ITOL=3  Means the  iteration stops   when  the  third test */
/*                 described below on Minv*Residual is satisfied, and */
/*                 there is either left  or no  preconditioning begin */
/*                 used. */
/*         ITOL=11 is    often  useful  for   checking  and comparing */
/*                 different routines.  For this case, the  user must */
/*                 supply  the  "exact" solution or  a  very accurate */
/*                 approximation (one with  an  error much less  than */
/*                 TOL) through a common block, */
/*                     COMMON /SSLBLK/ SOLN( ) */
/*                 If ITOL=11, iteration stops when the 2-norm of the */
/*                 difference between the iterative approximation and */
/*                 the user-supplied solution  divided by the  2-norm */
/*                 of the  user-supplied solution  is  less than TOL. */
/*                 Note that this requires  the  user to  set up  the */
/*                 "COMMON     /SSLBLK/ SOLN(LENGTH)"  in the calling */
/*                 routine.  The routine with this declaration should */
/*                 be loaded before the stop test so that the correct */
/*                 length is used by  the loader.  This procedure  is */
/*                 not standard Fortran and may not work correctly on */
/*                 your   system (although  it  has  worked  on every */
/*                 system the authors have tried).  If ITOL is not 11 */
/*                 then this common block is indeed standard Fortran. */
/* TOL    :IN       Real. */
/*         Convergence criterion, as described above. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :IN       Integer. */
/*         The iteration for which to check for convergence. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL.  Letting norm() denote the Euclidean */
/*         norm, ERR is defined as follows.. */

/*         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               for right or no preconditioning, and */
/*                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               for left preconditioning. */
/*         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               since right or no preconditioning */
/*                               being used. */
/*         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               since left preconditioning is being */
/*                               used. */
/*         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)| */
/*                               i=1,n */
/*         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN). */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* R      :INOUT    Real R(N). */
/*         Work array used in calling routine.  It contains */
/*         information necessary to compute the residual RL = B-A*XL. */
/* Z      :WORK     Real Z(N). */
/*         Workspace used to hold the pseudo-residual M z = r. */
/* DZ     :WORK     Real DZ(N). */
/*         Workspace used to hold temporary vector(s). */
/* RWORK  :WORK     Real RWORK(USER DEFINED). */
/*         Real array that can be used by MSOLVE. */
/* IWORK  :WORK     Integer IWORK(USER DEFINED). */
/*         Integer array that can be used by MSOLVE. */
/* RNRM   :IN       Real. */
/*         Norm of the current residual.  Type of norm depends on ITOL. */
/* BNRM   :IN       Real. */
/*         Norm of the right hand side.  Type of norm depends on ITOL. */
/* SB     :IN       Real SB(N). */
/*         Scaling vector for B. */
/* SX     :IN       Real SX(N). */
/*         Scaling vector for X. */
/* JSCAL  :IN       Integer. */
/*         Flag indicating if scaling arrays SB and SX are being */
/*         used in the calling routine SPIGMR. */
/*         JSCAL=0 means SB and SX are not used and the */
/*                 algorithm will perform as if all */
/*                 SB(i) = 1 and SX(i) = 1. */
/*         JSCAL=1 means only SX is used, and the algorithm */
/*                 performs as if all SB(i) = 1. */
/*         JSCAL=2 means only SB is used, and the algorithm */
/*                 performs as if all SX(i) = 1. */
/*         JSCAL=3 means both SB and SX are used. */
/* KMP    :IN       Integer */
/*         The number of previous vectors the new vector VNEW */
/*         must be made orthogonal to.  (KMP .le. MAXL) */
/* LGMR   :IN       Integer */
/*         The number of GMRES iterations performed on the current call */
/*         to SPIGMR (i.e., # iterations since the last restart) and */
/*         the current order of the upper Hessenberg */
/*         matrix HES. */
/* MAXL   :IN       Integer */
/*         The maximum allowable order of the matrix H. */
/* MAXLP1 :IN       Integer */
/*         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES. */
/* V      :IN       Real V(N,MAXLP1) */
/*         The N by (LGMR+1) array containing the LGMR */
/*         orthogonal vectors V(*,1) to V(*,LGMR). */
/* Q      :IN       Real Q(2*MAXL) */
/*         A real array of length 2*MAXL containing the components */
/*         of the Givens rotations used in the QR decomposition */
/*         of HES. */
/* SNORMW :IN       Real */
/*         A scalar containing the scaled norm of VNEW before it */
/*         is renormalized in SPIGMR. */
/* PROD   :IN       Real */
/*         The product s1*s2*...*sl = the product of the sines of the */
/*         Givens rotations used in the QR factorization of the */
/*         Hessenberg matrix HES. */
/* R0NRM  :IN       Real */
/*         The scaled norm of initial residual R0. */
/* HES    :IN       Real HES(MAXLP1,MAXL) */
/*         The upper triangular factor of the QR decomposition */
/*         of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*         entries are the scaled inner-products of A*V(*,I) */
/*         and V(*,K). */
/* JPRE   :IN       Integer */
/*         Preconditioner type flag. */
/*         (See description of IGWK(4) in SGMRES.) */

/* *Description */
/*       When using the GMRES solver,  the preferred value  for ITOL */
/*       is 0.  This is due to the fact that when ITOL=0 the norm of */
/*       the residual required in the stopping test is  obtained for */
/*       free, since this value is already  calculated  in the GMRES */
/*       algorithm.   The  variable  RNRM contains the   appropriate */
/*       norm, which is equal to norm(SB*(RL - A*XL))  when right or */
/*       no   preconditioning is  being  performed,   and equal   to */
/*       norm(SB*Minv*(RL - A*XL))  when using left preconditioning. */
/*       Here, norm() is the Euclidean norm.  Nonzero values of ITOL */
/*       require  additional work  to  calculate the  actual  scaled */
/*       residual  or its scaled/preconditioned  form,  and/or   the */
/*       approximate solution XL.  Hence, these values of  ITOL will */
/*       not be as efficient as ITOL=0. */

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */

/*     This routine does not verify that ITOL has a valid value. */
/*     The calling routine should make such a test before calling */
/*     ISSGMR, as is done in SGMRES. */

/* ***SEE ALSO  SGMRES */
/* ***ROUTINES CALLED  R1MACH, SCOPY, SNRM2, SRLCAL, SSCAL, SXLCAL */
/* ***COMMON BLOCKS    SSLBLK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871211  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Corrected conversion errors, etc.  (FNF) */
/*   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF) */
/*   910506  Made subsidiary to SGMRES.  (FNF) */
/*   920407  COMMON BLOCK renamed SSLBLK.  (WRB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   921113  Corrected C***CATEGORY line.  (FNF) */
/* ***END PROLOGUE  ISSGMR */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Arrays in Common .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Common blocks .. */
/*     .. Save statement .. */
/* ***FIRST EXECUTABLE STATEMENT  ISSGMR */
    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --b;
    --x;
    --xl;
    --ia;
    --ja;
    --a;
    --r__;
    --z__;
    --dz;
    --rwork;
    --iwork;
    --sb;
    --sx;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --q;

    /* Function Body */
    ret_val = 0;
    if (*itol == 0) {

/*       Use input from SPIGMR to determine if stop conditions are met. */

	*err = *rnrm / *bnrm;
    }
    if (*itol > 0 && *itol <= 3) {

/*       Use SRLCAL to calculate the scaled residual vector. */
/*       Store answer in R. */

	if (*lgmr != 0) {
	    srlcal_(n, kmp, lgmr, maxl, &v[v_offset], &q[1], &r__[1], snormw, 
		    prod, r0nrm);
	}
	if (*itol <= 2) {
/*         err = ||Residual||/||RightHandSide||(2-Norms). */
	    *err = snrm2_(n, &r__[1], &c__1) / *bnrm;

/*         Unscale R by R0NRM*PROD when KMP < MAXL. */

	    if (*kmp < *maxl && *lgmr != 0) {
		tem = 1.f / (*r0nrm * *prod);
		sscal_(n, &tem, &r__[1], &c__1);
	    }
	} else if (*itol == 3) {
/*         err = Max |(Minv*Residual)(i)/x(i)| */
/*         When JPRE .lt. 0, R already contains Minv*Residual. */
	    if (*jpre > 0) {
		(*msolve)(n, &r__[1], &dz[1], nelt, &ia[1], &ja[1], &a[1], 
			isym, &rwork[1], &iwork[1]);
		++(*nmsl);
	    }

/*         Unscale R by R0NRM*PROD when KMP < MAXL. */

	    if (*kmp < *maxl && *lgmr != 0) {
		tem = 1.f / (*r0nrm * *prod);
		sscal_(n, &tem, &r__[1], &c__1);
	    }

	    fuzz = r1mach_(&c__1);
	    ielmax = 1;
/* Computing MAX */
	    r__1 = dabs(x[1]);
	    ratmax = dabs(dz[1]) / dmax(r__1,fuzz);
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		r__3 = (r__2 = x[i__], dabs(r__2));
		rat = (r__1 = dz[i__], dabs(r__1)) / dmax(r__3,fuzz);
		if (rat > ratmax) {
		    ielmax = i__;
		    ratmax = rat;
		}
/* L25: */
	    }
	    *err = ratmax;
	    if (ratmax <= *tol) {
		ret_val = 1;
	    }
	    if (*iunit > 0) {
		io___7.ciunit = *iunit;
		s_wsfe(&io___7);
		do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ielmax, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ratmax, (ftnlen)sizeof(real));
		e_wsfe();
	    }
	    return ret_val;
	}
    }
    if (*itol == 11) {

/*       Use SXLCAL to calculate the approximate solution XL. */

	if (*lgmr != 0 && *iter > 0) {
	    sxlcal_(n, lgmr, &x[1], &xl[1], &xl[1], &hes[hes_offset], maxlp1, 
		    &q[1], &v[v_offset], r0nrm, &dz[1], &sx[1], jscal, jpre, (
		    S_fp)msolve, nmsl, &rwork[1], &iwork[1], nelt, &ia[1], &
		    ja[1], &a[1], isym);
	} else if (*iter == 0) {
/*         Copy X to XL to check if initial guess is good enough. */
	    scopy_(n, &x[1], &c__1, &xl[1], &c__1);
	} else {
/*         Return since this is the first call to SPIGMR on a restart. */
	    return ret_val;
	}

	if (*jscal == 0 || *jscal == 2) {
/*         err = ||x-TrueSolution||/||TrueSolution||(2-Norms). */
	    if (*iter == 0) {
		solnrm = snrm2_(n, sslblk_1.soln, &c__1);
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dz[i__] = xl[i__] - sslblk_1.soln[i__ - 1];
/* L30: */
	    }
	    *err = snrm2_(n, &dz[1], &c__1) / solnrm;
	} else {
	    if (*iter == 0) {
		solnrm = 0.f;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		    r__1 = sx[i__] * sslblk_1.soln[i__ - 1];
		    solnrm += r__1 * r__1;
/* L40: */
		}
		solnrm = sqrt(solnrm);
	    }
	    dxnrm = 0.f;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		r__1 = sx[i__] * (xl[i__] - sslblk_1.soln[i__ - 1]);
		dxnrm += r__1 * r__1;
/* L50: */
	    }
	    dxnrm = sqrt(dxnrm);
/*         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms). */
	    *err = dxnrm / solnrm;
	}
    }

    if (*iunit != 0) {
	if (*iter == 0) {
	    io___10.ciunit = *iunit;
	    s_wsfe(&io___10);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itol), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*maxl), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*kmp), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	io___11.ciunit = *iunit;
	s_wsfe(&io___11);
	do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	r__1 = *rnrm / *bnrm;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*err), (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (*err <= *tol) {
	ret_val = 1;
    }

    return ret_val;
/* ------------- LAST LINE OF ISSGMR FOLLOWS ---------------------------- */
} /* issgmr_ */

