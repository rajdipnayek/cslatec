/* isdcgn.f -- translated by f2c (version 12.02.01).
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
    doublereal soln[1];
} dslblk_;

#define dslblk_1 dslblk_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* DECK ISDCGN */
integer isdcgn_(integer *n, doublereal *b, doublereal *x, integer *nelt, 
	integer *ia, integer *ja, doublereal *a, integer *isym, U_fp matvec, 
	S_fp mttvec, S_fp msolve, integer *itol, doublereal *tol, integer *
	itmax, integer *iter, doublereal *err, integer *ierr, integer *iunit, 
	doublereal *r__, doublereal *z__, doublereal *p, doublereal *atp, 
	doublereal *atz, doublereal *dz, doublereal *atdz, doublereal *rwork, 
	integer *iwork, doublereal *ak, doublereal *bk, doublereal *bnrm, 
	doublereal *solnrm)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 PCG Applied to the Normal Equations for"
	    " \002,\002N, ITOL = \002,i5,i5,/\002 ITER\002,\002   Error Estim"
	    "ate\002,\002            Alpha\002,\002             Beta\002)";
    static char fmt_1010[] = "(1x,i4,1x,d16.7,1x,d16.7,1x,d16.7)";

    /* System generated locals */
    integer ret_val, i__1;

    /* Local variables */
    static integer i__;
    extern doublereal dnrm2_(integer *, doublereal *, integer *), d1mach_(
	    integer *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_1010, 0 };


/* ***BEGIN PROLOGUE  ISDCGN */
/* ***SUBSIDIARY */
/* ***PURPOSE  Preconditioned CG on Normal Equations Stop Test. */
/*            This routine calculates the stop test for the Conjugate */
/*            Gradient iteration scheme applied to the normal equations. */
/*            It returns a non-zero if the error estimate (the type of */
/*            which is determined by ITOL) is less than the user */
/*            specified tolerance TOL. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      DOUBLE PRECISION (ISSCGN-S, ISDCGN-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM, */
/*             NORMAL EQUATIONS, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER */
/*     INTEGER  IERR, IUNIT, IWORK(USER DEFINED) */
/*     DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), Z(N), P(N) */
/*     DOUBLE PRECISION ATP(N), ATZ(N), DZ(N), ATDZ(N) */
/*     DOUBLE PRECISION RWORK(USER DEFINED), AK, BK, BNRM, SOLNRM */
/*     EXTERNAL MATVEC, MTTVEC, MSOLVE */

/*     IF( ISTPCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, */
/*    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, */
/*    $     ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM) */
/*    $     .NE. 0 ) THEN ITERATION DONE */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Double Precision B(N). */
/*         Right-hand side vector. */
/* X      :IN       Double Precision X(N). */
/*         The current approximate solution vector. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Double Precision A(NELT). */
/*         These arrays contain the matrix data structure for A. */
/*         It could take any form.  See "Description" in the */
/*         DCGN routine. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* MATVEC :EXT      External. */
/*         Name of a routine which performs the matrix vector multiply */
/*         Y = A*X given A and X.  The name of the MATVEC routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MATVEC is: */
/*             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM ) */
/*         Where N is the number of unknowns, Y is the product A*X */
/*         upon return X is an input vector, NELT is the number of */
/*         non-zeros in the SLAP-Column IA, JA, A storage for the matrix */
/*         A.  ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MTTVEC :EXT      External. */
/*         Name of a routine which performs the matrix transpose vector */
/*         multiply y = A'*X given A and X (where ' denotes transpose). */
/*         The name of the MTTVEC routine must be declared external in */
/*         the calling program.  The calling sequence to MTTVEC is the */
/*         same as that for MATVEC, viz.: */
/*             CALL MTTVEC( N, X, Y, NELT, IA, JA, A, ISYM ) */
/*         Where N is the number of unknowns, Y is the product A'*X */
/*         upon return X is an input vector, NELT is the number of */
/*         non-zeros in the SLAP-Column IA, JA, A storage for the matrix */
/*         A.  ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MSOLVE :EXT      External. */
/*         Name of a routine which solves a linear system MZ = R for */
/*         Z given R with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays).  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSOLVE is: */
/*             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) */
/*         Where N is the number of unknowns, R is the right-hand side */
/*         vector and Z is the solution upon return.  NELT, IA, JA, A and */
/*         ISYM are defined as above.  RWORK is a double precision array */
/*         that can be used to pass necessary preconditioning information */
/*         and/or workspace to MSOLVE.  IWORK is an integer work array */
/*         for the same purpose as RWORK. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate type of convergence criterion. */
/*         If ITOL=1, iteration stops when the 2-norm of the residual */
/*         divided by the 2-norm of the right-hand side is less than TOL. */
/*         If ITOL=2, iteration stops when the 2-norm of M-inv times the */
/*         residual divided by the 2-norm of M-inv times the right hand */
/*         side is less than TOL, where M-inv is the inverse of the */
/*         diagonal of A. */
/*         ITOL=11 is often useful for checking and comparing different */
/*         routines.  For this case, the user must supply the "exact" */
/*         solution or a very accurate approximation (one with an error */
/*         much less than TOL) through a common block, */
/*             COMMON /DSLBLK/ SOLN( ) */
/*         If ITOL=11, iteration stops when the 2-norm of the difference */
/*         between the iterative approximation and the user-supplied */
/*         solution divided by the 2-norm of the user-supplied solution */
/*         is less than TOL.  Note that this requires the user to set up */
/*         the "COMMON /DSLBLK/ SOLN(LENGTH)" in the calling routine. */
/*         The routine with this declaration should be loaded before the */
/*         stop test so that the correct length is used by the loader. */
/*         This procedure is not standard Fortran and may not work */
/*         correctly on your system (although it has worked on every */
/*         system the authors have tried).  If ITOL is not 11 then this */
/*         common block is indeed standard Fortran. */
/* TOL    :IN       Double Precision. */
/*         Convergence criterion, as described above. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :IN       Integer. */
/*         Current iteration count.  (Must be zero on first call.) */
/* ERR    :OUT      Double Precision. */
/*         Error estimate of error in the X(N) approximate solution, as */
/*         defined by ITOL. */
/* IERR   :OUT      Integer. */
/*         Error flag.  IERR is set to 3 if ITOL is not one of the */
/*         acceptable values, see above. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* R      :IN       Double Precision R(N). */
/*         The residual R = B-AX. */
/* Z      :WORK     Double Precision Z(N). */
/*         Double Precision array used for workspace. */
/* P      :IN       Double Precision P(N). */
/*         The conjugate direction vector. */
/* ATP    :IN       Double Precision ATP(N). */
/*         A-transpose times the conjugate direction vector. */
/* ATZ    :IN       Double Precision ATZ(N). */
/*         A-transpose times the pseudo-residual. */
/* DZ     :IN       Double Precision DZ(N). */
/*         Workspace used to hold temporary vector(s). */
/* ATDZ   :WORK     Double Precision ATDZ(N). */
/*         Workspace. */
/* RWORK  :WORK     Double Precision RWORK(USER DEFINED). */
/*         Double Precision array that can be used by MSOLVE. */
/* IWORK  :WORK     Integer IWORK(USER DEFINED). */
/*         Integer array that can be used by MSOLVE. */
/* AK     :IN       Double Precision. */
/* BK     :IN       Double Precision. */
/*         Current conjugate gradient parameters alpha and beta. */
/* BNRM   :INOUT    Double Precision. */
/*         Norm of the right hand side.  Type of norm depends on ITOL. */
/*         Calculated only on the first call. */
/* SOLNRM :INOUT    Double Precision. */
/*         2-Norm of the true solution, SOLN.  Only computed and used */
/*         if ITOL = 11. */

/* *Function Return Values: */
/*       0 : Error estimate (determined by ITOL) is *NOT* less than the */
/*           specified tolerance, TOL.  The iteration must continue. */
/*       1 : Error estimate (determined by ITOL) is less than the */
/*           specified tolerance, TOL.  The iteration can be considered */
/*           complete. */

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */

/* ***SEE ALSO  DCGN */
/* ***ROUTINES CALLED  D1MACH, DNRM2 */
/* ***COMMON BLOCKS    DSLBLK */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890404  DATE WRITTEN */
/*   890404  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   891003  Removed C***REFER TO line, per MKS. */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Removed MATVEC, MTTVEC and MSOLVE from ROUTINES CALLED */
/*           list.  (FNF) */
/*   910506  Made subsidiary to DCGN.  (FNF) */
/*   920407  COMMON BLOCK renamed DSLBLK.  (WRB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   920930  Corrected to not print AK,BK when ITER=0.  (FNF) */
/*   921026  Changed 1.0E10 to D1MACH(2) and corrected D to E in */
/*           output format.  (FNF) */
/*   921113  Corrected C***CATEGORY line.  (FNF) */
/* ***END PROLOGUE  ISDCGN */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Arrays in Common .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Common blocks .. */
/* ***FIRST EXECUTABLE STATEMENT  ISDCGN */
    /* Parameter adjustments */
    --atdz;
    --dz;
    --atz;
    --atp;
    --p;
    --z__;
    --r__;
    --a;
    --x;
    --b;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    ret_val = 0;

    if (*itol == 1) {
/*         err = ||Residual||/||RightHandSide|| (2-Norms). */
	if (*iter == 0) {
	    *bnrm = dnrm2_(n, &b[1], &c__1);
	}
	*err = dnrm2_(n, &r__[1], &c__1) / *bnrm;
    } else if (*itol == 2) {
/*                  -1              -1 */
/*         err = ||M  Residual||/||M  RightHandSide|| (2-Norms). */
	if (*iter == 0) {
	    (*msolve)(n, &b[1], &dz[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		    rwork[1], &iwork[1]);
	    (*mttvec)(n, &dz[1], &atdz[1], nelt, &ia[1], &ja[1], &a[1], isym);
	    *bnrm = dnrm2_(n, &atdz[1], &c__1);
	}
	*err = dnrm2_(n, &atz[1], &c__1) / *bnrm;
    } else if (*itol == 11) {
/*         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms). */
	if (*iter == 0) {
	    *solnrm = dnrm2_(n, dslblk_1.soln, &c__1);
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dz[i__] = x[i__] - dslblk_1.soln[i__ - 1];
/* L10: */
	}
	*err = dnrm2_(n, &dz[1], &c__1) / *solnrm;
    } else {

/*         If we get here ITOL is not one of the acceptable values. */
	*err = d1mach_(&c__2);
	*ierr = 3;
    }

    if (*iunit != 0) {
	if (*iter == 0) {
	    io___2.ciunit = *iunit;
	    s_wsfe(&io___2);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itol), (ftnlen)sizeof(integer));
	    e_wsfe();
	    io___3.ciunit = *iunit;
	    s_wsfe(&io___3);
	    do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*err), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	} else {
	    io___4.ciunit = *iunit;
	    s_wsfe(&io___4);
	    do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*err), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*ak), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*bk), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
    }
    if (*err <= *tol) {
	ret_val = 1;
    }

    return ret_val;
/* ------------- LAST LINE OF ISDCGN FOLLOWS ---------------------------- */
} /* isdcgn_ */

