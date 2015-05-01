/* somn.f -- translated by f2c (version 12.02.01).
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

static integer c__3 = 3;
static integer c__1 = 1;

/* DECK SOMN */
/* Subroutine */ int somn_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, S_fp matvec, S_fp 
	msolve, integer *nsave, integer *itol, real *tol, integer *itmax, 
	integer *iter, real *err, integer *ierr, integer *iunit, real *r__, 
	real *z__, real *p, real *ap, real *emap, real *dz, real *csav, real *
	rwork, integer *iwork)
{
    /* System generated locals */
    integer ap_dim1, ap_offset, emap_dim1, emap_offset, p_dim1, p_offset, 
	    i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static integer i__, k, l;
    static real ak;
    static integer ip;
    static real bkl;
    static integer ipo;
    static real bnrm;
    static integer lmax;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real fuzz, akden, aknum;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);
    extern doublereal r1mach_(integer *);
    extern integer issomn_(integer *, real *, real *, integer *, integer *, 
	    integer *, real *, integer *, S_fp, integer *, integer *, real *, 
	    integer *, integer *, real *, integer *, integer *, real *, real *
	    , real *, real *, real *, real *, real *, real *, integer *, real 
	    *, real *, real *);
    static real solnrm;

/* ***BEGIN PROLOGUE  SOMN */
/* ***PURPOSE  Preconditioned Orthomin Sparse Iterative Ax=b Solver. */
/*            Routine to solve a general linear system  Ax = b  using */
/*            the Preconditioned Orthomin method. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2A4, D2B4 */
/* ***TYPE      SINGLE PRECISION (SOMN-S, DOMN-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, NON-SYMMETRIC LINEAR SYSTEM, */
/*             ORTHOMIN, SLAP, SPARSE */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE, ITOL, ITMAX */
/*     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINED) */
/*     REAL     B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N) */
/*     REAL     P(N,0:NSAVE), AP(N,0:NSAVE), EMAP(N,0:NSAVE) */
/*     REAL     DZ(N), CSAV(NSAVE), RWORK(USER DEFINED) */
/*     EXTERNAL MATVEC, MSOLVE */

/*     CALL SOMN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, */
/*    $     NSAVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, */
/*    $     Z, P, AP, EMAP, DZ, CSAV, RWORK, IWORK) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand side vector. */
/* X      :INOUT    Real X(N). */
/*         On input X is your initial guess for solution vector. */
/*         On output X is the final approximate solution. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays contain the matrix data structure for A. */
/*         It could take any form.  See "Description", below, for more */
/*         details. */
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
/*         non-zeros in the SLAP IA, JA, A storage for the matrix A. */
/*         ISYM is a flag which, if non-zero, denotest that A is */
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
/*         ISYM are defined as above.  RWORK is a real array that can */
/*         be used to pass necessary preconditioning information and/or */
/*         workspace to MSOLVE.  IWORK is an integer work array for */
/*         the same purpose as RWORK. */
/* NSAVE  :IN       Integer. */
/*         Number of  direction vectors to save and orthogonalize */
/*         against.  NSAVE >= 0. */
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
/*             COMMON /SSLBLK/ SOLN( ) */
/*         If ITOL=11, iteration stops when the 2-norm of the difference */
/*         between the iterative approximation and the user-supplied */
/*         solution divided by the 2-norm of the user-supplied solution */
/*         is less than TOL.  Note that this requires the user to set up */
/*         the "COMMON /SSLBLK/ SOLN(LENGTH)" in the calling routine. */
/*         The routine with this declaration should be loaded before the */
/*         stop test so that the correct length is used by the loader. */
/*         This procedure is not standard Fortran and may not work */
/*         correctly on your system (although it has worked on every */
/*         system the authors have tried).  If ITOL is not 11 then this */
/*         common block is indeed standard Fortran. */
/* TOL    :INOUT    Real. */
/*         Convergence criterion, as described above.  (Reset if IERR=4.) */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX+1 if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL. */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*           IERR = 0 => All went well. */
/*           IERR = 1 => Insufficient space allocated for WORK or IWORK. */
/*           IERR = 2 => Method failed to converge in ITMAX steps. */
/*           IERR = 3 => Error in user input. */
/*                       Check input values of N, ITOL. */
/*           IERR = 4 => User error tolerance set too tight. */
/*                       Reset to 500*R1MACH(3).  Iteration proceeded. */
/*           IERR = 5 => Preconditioning matrix, M, is not positive */
/*                       definite.  (r,z) < 0. */
/*           IERR = 6 => Breakdown of method detected. */
/*                       (p,Ap) < epsilon**2. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* R      :WORK     Real R(N). */
/* Z      :WORK     Real Z(N). */
/* P      :WORK     Real P(N,0:NSAVE). */
/* AP     :WORK     Real AP(N,0:NSAVE). */
/* EMAP   :WORK     Real EMAP(N,0:NSAVE). */
/* DZ     :WORK     Real DZ(N). */
/* CSAV   :WORK     Real CSAV(NSAVE) */
/*         Real arrays used for workspace. */
/* RWORK  :WORK     Real RWORK(USER DEFINED). */
/*         Real array that can be used for workspace in MSOLVE. */
/* IWORK  :WORK     Integer IWORK(USER DEFINED). */
/*         Integer array that can be used for workspace in MSOLVE. */

/* *Description */
/*       This routine does  not care  what matrix data   structure is */
/*       used for  A and M.  It simply   calls  the MATVEC and MSOLVE */
/*       routines, with  the arguments as  described above.  The user */
/*       could write any type of structure and the appropriate MATVEC */
/*       and MSOLVE routines.  It is assumed  that A is stored in the */
/*       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is */
/*       stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP */
/*       routines SSDOMN and SSLUOM are examples of this procedure. */

/*       Two  examples  of  matrix  data structures  are the: 1) SLAP */
/*       Triad  format and 2) SLAP Column format. */

/*       =================== S L A P Triad format =================== */
/*       In  this   format only the  non-zeros are  stored.  They may */
/*       appear  in *ANY* order.   The user  supplies three arrays of */
/*       length NELT, where  NELT  is the number  of non-zeros in the */
/*       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero */
/*       the  user puts   the row  and  column index   of that matrix */
/*       element in the IA and JA arrays.  The  value of the non-zero */
/*       matrix  element is  placed in  the corresponding location of */
/*       the A  array.  This is  an extremely easy data  structure to */
/*       generate.  On  the other hand it  is  not too  efficient  on */
/*       vector  computers   for the  iterative  solution  of  linear */
/*       systems.  Hence, SLAP  changes this input  data structure to */
/*       the SLAP   Column  format for the  iteration (but   does not */
/*       change it back). */

/*       Here is an example of the  SLAP Triad   storage format for a */
/*       5x5 Matrix.  Recall that the entries may appear in any order. */

/*           5x5 Matrix      SLAP Triad format for 5x5 matrix on left. */
/*                              1  2  3  4  5  6  7  8  9 10 11 */
/*       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 */
/*       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 */
/*       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       =================== S L A P Column format ================== */

/*       In  this format   the non-zeros are    stored counting  down */
/*       columns (except  for the diagonal  entry, which must  appear */
/*       first in each "column") and are  stored in the real array A. */
/*       In other words,  for  each column    in the matrix   put the */
/*       diagonal  entry  in A.   Then   put  in the  other  non-zero */
/*       elements going   down the  column (except  the  diagonal) in */
/*       order.  The IA array holds the row index  for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning   of   each  column.      That is,   IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)  points to the */
/*       end of the ICOL-th column.  Note that we always have JA(N+1) */
/*       = NELT+1, where N is the number of columns in the matrix and */
/*       NELT is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */

/* ***SEE ALSO  SSDOMN, SSLUOM, ISSOMN */
/* ***REFERENCES  1. Mark K. Seager, A SLAP for the Masses, in */
/*                  G. F. Carey, Ed., Parallel Supercomputing: Methods, */
/*                  Algorithms and Applications, Wiley, 1989, pp.135-155. */
/* ***ROUTINES CALLED  ISSOMN, R1MACH, SAXPY, SCOPY, SDOT */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871119  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890921  Removed TeX from comments.  (FNF) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   891004  Added new reference. */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF) */
/*   920407  COMMON BLOCK renamed SSLBLK.  (WRB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   920929  Corrected format of reference.  (FNF) */
/*   921019  Changed 500.0 to 500 to reduce SP/DP differences.  (FNF) */
/*   921113  Corrected C***CATEGORY line.  (FNF) */
/*   930326  Removed unused variable.  (FNF) */
/* ***END PROLOGUE  SOMN */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Subroutine Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/* ***FIRST EXECUTABLE STATEMENT  SOMN */

/*         Check some of the input data. */

    /* Parameter adjustments */
    --dz;
    --z__;
    --r__;
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --csav;
    emap_dim1 = *n;
    emap_offset = 1 + emap_dim1 * 0;
    emap -= emap_offset;
    ap_dim1 = *n;
    ap_offset = 1 + ap_dim1 * 0;
    ap -= ap_offset;
    p_dim1 = *n;
    p_offset = 1 + p_dim1 * 0;
    p -= p_offset;
    --rwork;
    --iwork;

    /* Function Body */
    *iter = 0;
    *ierr = 0;
    if (*n < 1) {
	*ierr = 3;
	return 0;
    }
    fuzz = r1mach_(&c__3);
    if (*tol < fuzz * 500) {
	*tol = fuzz * 500;
	*ierr = 4;
    }
    fuzz *= fuzz;

/*         Calculate initial residual and pseudo-residual, and check */
/*         stopping criterion. */
    (*matvec)(n, &x[1], &r__[1], nelt, &ia[1], &ja[1], &a[1], isym);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = b[i__] - r__[i__];
/* L10: */
    }
    (*msolve)(n, &r__[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &rwork[
	    1], &iwork[1]);

    if (issomn_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)
	    msolve, nsave, itol, tol, itmax, iter, err, ierr, iunit, &r__[1], 
	    &z__[1], &p[p_offset], &ap[ap_offset], &emap[emap_offset], &dz[1],
	     &csav[1], &rwork[1], &iwork[1], &ak, &bnrm, &solnrm) != 0) {
	goto L200;
    }
    if (*ierr != 0) {
	return 0;
    }


/*         ***** iteration loop ***** */

/* VD$R NOVECTOR */
/* VD$R NOCONCUR */
    i__1 = *itmax;
    for (k = 1; k <= i__1; ++k) {
	*iter = k;
	ip = (*iter - 1) % (*nsave + 1);

/*         calculate direction vector p, a*p, and (m-inv)*a*p, */
/*         and save if desired. */
	scopy_(n, &z__[1], &c__1, &p[ip * p_dim1 + 1], &c__1);
	(*matvec)(n, &p[ip * p_dim1 + 1], &ap[ip * ap_dim1 + 1], nelt, &ia[1],
		 &ja[1], &a[1], isym);
	(*msolve)(n, &ap[ip * ap_dim1 + 1], &emap[ip * emap_dim1 + 1], nelt, &
		ia[1], &ja[1], &a[1], isym, &rwork[1], &iwork[1]);
	if (*nsave == 0) {
	    akden = sdot_(n, &emap[emap_offset], &c__1, &emap[emap_offset], &
		    c__1);
	} else {
	    if (*iter > 1) {
/* Computing MIN */
		i__2 = *nsave, i__3 = *iter - 1;
		lmax = min(i__2,i__3);
		i__2 = lmax;
		for (l = 1; l <= i__2; ++l) {
		    ipo = (ip + (*nsave + 1 - l)) % (*nsave + 1);
		    bkl = sdot_(n, &emap[ip * emap_dim1 + 1], &c__1, &emap[
			    ipo * emap_dim1 + 1], &c__1);
		    bkl *= csav[l];
		    r__1 = -bkl;
		    saxpy_(n, &r__1, &p[ipo * p_dim1 + 1], &c__1, &p[ip * 
			    p_dim1 + 1], &c__1);
		    r__1 = -bkl;
		    saxpy_(n, &r__1, &ap[ipo * ap_dim1 + 1], &c__1, &ap[ip * 
			    ap_dim1 + 1], &c__1);
		    r__1 = -bkl;
		    saxpy_(n, &r__1, &emap[ipo * emap_dim1 + 1], &c__1, &emap[
			    ip * emap_dim1 + 1], &c__1);
/* L20: */
		}
		if (*nsave > 1) {
		    for (l = *nsave - 1; l >= 1; --l) {
			csav[l + 1] = csav[l];
/* L30: */
		    }
		}
	    }
	    akden = sdot_(n, &emap[ip * emap_dim1 + 1], &c__1, &emap[ip * 
		    emap_dim1 + 1], &c__1);
	    if (dabs(akden) < fuzz) {
		*ierr = 6;
		return 0;
	    }
	    csav[1] = 1.f / akden;

/*         calculate coefficient ak, new iterate x, new residual r, and */
/*         new pseudo-residual z. */
	}
	aknum = sdot_(n, &z__[1], &c__1, &emap[ip * emap_dim1 + 1], &c__1);
	ak = aknum / akden;
	saxpy_(n, &ak, &p[ip * p_dim1 + 1], &c__1, &x[1], &c__1);
	r__1 = -ak;
	saxpy_(n, &r__1, &ap[ip * ap_dim1 + 1], &c__1, &r__[1], &c__1);
	r__1 = -ak;
	saxpy_(n, &r__1, &emap[ip * emap_dim1 + 1], &c__1, &z__[1], &c__1);

/*         check stopping criterion. */
	if (issomn_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)
		msolve, nsave, itol, tol, itmax, iter, err, ierr, iunit, &r__[
		1], &z__[1], &p[p_offset], &ap[ap_offset], &emap[emap_offset],
		 &dz[1], &csav[1], &rwork[1], &iwork[1], &ak, &bnrm, &solnrm) 
		!= 0) {
	    goto L200;
	}

/* L100: */
    }

/*         *****   end of loop  ***** */

/*         Stopping criterion not satisfied. */
    *iter = *itmax + 1;
    *ierr = 2;

L200:
    return 0;
/* ------------- LAST LINE OF SOMN FOLLOWS ---------------------------- */
} /* somn_ */

