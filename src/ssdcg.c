/* ssdcg.f -- translated by f2c (version 12.02.01).
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

/* DECK SSDCG */
/* Subroutine */ int ssdcg_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *itol, real 
	*tol, integer *itmax, integer *iter, real *err, integer *ierr, 
	integer *iunit, real *rwork, integer *lenw, integer *iwork, integer *
	leniw)
{
    extern /* Subroutine */ int scg_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, U_fp, U_fp, integer *, 
	    real *, integer *, integer *, real *, integer *, integer *, real *
	    , real *, real *, real *, real *, integer *), ss2y_(integer *, 
	    integer *, integer *, integer *, real *, integer *);
    static integer locd, locp, locr, locw;
    extern /* Subroutine */ int ssdi_();
    static integer locz;
    extern /* Subroutine */ int ssds_(integer *, integer *, integer *, 
	    integer *, real *, integer *, real *);
    extern /* Subroutine */ int ssmv_();
    static integer locdz, lociw;
    extern /* Subroutine */ int schkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, ftnlen);

/* ***BEGIN PROLOGUE  SSDCG */
/* ***PURPOSE  Diagonally Scaled Conjugate Gradient Sparse Ax=b Solver. */
/*            Routine to solve a symmetric positive definite linear */
/*            system  Ax = b  using the Preconditioned Conjugate */
/*            Gradient method.  The preconditioner is diagonal scaling. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2B4 */
/* ***TYPE      SINGLE PRECISION (SSDCG-S, DSDCG-D) */
/* ***KEYWORDS  ITERATIVE PRECONDITION, SLAP, SPARSE, */
/*             SYMMETRIC LINEAR SYSTEM */
/* ***AUTHOR  Greenbaum, Anne, (Courant Institute) */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-60 */
/*             Livermore, CA 94550 (510) 423-3141 */
/*             seager@llnl.gov */
/* ***DESCRIPTION */

/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*     INTEGER ITER, IERR, IUNIT, LENW, IWORK(10), LENIW */
/*     REAL B(N), X(N), A(NELT), TOL, ERR, RWORK(5*N) */

/*     CALL SSDCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, */
/*    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW ) */

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
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Real A(NELT). */
/*         These arrays should hold the matrix A in either the SLAP */
/*         Triad format or the SLAP Column format.  See "Description", */
/*         below.  If the SLAP Triad format is chosen it is changed */
/*         internally to the SLAP Column format. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all non-zero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
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
/*           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* RWORK  :WORK     Real RWORK(LENW). */
/*         Real array used for workspace. */
/* LENW   :IN       Integer. */
/*         Length of the real workspace, RWORK.  LENW >= 5*N. */
/* IWORK  :WORK     Integer IWORK(LENIW). */
/*         Used to hold pointers into the real workspace, RWORK. */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Real workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace, IWORK.  LENIW >= 10. */

/* *Description: */
/*       This  routine   performs preconditioned conjugate   gradient */
/*       method on  the  symmetric positive definite   linear  system */
/*       Ax=b.   The preconditioner is  M = DIAG(A), the  diagonal of */
/*       the matrix A.  This is the  simplest of preconditioners  and */
/*       vectorizes very well.   This routine is  simply a driver for */
/*       the SCG routine.  It  calls the SSDS  routine to  set up the */
/*       preconditioning  and  then  calls  SCG  with the appropriate */
/*       MATVEC and MSOLVE routines. */

/*       The Sparse Linear Algebra Package (SLAP) utilizes two matrix */
/*       data structures: 1) the  SLAP Triad  format or  2)  the SLAP */
/*       Column format.  The user can hand this routine either of the */
/*       of these data structures and SLAP  will figure out  which on */
/*       is being used and act accordingly. */

/*       =================== S L A P Triad format =================== */

/*       This routine requires that the  matrix A be   stored in  the */
/*       SLAP  Triad format.  In  this format only the non-zeros  are */
/*       stored.  They may appear in  *ANY* order.  The user supplies */
/*       three arrays of  length NELT, where  NELT is  the number  of */
/*       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For */
/*       each non-zero the user puts the row and column index of that */
/*       matrix element  in the IA and  JA arrays.  The  value of the */
/*       non-zero   matrix  element is  placed  in  the corresponding */
/*       location of the A array.   This is  an  extremely  easy data */
/*       structure to generate.  On  the  other hand it   is  not too */
/*       efficient on vector computers for  the iterative solution of */
/*       linear systems.  Hence,   SLAP changes   this  input    data */
/*       structure to the SLAP Column format  for  the iteration (but */
/*       does not change it back). */

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

/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except for  the diagonal entry, which */
/*       must appear first in each  "column")  and are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

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


/* *Side Effects: */
/*       The SLAP Triad format (IA, JA, A) is modified internally  to */
/*       be the SLAP Column format.  See above. */

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */

/* ***SEE ALSO  SCG, SSICCG */
/* ***REFERENCES  1. Louis Hageman and David Young, Applied Iterative */
/*                  Methods, Academic Press, New York, 1981. */
/*               2. Concus, Golub and O'Leary, A Generalized Conjugate */
/*                  Gradient Method for the Numerical Solution of */
/*                  Elliptic Partial Differential Equations, in Sparse */
/*                  Matrix Computations, Bunch and Rose, Eds., Academic */
/*                  Press, New York, 1979. */
/* ***ROUTINES CALLED  SCG, SCHKW, SS2Y, SSDI, SSDS, SSMV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   871119  DATE WRITTEN */
/*   881213  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890921  Removed TeX from comments.  (FNF) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   920407  COMMON BLOCK renamed SSLBLK.  (WRB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   920929  Corrected format of references.  (FNF) */
/* ***END PROLOGUE  SSDCG */
/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  SSDCG */

    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    *ierr = 0;
    if (*n < 1 || *nelt < 1) {
	*ierr = 3;
	return 0;
    }

/*         Modify the SLAP matrix data structure to YSMP-Column. */
    ss2y_(n, nelt, &ia[1], &ja[1], &a[1], isym);

/*         Set up the work arrays. */
    lociw = 11;

    locd = 1;
    locr = locd + *n;
    locz = locr + *n;
    locp = locz + *n;
    locdz = locp + *n;
    locw = locdz + *n;

/*         Check the workspace allocations. */
    schkw_("SSDCG", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)5);
    if (*ierr != 0) {
	return 0;
    }

    iwork[4] = locd;
    iwork[9] = lociw;
    iwork[10] = locw;

/*         Compute the inverse of the diagonal of the matrix.  This */
/*         will be used as the preconditioner. */
    ssds_(n, nelt, &ia[1], &ja[1], &a[1], isym, &rwork[locd]);

/*         Do the Preconditioned Conjugate Gradient. */
    scg_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (U_fp)ssmv_, (
	    U_fp)ssdi_, itol, tol, itmax, iter, err, ierr, iunit, &rwork[locr]
	    , &rwork[locz], &rwork[locp], &rwork[locdz], &rwork[1], &iwork[1])
	    ;
    return 0;
/* ------------- LAST LINE OF SSDCG FOLLOWS ----------------------------- */
} /* ssdcg_ */

