/* dsiccg.f -- translated by f2c (version 12.02.01).
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
static integer c__3 = 3;

/* DECK DSICCG */
/* Subroutine */ int dsiccg_(integer *n, doublereal *b, doublereal *x, 
	integer *nelt, integer *ia, integer *ja, doublereal *a, integer *isym,
	 integer *itol, doublereal *tol, integer *itmax, integer *iter, 
	doublereal *err, integer *ierr, integer *iunit, doublereal *rwork, 
	integer *lenw, integer *iwork, integer *leniw)
{
    /* System generated locals */
    address a__1[3];
    integer i__1[3];
    char ch__1[101];

    /* Local variables */
    static integer nl;
    extern /* Subroutine */ int dcg_(integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, U_fp, 
	    U_fp, integer *, doublereal *, integer *, integer *, doublereal *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), ds2y_(integer *, integer *
	    , integer *, integer *, doublereal *, integer *);
    static integer locp, locr, locw, locz;
    extern /* Subroutine */ int dsmv_();
    static char xern1[8];
    static integer locel;
    extern /* Subroutine */ int dchkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, doublereal *, ftnlen), dsics_(
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *);
    static integer locdz, lociw, lociel, locdin, locjel;
    extern /* Subroutine */ int dsllti_();
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

    /* Fortran I/O blocks */
    static icilist io___13 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  DSICCG */
/* ***PURPOSE  Incomplete Cholesky Conjugate Gradient Sparse Ax=b Solver. */
/*            Routine to solve a symmetric positive definite linear */
/*            system  Ax = b  using the incomplete Cholesky */
/*            Preconditioned Conjugate Gradient method. */
/* ***LIBRARY   SLATEC (SLAP) */
/* ***CATEGORY  D2B4 */
/* ***TYPE      DOUBLE PRECISION (SSICCG-S, DSICCG-D) */
/* ***KEYWORDS  INCOMPLETE CHOLESKY, ITERATIVE PRECONDITION, SLAP, SPARSE, */
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
/*     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NL+2*N+1), LENIW */
/*     DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(NL+5*N) */

/*     CALL DSICCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, */
/*    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Double Precision B(N). */
/*         Right-hand side vector. */
/* X      :INOUT    Double Precision X(N). */
/*         On input X is your initial guess for solution vector. */
/*         On output X is the final approximate solution. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Double Precision A(NELT). */
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
/* TOL    :INOUT    Double Precision. */
/*         Convergence criterion, as described above.  (Reset if IERR=4.) */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX+1 if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Double Precision. */
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
/*                       Reset to 500*D1MACH(3).  Iteration proceeded. */
/*           IERR = 5 => Preconditioning matrix, M, is not positive */
/*                       definite.  (r,z) < 0. */
/*           IERR = 6 => Matrix A is not positive definite.  (p,Ap) < 0. */
/*           IERR = 7 => Incomplete factorization broke down and was */
/*                       fudged.  Resulting preconditioning may be less */
/*                       than the best. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* RWORK  :WORK     Double Precision RWORK(LENW). */
/*         Double Precision array used for workspace. */
/* LENW   :IN       Integer. */
/*         Length of the double precision workspace, RWORK. */
/*         LENW >= NL+5*N. */
/*         NL is the number of non-zeros in the lower triangle of the */
/*         matrix (including the diagonal). */
/* IWORK  :WORK     Integer IWORK(LENIW). */
/*         Integer array used for workspace. */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Double Precision workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace, IWORK.  LENIW >= NL+N+11. */
/*         NL is the number of non-zeros in the lower triangle of the */
/*         matrix (including the diagonal). */

/* *Description: */
/*       This routine  performs  preconditioned  conjugate   gradient */
/*       method on the   symmetric positive  definite  linear  system */
/*       Ax=b.   The preconditioner  is  the incomplete Cholesky (IC) */
/*       factorization of the matrix A.  See  DSICS for details about */
/*       the incomplete   factorization algorithm.  One   should note */
/*       here however, that the  IC factorization is a  slow  process */
/*       and  that  one should   save  factorizations  for  reuse, if */
/*       possible.  The   MSOLVE operation (handled  in  DSLLTI) does */
/*       vectorize on machines  with  hardware  gather/scatter and is */
/*       quite fast. */

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
/*       double precision array A.   In other words,  for each column */
/*       in the matrix put the diagonal entry in  A.  Then put in the */
/*       other non-zero  elements going down  the column (except  the */
/*       diagonal) in order.   The  IA array holds the  row index for */
/*       each non-zero.  The JA array holds the offsets  into the IA, */
/*       A arrays  for  the  beginning  of each   column.   That  is, */
/*       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the */
/*       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1), */
/*       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column. */
/*       Note that we always have  JA(N+1) = NELT+1,  where N is  the */
/*       number of columns in  the matrix and NELT  is the number  of */
/*       non-zeros in the matrix. */

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
/*       The SLAP Triad format (IA, JA, A) is modified internally to be */
/*       the SLAP Column format.  See above. */

/* *Cautions: */
/*     This routine will attempt to write to the Fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this logical unit is attached to a file or terminal before calling */
/*     this routine with a non-zero value for IUNIT.  This routine does */
/*     not check for the validity of a non-zero IUNIT unit number. */

/* ***SEE ALSO  DCG, DSLLTI */
/* ***REFERENCES  1. Louis Hageman and David Young, Applied Iterative */
/*                  Methods, Academic Press, New York, 1981. */
/*               2. Concus, Golub and O'Leary, A Generalized Conjugate */
/*                  Gradient Method for the Numerical Solution of */
/*                  Elliptic Partial Differential Equations, in Sparse */
/*                  Matrix Computations, Bunch and Rose, Eds., Academic */
/*                  Press, New York, 1979. */
/* ***ROUTINES CALLED  DCG, DCHKW, DS2Y, DSICS, DSLLTI, DSMV, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   890404  DATE WRITTEN */
/*   890404  Previous REVISION DATE */
/*   890915  Made changes requested at July 1989 CML Meeting.  (MKS) */
/*   890921  Removed TeX from comments.  (FNF) */
/*   890922  Numerous changes to prologue to make closer to SLATEC */
/*           standard.  (FNF) */
/*   890929  Numerous changes to reduce SP/DP differences.  (FNF) */
/*   900805  Changed XERRWV calls to calls to XERMSG.  (RWC) */
/*   910411  Prologue converted to Version 4.0 format.  (BAB) */
/*   920407  COMMON BLOCK renamed DSLBLK.  (WRB) */
/*   920511  Added complete declaration section.  (WRB) */
/*   920929  Corrected format of references.  (FNF) */
/*   921019  Corrected NEL to NL.  (FNF) */
/* ***END PROLOGUE  DSICCG */
/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/* ***FIRST EXECUTABLE STATEMENT  DSICCG */

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

/*         Change the SLAP input matrix IA, JA, A to SLAP-Column format. */
    ds2y_(n, nelt, &ia[1], &ja[1], &a[1], isym);

/*         Count number of elements in lower triangle of the matrix. */
/*         Then set up the work arrays. */
    if (*isym == 0) {
	nl = (*nelt + *n) / 2;
    } else {
	nl = *nelt;
    }

    locjel = 11;
    lociel = locjel + nl;
    lociw = lociel + *n + 1;

    locel = 1;
    locdin = locel + nl;
    locr = locdin + *n;
    locz = locr + *n;
    locp = locz + *n;
    locdz = locp + *n;
    locw = locdz + *n;

/*         Check the workspace allocations. */
    dchkw_("DSICCG", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)6);
    if (*ierr != 0) {
	return 0;
    }

    iwork[1] = nl;
    iwork[2] = locjel;
    iwork[3] = lociel;
    iwork[4] = locel;
    iwork[5] = locdin;
    iwork[9] = lociw;
    iwork[10] = locw;

/*         Compute the Incomplete Cholesky decomposition. */

    dsics_(n, nelt, &ia[1], &ja[1], &a[1], isym, &nl, &iwork[lociel], &iwork[
	    locjel], &rwork[locel], &rwork[locdin], &rwork[locr], ierr);
    if (*ierr != 0) {
	s_wsfi(&io___13);
	do_fio(&c__1, (char *)&(*ierr), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 36, a__1[0] = "IC factorization broke down on step ";
	i__1[1] = 8, a__1[1] = xern1;
	i__1[2] = 57, a__1[2] = ".  Diagonal was set to unity and factorizat"
		"ion proceeded.";
	s_cat(ch__1, a__1, i__1, &c__3, (ftnlen)101);
	xermsg_("SLATEC", "DSICCG", ch__1, &c__1, &c__1, (ftnlen)6, (ftnlen)6,
		 (ftnlen)101);
	*ierr = 7;
    }

/*         Do the Preconditioned Conjugate Gradient. */
    dcg_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (U_fp)dsmv_, (
	    U_fp)dsllti_, itol, tol, itmax, iter, err, ierr, iunit, &rwork[
	    locr], &rwork[locz], &rwork[locp], &rwork[locdz], &rwork[1], &
	    iwork[1]);
    return 0;
/* ------------- LAST LINE OF DSICCG FOLLOWS ---------------------------- */
} /* dsiccg_ */

