/* sos.f -- translated by f2c (version 12.02.01).
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
static integer c__2 = 2;
static integer c__4 = 4;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__6 = 6;

/* DECK SOS */
/* Subroutine */ int sos_(U_fp fnc, integer *neq, real *x, real *rtolx, real *
	atolx, real *tolf, integer *iflag, real *rw, integer *lrw, integer *
	iw, integer *liw)
{
    /* System generated locals */
    address a__1[2], a__2[4];
    integer i__1[2], i__2[4];
    char ch__1[97], ch__2[151], ch__3[105], ch__4[179], ch__5[116], ch__6[100]
	    ;

    /* Local variables */
    static integer k1, k2, k3, k4, k5, k6, nc, ncjs, nsri, mxit;
    static char xern1[8], xern3[16], xern4[16];
    static integer nsrrc, inpflg;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static integer iprint;
    extern /* Subroutine */ int soseqs_(U_fp, integer *, real *, real *, real 
	    *, real *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, real *, real *, integer *, real *, real *, real *, 
	    real *, real *, real *, integer *);

    /* Fortran I/O blocks */
    static icilist io___3 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___5 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___7 = { 0, xern4, 0, "(1PE15.6)", 16, 1 };
    static icilist io___8 = { 0, xern3, 0, "(1PE15.6)", 16, 1 };
    static icilist io___11 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___13 = { 0, xern1, 0, "(I8)", 8, 1 };
    static icilist io___14 = { 0, xern1, 0, "(I8)", 8, 1 };


/* ***BEGIN PROLOGUE  SOS */
/* ***PURPOSE  Solve a square system of nonlinear equations. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F2A */
/* ***TYPE      SINGLE PRECISION (SOS-S, DSOS-D) */
/* ***KEYWORDS  BROWN'S METHOD, NEWTON'S METHOD, NONLINEAR EQUATIONS, */
/*             ROOTS, SOLUTIONS */
/* ***AUTHOR  Watts, H. A., (SNLA) */
/* ***DESCRIPTION */

/*     SOS solves a system of NEQ simultaneous nonlinear equations in */
/*     NEQ unknowns.  That is, it solves the problem   F(X)=0 */
/*     where X is a vector with components  X(1),...,X(NEQ)  and  F */
/*     is a vector of nonlinear functions.  Each equation is of the form */

/*               F (X(1),...,X(NEQ))=0     for K=1,...,NEQ. */
/*                K */

/*     The algorithm is based on an iterative method which is a */
/*     variation of Newton's method using Gaussian elimination */
/*     in a manner similar to the Gauss-Seidel process.  Convergence */
/*     is roughly quadratic.  All partial derivatives required by */
/*     the algorithm are approximated by first difference quotients. */
/*     The convergence behavior of this code is affected by the */
/*     ordering of the equations, and it is advantageous to place linear */
/*     and mildly nonlinear equations first in the ordering. */

/*     Actually, SOS is merely an interfacing routine for */
/*     calling subroutine SOSEQS which embodies the solution */
/*     algorithm.  The purpose of this is to add greater */
/*     flexibility and ease of use for the prospective user. */

/*     SOSEQS calls the accompanying routine SOSSOL, which solves special */
/*     triangular linear systems by back-substitution. */

/*     The user must supply a function subprogram which evaluates the */
/*     K-th equation only (K specified by SOSEQS) for each call */
/*     to the subprogram. */

/*     SOS represents an implementation of the mathematical algorithm */
/*     described in the references below.  It is a modification of the */
/*     code SOSNLE written by H. A. Watts in 1973. */

/* ********************************************************************** */
/*   -Input- */

/*     FNC -Name of the function program which evaluates the equations. */
/*          This name must be in an EXTERNAL statement in the calling */
/*          program.  The user must supply FNC in the form FNC(X,K), */
/*          where X is the solution vector (which must be dimensioned */
/*          in FNC) and FNC returns the value of the K-th function. */

/*     NEQ -Number of equations to be solved. */

/*     X   -Solution vector.  Initial guesses must be supplied. */

/*     RTOLX -Relative error tolerance used in the convergence criteria. */
/*          Each solution component X(I) is checked by an accuracy test */
/*          of the form   ABS(X(I)-XOLD(I)) .LE. RTOLX*ABS(X(I))+ATOLX, */
/*          where XOLD(I) represents the previous iteration value. */
/*          RTOLX must be non-negative. */

/*     ATOLX -Absolute error tolerance used in the convergence criteria. */
/*          ATOLX must be non-negative.  If the user suspects some */
/*          solution component may be zero, he should set ATOLX to an */
/*          appropriate (depends on the scale of the remaining variables) */
/*          positive value for better efficiency. */

/*     TOLF -Residual error tolerance used in the convergence criteria. */
/*          Convergence will be indicated if all residuals (values of the */
/*          functions or equations) are not bigger than TOLF in */
/*          magnitude.  Note that extreme care must be given in assigning */
/*          an appropriate value for TOLF because this convergence test */
/*          is dependent on the scaling of the equations.  An */
/*          inappropriate value can cause premature termination of the */
/*          iteration process. */

/*     IFLAG -Optional input indicator.  You must set  IFLAG=-1  if you */
/*          want to use any of the optional input items listed below. */
/*          Otherwise set it to zero. */

/*     RW  -A REAL work array which is split apart by SOS and used */
/*          internally by SOSEQS. */

/*     LRW -Dimension of the RW array.  LRW must be at least */
/*                    1 + 6*NEQ + NEQ*(NEQ+1)/2 */

/*     IW  -An INTEGER work array which is split apart by SOS and used */
/*          internally by SOSEQS. */

/*     LIW -Dimension of the IW array. LIW must be at least  3 + NEQ. */

/*   -Optional Input- */

/*     IW(1) -Internal printing parameter.  You must set  IW(1)=-1  if */
/*          you want the intermediate solution iterates to be printed. */

/*     IW(2) -Iteration limit.  The maximum number of allowable */
/*          iterations can be specified, if desired.  To override the */
/*          default value of 50, set IW(2) to the number wanted. */

/*     Remember, if you tell the code that you are using one of the */
/*               options (by setting IFLAG=-1), you must supply values */
/*               for both IW(1) and IW(2). */

/* ********************************************************************** */
/*   -Output- */

/*     X   -Solution vector. */

/*     IFLAG -Status indicator */

/*                         *** Convergence to a Solution *** */

/*          1 Means satisfactory convergence to a solution was achieved. */
/*            Each solution component X(I) satisfies the error tolerance */
/*            test   ABS(X(I)-XOLD(I)) .LE. RTOLX*ABS(X(I))+ATOLX. */

/*          2 Means procedure converged to a solution such that all */
/*            residuals are at most TOLF in magnitude, */
/*            ABS(FNC(X,I)) .LE. TOLF. */

/*          3 Means that conditions for both IFLAG=1 and IFLAG=2 hold. */

/*          4 Means possible numerical convergence.  Behavior indicates */
/*            limiting precision calculations as a result of user asking */
/*            for too much accuracy or else convergence is very slow. */
/*            Residual norms and solution increment norms have */
/*            remained roughly constant over several consecutive */
/*            iterations. */

/*                         *** Task Interrupted *** */

/*          5 Means the allowable number of iterations has been met */
/*            without obtaining a solution to the specified accuracy. */
/*            Very slow convergence may be indicated.  Examine the */
/*            approximate solution returned and see if the error */
/*            tolerances seem appropriate. */

/*          6 Means the allowable number of iterations has been met and */
/*            the iterative process does not appear to be converging. */
/*            A local minimum may have been encountered or there may be */
/*            limiting precision difficulties. */

/*          7 Means that the iterative scheme appears to be diverging. */
/*            Residual norms and solution increment norms have */
/*            increased over several consecutive iterations. */

/*                         *** Task Cannot Be Continued *** */

/*          8 Means that a Jacobian-related matrix was singular. */

/*          9 Means improper input parameters. */

/*          *** IFLAG should be examined after each call to   *** */
/*          *** SOS with the appropriate action being taken.  *** */


/*     RW(1) -Contains a norm of the residual. */

/*     IW(3) -Contains the number of iterations used by the process. */

/* ********************************************************************** */
/* ***REFERENCES  K. M. Brown, Solution of simultaneous nonlinear */
/*                 equations, Algorithm 316, Communications of the */
/*                 A.C.M. 10, (1967), pp. 728-729. */
/*               K. M. Brown, A quadratically convergent Newton-like */
/*                 method based upon Gaussian elimination, SIAM Journal */
/*                 on Numerical Analysis 6, (1969), pp. 560-569. */
/* ***ROUTINES CALLED  SOSEQS, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   801001  DATE WRITTEN */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900510  Convert XERRWV calls to XERMSG calls, changed Prologue */
/*           comments to agree with DSOS.  (RWC) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SOS */
/* ***FIRST EXECUTABLE STATEMENT  SOS */
    /* Parameter adjustments */
    --iw;
    --rw;
    --x;

    /* Function Body */
    inpflg = *iflag;

/*     CHECK FOR VALID INPUT */

    if (*neq <= 0) {
	s_wsfi(&io___3);
	do_fio(&c__1, (char *)&(*neq), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 89, a__1[0] = "THE NUMBER OF EQUATIONS MUST BE A POSITIVE "
		"INTEGER.  YOU HAVE CALLED THE CODE WITH NEQ = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)97);
	xermsg_("SLATEC", "SOS", ch__1, &c__1, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)97);
	*iflag = 9;
    }

    if (*rtolx < 0. || *atolx < 0.) {
	s_wsfi(&io___5);
	do_fio(&c__1, (char *)&(*atolx), (ftnlen)sizeof(real));
	e_wsfi();
	s_wsfi(&io___7);
	do_fio(&c__1, (char *)&(*rtolx), (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__2[0] = 106, a__2[0] = "THE ERROR TOLERANCES FOR THE SOLUTION ITER"
		"ATES CANNOT BE NEGATIVE. YOU HAVE CALLED THE CODE WITH  RTOL"
		"X = ";
	i__2[1] = 16, a__2[1] = xern3;
	i__2[2] = 13, a__2[2] = " AND ATOLX = ";
	i__2[3] = 16, a__2[3] = xern4;
	s_cat(ch__2, a__2, i__2, &c__4, (ftnlen)151);
	xermsg_("SLATEC", "SOS", ch__2, &c__2, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)151);
	*iflag = 9;
    }

    if (*tolf < 0.) {
	s_wsfi(&io___8);
	do_fio(&c__1, (char *)&(*tolf), (ftnlen)sizeof(real));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 89, a__1[0] = "THE RESIDUAL ERROR TOLERANCE MUST BE NON-NE"
		"GATIVE.  YOU HAVE CALLED THE CODE WITH TOLF = ";
	i__1[1] = 16, a__1[1] = xern3;
	s_cat(ch__3, a__1, i__1, &c__2, (ftnlen)105);
	xermsg_("SLATEC", "SOS", ch__3, &c__3, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)105);
	*iflag = 9;
    }

    iprint = 0;
    mxit = 50;
    if (inpflg == -1) {
	if (iw[1] == -1) {
	    iprint = -1;
	}
	mxit = iw[2];
	if (mxit <= 0) {
	    s_wsfi(&io___11);
	    do_fio(&c__1, (char *)&mxit, (ftnlen)sizeof(integer));
	    e_wsfi();
/* Writing concatenation */
	    i__1[0] = 171, a__1[0] = "YOU HAVE TOLD THE CODE TO USE OPTIONAL"
		    " IN PUT ITEMS BY SETTING  IFLAG=-1. HOWEVER YOU HAVE CAL"
		    "LED THE CODE WITH THE MAXIMUM ALLOWABLE NUMBER OF ITERAT"
		    "IONS SET TO  IW(2) = ";
	    i__1[1] = 8, a__1[1] = xern1;
	    s_cat(ch__4, a__1, i__1, &c__2, (ftnlen)179);
	    xermsg_("SLATEC", "SOS", ch__4, &c__4, &c__1, (ftnlen)6, (ftnlen)
		    3, (ftnlen)179);
	    *iflag = 9;
	}
    }

    nc = *neq * (*neq + 1) / 2;
    if (*lrw < *neq * 6 + 1 + nc) {
	s_wsfi(&io___13);
	do_fio(&c__1, (char *)&(*lrw), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 108, a__1[0] = "DIMENSION OF THE RW ARRAY MUST BE AT LEAST"
		" 1 + 6*NEQ + NEQ*(NEQ+1)/2 .  YOU HAVE CALLED THE CODE WITH "
		"LRW = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__5, a__1, i__1, &c__2, (ftnlen)116);
	xermsg_("SLATEC", "SOS", ch__5, &c__5, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)116);
	*iflag = 9;
    }

    if (*liw < *neq + 3) {
	s_wsfi(&io___14);
	do_fio(&c__1, (char *)&(*liw), (ftnlen)sizeof(integer));
	e_wsfi();
/* Writing concatenation */
	i__1[0] = 92, a__1[0] = "DIMENSION OF THE IW ARRAY MUST BE AT LEAST "
		"  3 + NEQ.  YOU HAVE CALLED THE CODE WITH  LIW = ";
	i__1[1] = 8, a__1[1] = xern1;
	s_cat(ch__6, a__1, i__1, &c__2, (ftnlen)100);
	xermsg_("SLATEC", "SOS", ch__6, &c__6, &c__1, (ftnlen)6, (ftnlen)3, (
		ftnlen)100);
	*iflag = 9;
    }

    if (*iflag != 9) {
	ncjs = 6;
	nsrrc = 4;
	nsri = 5;

	k1 = nc + 2;
	k2 = k1 + *neq;
	k3 = k2 + *neq;
	k4 = k3 + *neq;
	k5 = k4 + *neq;
	k6 = k5 + *neq;

	soseqs_((U_fp)fnc, neq, &x[1], rtolx, atolx, tolf, iflag, &mxit, &
		ncjs, &nsrrc, &nsri, &iprint, &rw[1], &rw[2], &nc, &rw[k1], &
		rw[k2], &rw[k3], &rw[k4], &rw[k5], &rw[k6], &iw[4]);

	iw[3] = mxit;
    }
    return 0;
} /* sos_ */

