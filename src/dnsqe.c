/* dnsqe.f -- translated by f2c (version 12.02.01).
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

static integer c__2 = 2;
static integer c__1 = 1;

/* DECK DNSQE */
/* Subroutine */ int dnsqe_(U_fp fcn, U_fp jac, integer *iopt, integer *n, 
	doublereal *x, doublereal *fvec, doublereal *tol, integer *nprint, 
	integer *info, doublereal *wa, integer *lwa)
{
    /* Initialized data */

    static doublereal factor = 100.;
    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j, ml, lr, mu, mode, nfev, njev;
    extern /* Subroutine */ int dnsq_(U_fp, U_fp, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal xtol;
    static integer index;
    static doublereal epsfcn;
    static integer maxfev;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DNSQE */
/* ***PURPOSE  An easy-to-use code to find a zero of a system of N */
/*            nonlinear functions in N variables by a modification of */
/*            the Powell hybrid method. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F2A */
/* ***TYPE      DOUBLE PRECISION (SNSQE-S, DNSQE-D) */
/* ***KEYWORDS  EASY-TO-USE, NONLINEAR SQUARE SYSTEM, */
/*             POWELL HYBRID METHOD, ZEROS */
/* ***AUTHOR  Hiebert, K. L. (SNLA) */
/* ***DESCRIPTION */

/* 1. Purpose. */

/*       The purpose of DNSQE is to find a zero of a system of N */
/*       nonlinear functions in N variables by a modification of the */
/*       Powell hybrid method.  This is done by using the more general */
/*       nonlinear equation solver DNSQ.  The user must provide a */
/*       subroutine which calculates the functions.  The user has the */
/*       option of either to provide a subroutine which calculates the */
/*       Jacobian or to let the code calculate it by a forward-difference */
/*       approximation.  This code is the combination of the MINPACK */
/*       codes (Argonne) HYBRD1 and HYBRJ1. */

/* 2. Subroutine and Type Statements. */

/*       SUBROUTINE DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO, */
/*      *                  WA,LWA) */
/*       INTEGER IOPT,N,NPRINT,INFO,LWA */
/*       DOUBLE PRECISION TOL */
/*       DOUBLE PRECISION X(N),FVEC(N),WA(LWA) */
/*       EXTERNAL FCN,JAC */

/* 3. Parameters. */

/*       Parameters designated as input parameters must be specified on */
/*       entry to DNSQE and are not changed on exit, while parameters */
/*       designated as output parameters need not be specified on entry */
/*       and are set to appropriate values on exit from DNSQE. */

/*       FCN is the name of the user-supplied subroutine which calculates */
/*         the functions.  FCN must be declared in an external statement */
/*         in the user calling program, and should be written as follows. */

/*         SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*         INTEGER N,IFLAG */
/*         DOUBLE PRECISION X(N),FVEC(N) */
/*         ---------- */
/*         Calculate the functions at X and */
/*         return this vector in FVEC. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by FCN unless the */
/*         user wants to terminate execution of DNSQE.  In this case set */
/*         IFLAG to a negative integer. */

/*       JAC is the name of the user-supplied subroutine which calculates */
/*         the Jacobian.  If IOPT=1, then JAC must be declared in an */
/*         external statement in the user calling program, and should be */
/*         written as follows. */

/*         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG) */
/*         INTEGER N,LDFJAC,IFLAG */
/*         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N) */
/*         ---------- */
/*         Calculate the Jacobian at X and return this */
/*         matrix in FJAC.  FVEC contains the function */
/*         values at X and should not be altered. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by JAC unless the */
/*         user wants to terminate execution of DNSQE. In this case set */
/*         IFLAG to a negative integer. */

/*         If IOPT=2, JAC can be ignored (treat it as a dummy argument). */

/*       IOPT is an input variable which specifies how the Jacobian will */
/*         be calculated.  If IOPT=1, then the user must supply the */
/*         Jacobian through the subroutine JAC.  If IOPT=2, then the */
/*         code will approximate the Jacobian by forward-differencing. */

/*       N is a positive integer input variable set to the number of */
/*         functions and variables. */

/*       X is an array of length N.  On input X must contain an initial */
/*         estimate of the solution vector.  On output X contains the */
/*         final estimate of the solution vector. */

/*       FVEC is an output array of length N which contains the functions */
/*         evaluated at the output X. */

/*       TOL is a nonnegative input variable.  Termination occurs when */
/*         the algorithm estimates that the relative error between X and */
/*         the solution is at most TOL.  Section 4 contains more details */
/*         about TOL. */

/*       NPRINT is an integer input variable that enables controlled */
/*         printing of iterates if it is positive.  In this case, FCN is */
/*         called with IFLAG = 0 at the beginning of the first iteration */
/*         and every NPRINT iterations thereafter and immediately prior */
/*         to return, with X and FVEC available for printing. Appropriate */
/*         print statements must be added to FCN(see example).  If NPRINT */
/*         is not positive, no special calls of FCN with IFLAG = 0 are */
/*         made. */

/*       INFO is an integer output variable.  If the user has terminated */
/*         execution, INFO is set to the (negative) value of IFLAG.  See */
/*         description of FCN and JAC. Otherwise, INFO is set as follows. */

/*         INFO = 0  Improper input parameters. */

/*         INFO = 1  Algorithm estimates that the relative error between */
/*                   X and the solution is at most TOL. */

/*         INFO = 2  Number of calls to FCN has reached or exceeded */
/*                   100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2. */

/*         INFO = 3  TOL is too small.  No further improvement in the */
/*                   approximate solution X is possible. */

/*         INFO = 4  Iteration is not making good progress. */

/*         Sections 4 and 5 contain more details about INFO. */

/*       WA is a work array of length LWA. */

/*       LWA is a positive integer input variable not less than */
/*         (3*N**2+13*N))/2. */

/* 4. Successful Completion. */

/*       The accuracy of DNSQE is controlled by the convergence parameter */
/*       TOL.  This parameter is used in a test which makes a comparison */
/*       between the approximation X and a solution XSOL.  DNSQE */
/*       terminates when the test is satisfied.  If TOL is less than the */
/*       machine precision (as defined by the  function D1MACH(4)), then */
/*       DNSQE only attempts to satisfy the test defined by the machine */
/*       precision.  Further progress is not usually possible.  Unless */
/*       high precision solutions are required, the recommended value */
/*       for TOL is the square root of the machine precision. */

/*       The test assumes that the functions are reasonably well behaved, */
/*       and, if the Jacobian is supplied by the user, that the functions */
/*       and the Jacobian are coded consistently. If these conditions are */
/*       not satisfied, then DNSQE may incorrectly indicate convergence. */
/*       The coding of the Jacobian can be checked by the subroutine */
/*       DCKDER.  If the Jacobian is coded correctly or IOPT=2, then */
/*       the validity of the answer can be checked, for example, by */
/*       rerunning DNSQE with a tighter tolerance. */

/*       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a */
/*         vector Z, then this test attempts to guarantee that */

/*               DENORM(X-XSOL) .LE. TOL*DENORM(XSOL). */

/*         If this condition is satisfied with TOL = 10**(-K), then the */
/*         larger components of X have K significant decimal digits and */
/*         INFO is set to 1.  There is a danger that the smaller */
/*         components of X may have large relative errors, but the fast */
/*         rate of convergence of DNSQE usually avoids this possibility. */

/* 5. Unsuccessful Completion. */

/*       Unsuccessful termination of DNSQE can be due to improper input */
/*       parameters, arithmetic interrupts, an excessive number of */
/*       function evaluations, errors in the functions, or lack of good */
/*       progress. */

/*       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1, or */
/*         IOPT .GT. 2, or N .LE. 0, or TOL .LT. 0.E0, or */
/*         LWA .LT. (3*N**2+13*N)/2. */

/*       Arithmetic Interrupts.  If these interrupts occur in the FCN */
/*         subroutine during an early stage of the computation, they may */
/*         be caused by an unacceptable choice of X by DNSQE.  In this */
/*         case, it may be possible to remedy the situation by not */
/*         evaluating the functions here, but instead setting the */
/*         components of FVEC to numbers that exceed those in the initial */
/*         FVEC. */

/*       Excessive Number of Function Evaluations.  If the number of */
/*         calls to FCN reaches 100*(N+1) for IOPT=1 or 200*(N+1) for */
/*         IOPT=2, then this indicates that the routine is converging */
/*         very slowly as measured by the progress of FVEC, and INFO is */
/*         set to 2.  This situation should be unusual because, as */
/*         indicated below, lack of good progress is usually diagnosed */
/*         earlier by DNSQE, causing termination with INFO = 4. */

/*       Errors In the Functions.  When IOPT=2, the choice of step length */
/*         in the forward-difference approximation to the Jacobian */
/*         assumes that the relative errors in the functions are of the */
/*         order of the machine precision.  If this is not the case, */
/*         DNSQE may fail (usually with INFO = 4).  The user should */
/*         then either use DNSQ and set the step length or use IOPT=1 */
/*         and supply the Jacobian. */

/*       Lack of Good Progress.  DNSQE searches for a zero of the system */
/*         by minimizing the sum of the squares of the functions.  In so */
/*         doing, it can become trapped in a region where the minimum */
/*         does not correspond to a zero of the system and, in this */
/*         situation, the iteration eventually fails to make good */
/*         progress.  In particular, this will happen if the system does */
/*         not have a zero.  If the system has a zero, rerunning DNSQE */
/*         from a different starting point may be helpful. */

/* 6. Characteristics of The Algorithm. */

/*       DNSQE is a modification of the Powell Hybrid method.  Two of */
/*       its main characteristics involve the choice of the correction as */
/*       a convex combination of the Newton and scaled gradient */
/*       directions, and the updating of the Jacobian by the rank-1 */
/*       method of Broyden.  The choice of the correction guarantees */
/*       (under reasonable conditions) global convergence for starting */
/*       points far from the solution and a fast rate of convergence. */
/*       The Jacobian is calculated at the starting point by either the */
/*       user-supplied subroutine or a forward-difference approximation, */
/*       but it is not recalculated until the rank-1 method fails to */
/*       produce satisfactory progress. */

/*       Timing.  The time required by DNSQE to solve a given problem */
/*         depends on N, the behavior of the functions, the accuracy */
/*         requested, and the starting point.  The number of arithmetic */
/*         operations needed by DNSQE is about 11.5*(N**2) to process */
/*         each evaluation of the functions (call to FCN) and 1.3*(N**3) */
/*         to process each evaluation of the Jacobian (call to JAC, */
/*         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly, */
/*         the timing of DNSQE will be strongly influenced by the time */
/*         spent in FCN and JAC. */

/*       Storage.  DNSQE requires (3*N**2 + 17*N)/2 single precision */
/*         storage locations, in addition to the storage required by the */
/*         program.  There are no internally declared storage arrays. */

/* *Long Description: */

/* 7. Example. */

/*       The problem is to determine the values of X(1), X(2), ..., X(9), */
/*       which solve the system of tridiagonal equations */

/*       (3-2*X(1))*X(1)           -2*X(2)                   = -1 */
/*               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8 */
/*                                   -X(8) + (3-2*X(9))*X(9) = -1 */

/*       ********** */

/*       PROGRAM TEST */
/* C */
/* C     DRIVER FOR DNSQE EXAMPLE. */
/* C */
/*       INTEGER J,N,IOPT,NPRINT,INFO,LWA,NWRITE */
/*       DOUBLE PRECISION TOL,FNORM */
/*       DOUBLE PRECISION X(9),FVEC(9),WA(180) */
/*       DOUBLE PRECISION DENORM,D1MACH */
/*       EXTERNAL FCN */
/*       DATA NWRITE /6/ */
/* C */
/*       IOPT = 2 */
/*       N = 9 */
/* C */
/* C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION. */
/* C */
/*       DO 10 J = 1, 9 */
/*          X(J) = -1.E0 */
/*    10    CONTINUE */

/*       LWA = 180 */
/*       NPRINT = 0 */
/* C */
/* C     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION. */
/* C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED, */
/* C     THIS IS THE RECOMMENDED SETTING. */
/* C */
/*       TOL = SQRT(D1MACH(4)) */
/* C */
/*       CALL DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA) */
/*       FNORM = DENORM(N,FVEC) */
/*       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N) */
/*       STOP */
/*  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 // */
/*      *        5X,' EXIT PARAMETER',16X,I10 // */
/*      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7)) */
/*       END */
/*       SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*       INTEGER N,IFLAG */
/*       DOUBLE PRECISION X(N),FVEC(N) */
/*       INTEGER K */
/*       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO */
/*       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/ */
/* C */
/*       DO 10 K = 1, N */
/*          TEMP = (THREE - TWO*X(K))*X(K) */
/*          TEMP1 = ZERO */
/*          IF (K .NE. 1) TEMP1 = X(K-1) */
/*          TEMP2 = ZERO */
/*          IF (K .NE. N) TEMP2 = X(K+1) */
/*          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE */
/*    10    CONTINUE */
/*       RETURN */
/*       END */

/*       RESULTS OBTAINED WITH DIFFERENT COMPILERS OR MACHINES */
/*       MAY BE SLIGHTLY DIFFERENT. */

/*       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07 */

/*       EXIT PARAMETER                         1 */

/*       FINAL APPROXIMATE SOLUTION */

/*       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00 */
/*       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00 */
/*       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00 */

/* ***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa- */
/*                 tions. In Numerical Methods for Nonlinear Algebraic */
/*                 Equations, P. Rabinowitz, Editor.  Gordon and Breach, */
/*                 1988. */
/* ***ROUTINES CALLED  DNSQ, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DNSQE */
    /* Parameter adjustments */
    --wa;
    --fvec;
    --x;

    /* Function Body */
/*     BEGIN BLOCK PERMITTING ...EXITS TO 20 */
/* ***FIRST EXECUTABLE STATEMENT  DNSQE */
    *info = 0;

/*        CHECK THE INPUT PARAMETERS FOR ERRORS. */

/*     ...EXIT */
/* Computing 2nd power */
    i__1 = *n;
    if (*iopt < 1 || *iopt > 2 || *n <= 0 || *tol < zero || *lwa < (i__1 * 
	    i__1 * 3 + *n * 13) / 2) {
	goto L20;
    }

/*        CALL DNSQ. */

    maxfev = (*n + 1) * 100;
    if (*iopt == 2) {
	maxfev <<= 1;
    }
    xtol = *tol;
    ml = *n - 1;
    mu = *n - 1;
    epsfcn = zero;
    mode = 2;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa[j] = one;
/* L10: */
    }
    lr = *n * (*n + 1) / 2;
    index = *n * 6 + lr;
    dnsq_((U_fp)fcn, (U_fp)jac, iopt, n, &x[1], &fvec[1], &wa[index + 1], n, &
	    xtol, &maxfev, &ml, &mu, &epsfcn, &wa[1], &mode, &factor, nprint, 
	    info, &nfev, &njev, &wa[*n * 6 + 1], &lr, &wa[*n + 1], &wa[(*n << 
	    1) + 1], &wa[*n * 3 + 1], &wa[(*n << 2) + 1], &wa[*n * 5 + 1]);
    if (*info == 5) {
	*info = 4;
    }
L20:
    if (*info == 0) {
	xermsg_("SLATEC", "DNSQE", "INVALID INPUT PARAMETER.", &c__2, &c__1, (
		ftnlen)6, (ftnlen)5, (ftnlen)24);
    }
    return 0;

/*     LAST CARD OF SUBROUTINE DNSQE. */

} /* dnsqe_ */

