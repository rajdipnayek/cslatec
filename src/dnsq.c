/* dnsq.f -- translated by f2c (version 12.02.01).
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

static integer c__4 = 4;
static logical c_false = FALSE_;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__9 = 9;
static integer c__3 = 3;

/* DECK DNSQ */
/* Subroutine */ int dnsq_(S_fp fcn, S_fp jac, integer *iopt, integer *n, 
	doublereal *x, doublereal *fvec, doublereal *fjac, integer *ldfjac, 
	doublereal *xtol, integer *maxfev, integer *ml, integer *mu, 
	doublereal *epsfcn, doublereal *diag, integer *mode, doublereal *
	factor, integer *nprint, integer *info, integer *nfev, integer *njev, 
	doublereal *r__, integer *lr, doublereal *qtf, doublereal *wa1, 
	doublereal *wa2, doublereal *wa3, doublereal *wa4)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal p1 = .1;
    static doublereal p5 = .5;
    static doublereal p001 = .001;
    static doublereal p0001 = 1e-4;
    static doublereal zero = 0.;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, l, jm1, iwa[1];
    static doublereal sum;
    static logical sing;
    static integer iter;
    static doublereal temp;
    static integer iflag;
    static doublereal delta;
    static logical jeval;
    static integer ncsuc;
    static doublereal ratio, fnorm, pnorm;
    extern /* Subroutine */ int dfdjc1_(S_fp, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    extern doublereal d1mach_(integer *);
    static doublereal xnorm;
    extern /* Subroutine */ int d1updt_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, logical *);
    static doublereal fnorm1;
    extern /* Subroutine */ int d1mpyq_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);
    static integer nslow1, nslow2, ncfail;
    extern /* Subroutine */ int dqrfac_(integer *, integer *, doublereal *, 
	    integer *, logical *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), ddoglg_(integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal actred, epsmch, prered;
    extern doublereal denorm_(integer *, doublereal *);
    extern /* Subroutine */ int dqform_(integer *, integer *, doublereal *, 
	    integer *, doublereal *), xermsg_(char *, char *, char *, integer 
	    *, integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  DNSQ */
/* ***PURPOSE  Find a zero of a system of a N nonlinear functions in N */
/*            variables by a modification of the Powell hybrid method. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F2A */
/* ***TYPE      DOUBLE PRECISION (SNSQ-S, DNSQ-D) */
/* ***KEYWORDS  NONLINEAR SQUARE SYSTEM, POWELL HYBRID METHOD, ZEROS */
/* ***AUTHOR  Hiebert, K. L. (SNLA) */
/* ***DESCRIPTION */

/* 1. Purpose. */

/*       The purpose of DNSQ is to find a zero of a system of N nonlinear */
/*       functions in N variables by a modification of the Powell */
/*       hybrid method.  The user must provide a subroutine which */
/*       calculates the functions.  The user has the option of either to */
/*       provide a subroutine which calculates the Jacobian or to let the */
/*       code calculate it by a forward-difference approximation. */
/*       This code is the combination of the MINPACK codes (Argonne) */
/*       HYBRD and HYBRDJ. */

/* 2. Subroutine and Type Statements. */

/*       SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV, */
/*      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV, */
/*      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4) */
/*       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR */
/*       DOUBLE PRECISION XTOL,EPSFCN,FACTOR */
/*       DOUBLE PRECISION */
/*       X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N), */
/*      *     WA1(N),WA2(N),WA3(N),WA4(N) */
/*       EXTERNAL FCN,JAC */

/* 3. Parameters. */

/*       Parameters designated as input parameters must be specified on */
/*       entry to DNSQ and are not changed on exit, while parameters */
/*       designated as output parameters need not be specified on entry */
/*       and are set to appropriate values on exit from DNSQ. */

/*       FCN is the name of the user-supplied subroutine which calculates */
/*         the functions.  FCN must be declared in an EXTERNAL statement */
/*         in the user calling program, and should be written as follows. */

/*         SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*         INTEGER N,IFLAG */
/*         DOUBLE PRECISION X(N),FVEC(N) */
/*         ---------- */
/*         CALCULATE THE FUNCTIONS AT X AND */
/*         RETURN THIS VECTOR IN FVEC. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by FCN unless the */
/*         user wants to terminate execution of DNSQ.  In this case set */
/*         IFLAG to a negative integer. */

/*       JAC is the name of the user-supplied subroutine which calculates */
/*         the Jacobian.  If IOPT=1, then JAC must be declared in an */
/*         EXTERNAL statement in the user calling program, and should be */
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
/*         user wants to terminate execution of DNSQ.  In this case set */
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

/*       FJAC is an output N by N array which contains the orthogonal */
/*         matrix Q produced by the QR factorization of the final */
/*         approximate Jacobian. */

/*       LDFJAC is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array FJAC. */

/*       XTOL is a nonnegative input variable.  Termination occurs when */
/*         the relative error between two consecutive iterates is at most */
/*         XTOL.  Therefore, XTOL measures the relative error desired in */
/*         the approximate solution.  Section 4 contains more details */
/*         about XTOL. */

/*       MAXFEV is a positive integer input variable.  Termination occurs */
/*         when the number of calls to FCN is at least MAXFEV by the end */
/*         of an iteration. */

/*       ML is a nonnegative integer input variable which specifies the */
/*         number of subdiagonals within the band of the Jacobian matrix. */
/*         If the Jacobian is not banded or IOPT=1, set ML to at */
/*         least N - 1. */

/*       MU is a nonnegative integer input variable which specifies the */
/*         number of superdiagonals within the band of the Jacobian */
/*         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at */
/*         least N - 1. */

/*       EPSFCN is an input variable used in determining a suitable step */
/*         for the forward-difference approximation.  This approximation */
/*         assumes that the relative errors in the functions are of the */
/*         order of EPSFCN.  If EPSFCN is less than the machine */
/*         precision, it is assumed that the relative errors in the */
/*         functions are of the order of the machine precision.  If */
/*         IOPT=1, then EPSFCN can be ignored (treat it as a dummy */
/*         argument). */

/*       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is */
/*         internally set.  If MODE = 2, DIAG must contain positive */
/*         entries that serve as implicit (multiplicative) scale factors */
/*         for the variables. */

/*       MODE is an integer input variable.  If MODE = 1, the variables */
/*         will be scaled internally.  If MODE = 2, the scaling is */
/*         specified by the input DIAG.  Other values of MODE are */
/*         equivalent to MODE = 1. */

/*       FACTOR is a positive input variable used in determining the */
/*         initial step bound.  This bound is set to the product of */
/*         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else to */
/*         FACTOR itself.  In most cases FACTOR should lie in the */
/*         interval (.1,100.).  100. is a generally recommended value. */

/*       NPRINT is an integer input variable that enables controlled */
/*         printing of iterates if it is positive.  In this case, FCN is */
/*         called with IFLAG = 0 at the beginning of the first iteration */
/*         and every NPRINT iterations thereafter and immediately prior */
/*         to return, with X and FVEC available for printing. appropriate */
/*         print statements must be added to FCN(see example).  If NPRINT */
/*         is not positive, no special calls of FCN with IFLAG = 0 are */
/*         made. */

/*       INFO is an integer output variable.  If the user has terminated */
/*         execution, INFO is set to the (negative) value of IFLAG.  See */
/*         description of FCN and JAC. Otherwise, INFO is set as follows. */

/*         INFO = 0  Improper input parameters. */

/*         INFO = 1  Relative error between two consecutive iterates is */
/*                   at most XTOL. */

/*         INFO = 2  Number of calls to FCN has reached or exceeded */
/*                   MAXFEV. */

/*         INFO = 3  XTOL is too small.  No further improvement in the */
/*                   approximate solution X is possible. */

/*         INFO = 4  Iteration is not making good progress, as measured */
/*                   by the improvement from the last five Jacobian */
/*                   evaluations. */

/*         INFO = 5  Iteration is not making good progress, as measured */
/*                   by the improvement from the last ten iterations. */

/*         Sections 4 and 5 contain more details about INFO. */

/*       NFEV is an integer output variable set to the number of calls to */
/*         FCN. */

/*       NJEV is an integer output variable set to the number of calls to */
/*         JAC. (If IOPT=2, then NJEV is set to zero.) */

/*       R is an output array of length LR which contains the upper */
/*         triangular matrix produced by the QR factorization of the */
/*         final approximate Jacobian, stored rowwise. */

/*       LR is a positive integer input variable not less than */
/*         (N*(N+1))/2. */

/*       QTF is an output array of length N which contains the vector */
/*         (Q transpose)*FVEC. */

/*       WA1, WA2, WA3, and WA4 are work arrays of length N. */


/* 4. Successful completion. */

/*       The accuracy of DNSQ is controlled by the convergence parameter */
/*       XTOL.  This parameter is used in a test which makes a comparison */
/*       between the approximation X and a solution XSOL.  DNSQ */
/*       terminates when the test is satisfied.  If the convergence */
/*       parameter is less than the machine precision (as defined by the */
/*       function D1MACH(4)), then DNSQ only attempts to satisfy the test */
/*       defined by the machine precision.  Further progress is not */
/*       usually possible. */

/*       The test assumes that the functions are reasonably well behaved, */
/*       and, if the Jacobian is supplied by the user, that the functions */
/*       and the Jacobian are coded consistently.  If these conditions */
/*       are not satisfied, then DNSQ may incorrectly indicate */
/*       convergence.  The coding of the Jacobian can be checked by the */
/*       subroutine DCKDER. If the Jacobian is coded correctly or IOPT=2, */
/*       then the validity of the answer can be checked, for example, by */
/*       rerunning DNSQ with a tighter tolerance. */

/*       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a */
/*         vector Z and D is the diagonal matrix whose entries are */
/*         defined by the array DIAG, then this test attempts to */
/*         guarantee that */

/*               DENORM(D*(X-XSOL)) .LE. XTOL*DENORM(D*XSOL). */

/*         If this condition is satisfied with XTOL = 10**(-K), then the */
/*         larger components of D*X have K significant decimal digits and */
/*         INFO is set to 1.  There is a danger that the smaller */
/*         components of D*X may have large relative errors, but the fast */
/*         rate of convergence of DNSQ usually avoids this possibility. */
/*         Unless high precision solutions are required, the recommended */
/*         value for XTOL is the square root of the machine precision. */


/* 5. Unsuccessful Completion. */

/*       Unsuccessful termination of DNSQ can be due to improper input */
/*       parameters, arithmetic interrupts, an excessive number of */
/*       function evaluations, or lack of good progress. */

/*       Improper Input Parameters.  INFO is set to 0 if IOPT .LT .1, */
/*         or IOPT .GT. 2, or N .LE. 0, or LDFJAC .LT. N, or */
/*         XTOL .LT. 0.E0, or MAXFEV .LE. 0, or ML .LT. 0, or MU .LT. 0, */
/*         or FACTOR .LE. 0.E0, or LR .LT. (N*(N+1))/2. */

/*       Arithmetic Interrupts.  If these interrupts occur in the FCN */
/*         subroutine during an early stage of the computation, they may */
/*         be caused by an unacceptable choice of X by DNSQ.  In this */
/*         case, it may be possible to remedy the situation by rerunning */
/*         DNSQ with a smaller value of FACTOR. */

/*       Excessive Number of Function Evaluations.  A reasonable value */
/*         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2. */
/*         If the number of calls to FCN reaches MAXFEV, then this */
/*         indicates that the routine is converging very slowly as */
/*         measured by the progress of FVEC, and INFO is set to 2. This */
/*         situation should be unusual because, as indicated below, lack */
/*         of good progress is usually diagnosed earlier by DNSQ, */
/*         causing termination with info = 4 or INFO = 5. */

/*       Lack of Good Progress.  DNSQ searches for a zero of the system */
/*         by minimizing the sum of the squares of the functions.  In so */
/*         doing, it can become trapped in a region where the minimum */
/*         does not correspond to a zero of the system and, in this */
/*         situation, the iteration eventually fails to make good */
/*         progress.  In particular, this will happen if the system does */
/*         not have a zero.  If the system has a zero, rerunning DNSQ */
/*         from a different starting point may be helpful. */


/* 6. Characteristics of The Algorithm. */

/*       DNSQ is a modification of the Powell Hybrid method.  Two of its */
/*       main characteristics involve the choice of the correction as a */
/*       convex combination of the Newton and scaled gradient directions, */
/*       and the updating of the Jacobian by the rank-1 method of */
/*       Broyden.  The choice of the correction guarantees (under */
/*       reasonable conditions) global convergence for starting points */
/*       far from the solution and a fast rate of convergence.  The */
/*       Jacobian is calculated at the starting point by either the */
/*       user-supplied subroutine or a forward-difference approximation, */
/*       but it is not recalculated until the rank-1 method fails to */
/*       produce satisfactory progress. */

/*       Timing.  The time required by DNSQ to solve a given problem */
/*         depends on N, the behavior of the functions, the accuracy */
/*         requested, and the starting point.  The number of arithmetic */
/*         operations needed by DNSQ is about 11.5*(N**2) to process */
/*         each evaluation of the functions (call to FCN) and 1.3*(N**3) */
/*         to process each evaluation of the Jacobian (call to JAC, */
/*         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly, */
/*         the timing of DNSQ will be strongly influenced by the time */
/*         spent in FCN and JAC. */

/*       Storage.  DNSQ requires (3*N**2 + 17*N)/2 single precision */
/*         storage locations, in addition to the storage required by the */
/*         program.  There are no internally declared storage arrays. */

/* *Long Description: */

/* 7. Example. */

/*       The problem is to determine the values of X(1), X(2), ..., X(9), */
/*       which solve the system of tridiagonal equations */

/*       (3-2*X(1))*X(1)           -2*X(2)                   = -1 */
/*               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8 */
/*                                   -X(8) + (3-2*X(9))*X(9) = -1 */
/* C     ********** */

/*       PROGRAM TEST */
/* C */
/* C     Driver for DNSQ example. */
/* C */
/*       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR, */
/*      *        NWRITE */
/*       DOUBLE PRECISION XTOL,EPSFCN,FACTOR,FNORM */
/*       DOUBLE PRECISION X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9), */
/*      *     WA1(9),WA2(9),WA3(9),WA4(9) */
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
/* C */
/*       LDFJAC = 9 */
/*       LR = 45 */
/* C */
/* C     SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION. */
/* C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED, */
/* C     THIS IS THE RECOMMENDED SETTING. */
/* C */
/*       XTOL = SQRT(D1MACH(4)) */
/* C */
/*       MAXFEV = 2000 */
/*       ML = 1 */
/*       MU = 1 */
/*       EPSFCN = 0.E0 */
/*       MODE = 2 */
/*       DO 20 J = 1, 9 */
/*          DIAG(J) = 1.E0 */
/*    20    CONTINUE */
/*       FACTOR = 1.E2 */
/*       NPRINT = 0 */
/* C */
/*       CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU, */
/*      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV, */
/*      *           R,LR,QTF,WA1,WA2,WA3,WA4) */
/*       FNORM = DENORM(N,FVEC) */
/*       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N) */
/*       STOP */
/*  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 // */
/*      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 // */
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
/*       IF (IFLAG .NE. 0) GO TO 5 */
/* C */
/* C     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE. */
/* C */
/*       RETURN */
/*     5 CONTINUE */
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

/*       Results obtained with different compilers or machines */
/*       may be slightly different. */

/*       Final L2 norm of the residuals  0.1192636E-07 */

/*       Number of function evaluations        14 */

/*       Exit parameter                         1 */

/*       Final approximate solution */

/*       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00 */
/*       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00 */
/*       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00 */

/* ***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa- */
/*                 tions. In Numerical Methods for Nonlinear Algebraic */
/*                 Equations, P. Rabinowitz, Editor.  Gordon and Breach, */
/*                 1988. */
/* ***ROUTINES CALLED  D1MACH, D1MPYQ, D1UPDT, DDOGLG, DENORM, DFDJC1, */
/*                    DQFORM, DQRFAC, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  DNSQ */
    /* Parameter adjustments */
    --x;
    --fvec;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --diag;
    --r__;
    --qtf;
    --wa1;
    --wa2;
    --wa3;
    --wa4;

    /* Function Body */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 320 */
/* ***FIRST EXECUTABLE STATEMENT  DNSQ */
    epsmch = d1mach_(&c__4);

    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

/*        CHECK THE INPUT PARAMETERS FOR ERRORS. */

/*     ...EXIT */
    if (*iopt < 1 || *iopt > 2 || *n <= 0 || *xtol < zero || *maxfev <= 0 || *
	    ml < 0 || *mu < 0 || *factor <= zero || *ldfjac < *n || *lr < *n *
	     (*n + 1) / 2) {
	goto L320;
    }
    if (*mode != 2) {
	goto L20;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/*     .........EXIT */
	if (diag[j] <= zero) {
	    goto L320;
	}
/* L10: */
    }
L20:

/*        EVALUATE THE FUNCTION AT THE STARTING POINT */
/*        AND CALCULATE ITS NORM. */

    iflag = 1;
    (*fcn)(n, &x[1], &fvec[1], &iflag);
    *nfev = 1;
/*     ...EXIT */
    if (iflag < 0) {
	goto L320;
    }
    fnorm = denorm_(n, &fvec[1]);

/*        INITIALIZE ITERATION COUNTER AND MONITORS. */

    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;

/*        BEGINNING OF THE OUTER LOOP. */

L30:
/*           BEGIN BLOCK PERMITTING ...EXITS TO 90 */
    jeval = TRUE_;

/*              CALCULATE THE JACOBIAN MATRIX. */

    if (*iopt == 2) {
	goto L40;
    }

/*                 USER SUPPLIES JACOBIAN */

    (*jac)(n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    ++(*njev);
    goto L50;
L40:

/*                 CODE APPROXIMATES THE JACOBIAN */

    iflag = 2;
    dfdjc1_((S_fp)fcn, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag,
	     ml, mu, epsfcn, &wa1[1], &wa2[1]);
/* Computing MIN */
    i__1 = *ml + *mu + 1;
    *nfev += min(i__1,*n);
L50:

/*     .........EXIT */
    if (iflag < 0) {
	goto L320;
    }

/*              COMPUTE THE QR FACTORIZATION OF THE JACOBIAN. */

    dqrfac_(n, n, &fjac[fjac_offset], ldfjac, &c_false, iwa, &c__1, &wa1[1], &
	    wa2[1], &wa3[1]);

/*              ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING */
/*              TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN. */

/*           ...EXIT */
    if (iter != 1) {
	goto L90;
    }
    if (*mode == 2) {
	goto L70;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	diag[j] = wa2[j];
	if (wa2[j] == zero) {
	    diag[j] = one;
	}
/* L60: */
    }
L70:

/*              ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED */
/*              X AND INITIALIZE THE STEP BOUND DELTA. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa3[j] = diag[j] * x[j];
/* L80: */
    }
    xnorm = denorm_(n, &wa3[1]);
    delta = *factor * xnorm;
    if (delta == zero) {
	delta = *factor;
    }
L90:

/*           FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	qtf[i__] = fvec[i__];
/* L100: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (fjac[j + j * fjac_dim1] == zero) {
	    goto L130;
	}
	sum = zero;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum += fjac[i__ + j * fjac_dim1] * qtf[i__];
/* L110: */
	}
	temp = -sum / fjac[j + j * fjac_dim1];
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    qtf[i__] += fjac[i__ + j * fjac_dim1] * temp;
/* L120: */
	}
L130:
/* L140: */
	;
    }

/*           COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R. */

    sing = FALSE_;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L160;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[l] = fjac[i__ + j * fjac_dim1];
	    l = l + *n - i__;
/* L150: */
	}
L160:
	r__[l] = wa1[j];
	if (wa1[j] == zero) {
	    sing = TRUE_;
	}
/* L170: */
    }

/*           ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC. */

    dqform_(n, n, &fjac[fjac_offset], ldfjac, &wa1[1]);

/*           RESCALE IF NECESSARY. */

    if (*mode == 2) {
	goto L190;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = diag[j], d__2 = wa2[j];
	diag[j] = max(d__1,d__2);
/* L180: */
    }
L190:

/*           BEGINNING OF THE INNER LOOP. */

L200:

/*              IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES. */

    if (*nprint <= 0) {
	goto L210;
    }
    iflag = 0;
    if ((iter - 1) % *nprint == 0) {
	(*fcn)(n, &x[1], &fvec[1], &iflag);
    }
/*     ............EXIT */
    if (iflag < 0) {
	goto L320;
    }
L210:

/*              DETERMINE THE DIRECTION P. */

    ddoglg_(n, &r__[1], lr, &diag[1], &qtf[1], &delta, &wa1[1], &wa2[1], &wa3[
	    1]);

/*              STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = -wa1[j];
	wa2[j] = x[j] + wa1[j];
	wa3[j] = diag[j] * wa1[j];
/* L220: */
    }
    pnorm = denorm_(n, &wa3[1]);

/*              ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND. */

    if (iter == 1) {
	delta = min(delta,pnorm);
    }

/*              EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM. */

    iflag = 1;
    (*fcn)(n, &wa2[1], &wa4[1], &iflag);
    ++(*nfev);
/*     .........EXIT */
    if (iflag < 0) {
	goto L320;
    }
    fnorm1 = denorm_(n, &wa4[1]);

/*              COMPUTE THE SCALED ACTUAL REDUCTION. */

    actred = -one;
    if (fnorm1 < fnorm) {
/* Computing 2nd power */
	d__1 = fnorm1 / fnorm;
	actred = one - d__1 * d__1;
    }

/*              COMPUTE THE SCALED PREDICTED REDUCTION. */

    l = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = zero;
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    sum += r__[l] * wa1[j];
	    ++l;
/* L230: */
	}
	wa3[i__] = qtf[i__] + sum;
/* L240: */
    }
    temp = denorm_(n, &wa3[1]);
    prered = zero;
    if (temp < fnorm) {
/* Computing 2nd power */
	d__1 = temp / fnorm;
	prered = one - d__1 * d__1;
    }

/*              COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED */
/*              REDUCTION. */

    ratio = zero;
    if (prered > zero) {
	ratio = actred / prered;
    }

/*              UPDATE THE STEP BOUND. */

    if (ratio >= p1) {
	goto L250;
    }
    ncsuc = 0;
    ++ncfail;
    delta = p5 * delta;
    goto L260;
L250:
    ncfail = 0;
    ++ncsuc;
    if (ratio >= p5 || ncsuc > 1) {
/* Computing MAX */
	d__1 = delta, d__2 = pnorm / p5;
	delta = max(d__1,d__2);
    }
    if ((d__1 = ratio - one, abs(d__1)) <= p1) {
	delta = pnorm / p5;
    }
L260:

/*              TEST FOR SUCCESSFUL ITERATION. */

    if (ratio < p0001) {
	goto L280;
    }

/*                 SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = wa2[j];
	wa2[j] = diag[j] * x[j];
	fvec[j] = wa4[j];
/* L270: */
    }
    xnorm = denorm_(n, &wa2[1]);
    fnorm = fnorm1;
    ++iter;
L280:

/*              DETERMINE THE PROGRESS OF THE ITERATION. */

    ++nslow1;
    if (actred >= p001) {
	nslow1 = 0;
    }
    if (jeval) {
	++nslow2;
    }
    if (actred >= p1) {
	nslow2 = 0;
    }

/*              TEST FOR CONVERGENCE. */

    if (delta <= *xtol * xnorm || fnorm == zero) {
	*info = 1;
    }
/*     .........EXIT */
    if (*info != 0) {
	goto L320;
    }

/*              TESTS FOR TERMINATION AND STRINGENT TOLERANCES. */

    if (*nfev >= *maxfev) {
	*info = 2;
    }
/* Computing MAX */
    d__1 = p1 * delta;
    if (p1 * max(d__1,pnorm) <= epsmch * xnorm) {
	*info = 3;
    }
    if (nslow2 == 5) {
	*info = 4;
    }
    if (nslow1 == 10) {
	*info = 5;
    }
/*     .........EXIT */
    if (*info != 0) {
	goto L320;
    }

/*              CRITERION FOR RECALCULATING JACOBIAN */

/*           ...EXIT */
    if (ncfail == 2) {
	goto L310;
    }

/*              CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN */
/*              AND UPDATE QTF IF NECESSARY. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum = zero;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
/* L290: */
	}
	wa2[j] = (sum - wa3[j]) / pnorm;
	wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
	if (ratio >= p0001) {
	    qtf[j] = sum;
	}
/* L300: */
    }

/*              COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN. */

    d1updt_(n, n, &r__[1], lr, &wa1[1], &wa2[1], &wa3[1], &sing);
    d1mpyq_(n, n, &fjac[fjac_offset], ldfjac, &wa2[1], &wa3[1]);
    d1mpyq_(&c__1, n, &qtf[1], &c__1, &wa2[1], &wa3[1]);

/*              END OF THE INNER LOOP. */

    jeval = FALSE_;
    goto L200;
L310:

/*           END OF THE OUTER LOOP. */

    goto L30;
L320:

/*     TERMINATION, EITHER NORMAL OR USER IMPOSED. */

    if (iflag < 0) {
	*info = iflag;
    }
    iflag = 0;
    if (*nprint > 0) {
	(*fcn)(n, &x[1], &fvec[1], &iflag);
    }
    if (*info < 0) {
	xermsg_("SLATEC", "DNSQ", "EXECUTION TERMINATED BECAUSE USER SET IFL"
		"AG NEGATIVE.", &c__1, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)53)
		;
    }
    if (*info == 0) {
	xermsg_("SLATEC", "DNSQ", "INVALID INPUT PARAMETER.", &c__2, &c__1, (
		ftnlen)6, (ftnlen)4, (ftnlen)24);
    }
    if (*info == 2) {
	xermsg_("SLATEC", "DNSQ", "TOO MANY FUNCTION EVALUATIONS.", &c__9, &
		c__1, (ftnlen)6, (ftnlen)4, (ftnlen)30);
    }
    if (*info == 3) {
	xermsg_("SLATEC", "DNSQ", "XTOL TOO SMALL. NO FURTHER IMPROVEMENT PO"
		"SSIBLE.", &c__3, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)48);
    }
    if (*info > 4) {
	xermsg_("SLATEC", "DNSQ", "ITERATION NOT MAKING GOOD PROGRESS.", &
		c__1, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)35);
    }
    return 0;

/*     LAST CARD OF SUBROUTINE DNSQ. */

} /* dnsq_ */

