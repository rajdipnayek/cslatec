/* snsq.f -- translated by f2c (version 12.02.01).
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

/* DECK SNSQ */
/* Subroutine */ int snsq_(S_fp fcn, S_fp jac, integer *iopt, integer *n, 
	real *x, real *fvec, real *fjac, integer *ldfjac, real *xtol, integer 
	*maxfev, integer *ml, integer *mu, real *epsfcn, real *diag, integer *
	mode, real *factor, integer *nprint, integer *info, integer *nfev, 
	integer *njev, real *r__, integer *lr, real *qtf, real *wa1, real *
	wa2, real *wa3, real *wa4)
{
    /* Initialized data */

    static real one = 1.f;
    static real p1 = .1f;
    static real p5 = .5f;
    static real p001 = .001f;
    static real p0001 = 1e-4f;
    static real zero = 0.f;

    /* System generated locals */
    integer fjac_dim1, fjac_offset, i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    static integer i__, j, l, jm1, iwa[1];
    static real sum;
    static logical sing;
    static integer iter;
    static real temp;
    static integer iflag;
    static real delta;
    extern /* Subroutine */ int qrfac_(integer *, integer *, real *, integer *
	    , logical *, integer *, integer *, real *, real *, real *);
    static logical jeval;
    static integer ncsuc;
    static real ratio;
    extern doublereal enorm_(integer *, real *);
    static real fnorm;
    extern /* Subroutine */ int qform_(integer *, integer *, real *, integer *
	    , real *), fdjac1_(S_fp, integer *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, real *, real *, real *
	    );
    static real pnorm, xnorm;
    extern doublereal r1mach_(integer *);
    static real fnorm1;
    extern /* Subroutine */ int r1updt_(integer *, integer *, real *, integer 
	    *, real *, real *, real *, logical *);
    static integer nslow1, nslow2;
    extern /* Subroutine */ int r1mpyq_(integer *, integer *, real *, integer 
	    *, real *, real *);
    static integer ncfail;
    extern /* Subroutine */ int dogleg_(integer *, real *, integer *, real *, 
	    real *, real *, real *, real *, real *);
    static real actred, epsmch, prered;
    extern /* Subroutine */ int xermsg_(char *, char *, char *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);

/* ***BEGIN PROLOGUE  SNSQ */
/* ***PURPOSE  Find a zero of a system of a N nonlinear functions in N */
/*            variables by a modification of the Powell hybrid method. */
/* ***LIBRARY   SLATEC */
/* ***CATEGORY  F2A */
/* ***TYPE      SINGLE PRECISION (SNSQ-S, DNSQ-D) */
/* ***KEYWORDS  NONLINEAR SQUARE SYSTEM, POWELL HYBRID METHOD, ZEROS */
/* ***AUTHOR  Hiebert, K. L., (SNLA) */
/* ***DESCRIPTION */

/* 1. Purpose. */

/*       The purpose of SNSQ is to find a zero of a system of N non- */
/*       linear functions in N variables by a modification of the Powell */
/*       hybrid method.  The user must provide a subroutine which calcu- */
/*       lates the functions.  The user has the option of either to */
/*       provide a subroutine which calculates the Jacobian or to let the */
/*       code calculate it by a forward-difference approximation. */
/*       This code is the combination of the MINPACK codes (Argonne) */
/*       HYBRD and HYBRDJ. */


/* 2. Subroutine and Type Statements. */

/*       SUBROUTINE SNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV, */
/*      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV, */
/*      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4) */
/*       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR */
/*       REAL XTOL,EPSFCN,FACTOR */
/*       REAL X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N), */
/*      *     WA1(N),WA2(N),WA3(N),WA4(N) */
/*       EXTERNAL FCN,JAC */


/* 3. Parameters. */

/*       Parameters designated as input parameters must be specified on */
/*       entry to SNSQ and are not changed on exit, while parameters */
/*       designated as output parameters need not be specified on entry */
/*       and are set to appropriate values on exit from SNSQ. */

/*       FCN is the name of the user-supplied subroutine which calculates */
/*         the functions.  FCN must be declared in an EXTERNAL statement */
/*         in the user calling program, and should be written as follows. */

/*         SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*         INTEGER N,IFLAG */
/*         REAL X(N),FVEC(N) */
/*         ---------- */
/*         Calculate the functions at X and */
/*         return this vector in FVEC. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by FCN unless the */
/*         user wants to terminate execution of SNSQ.  In this case, set */
/*         IFLAG to a negative integer. */

/*       JAC is the name of the user-supplied subroutine which calculates */
/*         the Jacobian.  If IOPT=1, then JAC must be declared in an */
/*         EXTERNAL statement in the user calling program, and should be */
/*         written as follows. */

/*         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG) */
/*         INTEGER N,LDFJAC,IFLAG */
/*         REAL X(N),FVEC(N),FJAC(LDFJAC,N) */
/*         ---------- */
/*         Calculate the Jacobian at X and return this */
/*         matrix in FJAC.  FVEC contains the function */
/*         values at X and should not be altered. */
/*         ---------- */
/*         RETURN */
/*         END */

/*         The value of IFLAG should not be changed by JAC unless the */
/*         user wants to terminate execution of SNSQ.  In this case, set */
/*         IFLAG to a negative integer. */

/*         If IOPT=2, JAC can be ignored (treat it as a dummy argument). */

/*       IOPT is an input variable which specifies how the Jacobian will */
/*         be calculated.  If IOPT=1, then the user must supply the */
/*         Jacobian through the subroutine JAC.  If IOPT=2, then the */
/*         code will approximate the Jacobian by forward-differencing. */

/*       N is a positive integer input variable set to the number of */
/*         functions and variables. */

/*       X is an array of length N.  On input, X must contain an initial */
/*         estimate of the solution vector.  On output, X contains the */
/*         final estimate of the solution vector. */

/*       FVEC is an output array of length N which contains the functions */
/*         evaluated at the output X. */

/*       FJAC is an output N by N array which contains the orthogonal */
/*         matrix Q produced by the QR factorization of the final approx- */
/*         imate Jacobian. */

/*       LDFJAC is a positive integer input variable not less than N */
/*         which specifies the leading dimension of the array FJAC. */

/*       XTOL is a non-negative input variable.  Termination occurs when */
/*         the relative error between two consecutive iterates is at most */
/*         XTOL.  Therefore, XTOL measures the relative error desired in */
/*         the approximate solution.  Section 4 contains more details */
/*         about XTOL. */

/*       MAXFEV is a positive integer input variable.  Termination occurs */
/*         when the number of calls to FCN is at least MAXFEV by the end */
/*         of an iteration. */

/*       ML is a non-negative integer input variable which specifies the */
/*         number of subdiagonals within the band of the Jacobian matrix. */
/*         If the Jacobian is not banded or IOPT=1, set ML to at */
/*         least N - 1. */

/*       MU is a non-negative integer input variable which specifies the */
/*         number of superdiagonals within the band of the Jacobian */
/*         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at */
/*         least N - 1. */

/*       EPSFCN is an input variable used in determining a suitable step */
/*         for the forward-difference approximation.  This approximation */
/*         assumes that the relative errors in the functions are of the */
/*         order of EPSFCN.  If EPSFCN is less than the machine preci- */
/*         sion, it is assumed that the relative errors in the functions */
/*         are of the order of the machine precision.  If IOPT=1, then */
/*         EPSFCN can be ignored (treat it as a dummy argument). */

/*       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is */
/*         internally set.  If MODE = 2, DIAG must contain positive */
/*         entries that serve as implicit (multiplicative) scale factors */
/*         for the variables. */

/*       MODE is an integer input variable.  If MODE = 1, the variables */
/*         will be scaled internally.  If MODE = 2, the scaling is speci- */
/*         fied by the input DIAG.  Other values of MODE are equivalent */
/*         to MODE = 1. */

/*       FACTOR is a positive input variable used in determining the ini- */
/*         tial step bound.  This bound is set to the product of FACTOR */
/*         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR */
/*         itself.  In most cases FACTOR should lie in the interval */
/*         (.1,100.).  100. is a generally recommended value. */

/*       NPRINT is an integer input variable that enables controlled */
/*         printing of iterates if it is positive.  In this case, FCN is */
/*         called with IFLAG = 0 at the beginning of the first iteration */
/*         and every NPRINT iteration thereafter and immediately prior */
/*         to return, with X and FVEC available for printing. Appropriate */
/*         print statements must be added to FCN(see example).  If NPRINT */
/*         is not positive, no special calls of FCN with IFLAG = 0 are */
/*         made. */

/*       INFO is an integer output variable.  If the user has terminated */
/*         execution, INFO is set to the (negative) value of IFLAG.  See */
/*         description of FCN and JAC. Otherwise, INFO is set as follows. */

/*         INFO = 0  improper input parameters. */

/*         INFO = 1  relative error between two consecutive iterates is */
/*                   at most XTOL. */

/*         INFO = 2  number of calls to FCN has reached or exceeded */
/*                   MAXFEV. */

/*         INFO = 3  XTOL is too small.  No further improvement in the */
/*                   approximate solution X is possible. */

/*         INFO = 4  iteration is not making good progress, as measured */
/*                   by the improvement from the last five Jacobian eval- */
/*                   uations. */

/*         INFO = 5  iteration is not making good progress, as measured */
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
/*         (Q TRANSPOSE)*FVEC. */

/*       WA1, WA2, WA3, and WA4 are work arrays of length N. */


/* 4. Successful Completion. */

/*       The accuracy of SNSQ is controlled by the convergence parameter */
/*       XTOL.  This parameter is used in a test which makes a comparison */
/*       between the approximation X and a solution XSOL.  SNSQ termi- */
/*       nates when the test is satisfied.  If the convergence parameter */
/*       is less than the machine precision (as defined by the function */
/*       R1MACH(4)), then SNSQ only attempts to satisfy the test */
/*       defined by the machine precision.  Further progress is not */
/*       usually possible. */

/*       The test assumes that the functions are reasonably well behaved, */
/*       and, if the Jacobian is supplied by the user, that the functions */
/*       and the Jacobian are coded consistently.  If these conditions */
/*       are not satisfied, then SNSQ may incorrectly indicate conver- */
/*       gence.  The coding of the Jacobian can be checked by the */
/*       subroutine CHKDER. If the Jacobian is coded correctly or IOPT=2, */
/*       then the validity of the answer can be checked, for example, by */
/*       rerunning SNSQ with a tighter tolerance. */

/*       Convergence Test.  If ENORM(Z) denotes the Euclidean norm of a */
/*         vector Z and D is the diagonal matrix whose entries are */
/*         defined by the array DIAG, then this test attempts to guaran- */
/*         tee that */

/*               ENORM(D*(X-XSOL)) .LE. XTOL*ENORM(D*XSOL). */

/*         If this condition is satisfied with XTOL = 10**(-K), then the */
/*         larger components of D*X have K significant decimal digits and */
/*         INFO is set to 1.  There is a danger that the smaller compo- */
/*         nents of D*X may have large relative errors, but the fast rate */
/*         of convergence of SNSQ usually avoids this possibility. */
/*         Unless high precision solutions are required, the recommended */
/*         value for XTOL is the square root of the machine precision. */


/* 5. Unsuccessful Completion. */

/*       Unsuccessful termination of SNSQ can be due to improper input */
/*       parameters, arithmetic interrupts, an excessive number of func- */
/*       tion evaluations, or lack of good progress. */

/*       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1, */
/*         or IOPT .GT. 2, or N .LE. 0, or LDFJAC .LT. N, or */
/*         XTOL .LT. 0.E0, or MAXFEV .LE. 0, or ML .LT. 0, or MU .LT. 0, */
/*         or FACTOR .LE. 0.E0, or LR .LT. (N*(N+1))/2. */

/*       Arithmetic Interrupts.  If these interrupts occur in the FCN */
/*         subroutine during an early stage of the computation, they may */
/*         be caused by an unacceptable choice of X by SNSQ.  In this */
/*         case, it may be possible to remedy the situation by rerunning */
/*         SNSQ with a smaller value of FACTOR. */

/*       Excessive Number of Function Evaluations.  A reasonable value */
/*         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2. */
/*         If the number of calls to FCN reaches MAXFEV, then this */
/*         indicates that the routine is converging very slowly as */
/*         measured by the progress of FVEC, and INFO is set to 2.  This */
/*         situation should be unusual because, as indicated below, lack */
/*         of good progress is usually diagnosed earlier by SNSQ, */
/*         causing termination with INFO = 4 or INFO = 5. */

/*       Lack of Good Progress.  SNSQ searches for a zero of the system */
/*         by minimizing the sum of the squares of the functions.  In so */
/*         doing, it can become trapped in a region where the minimum */
/*         does not correspond to a zero of the system and, in this situ- */
/*         ation, the iteration eventually fails to make good progress. */
/*         In particular, this will happen if the system does not have a */
/*         zero.  If the system has a zero, rerunning SNSQ from a dif- */
/*         ferent starting point may be helpful. */


/* 6. Characteristics of the Algorithm. */

/*       SNSQ is a modification of the Powell hybrid method.  Two of its */
/*       main characteristics involve the choice of the correction as a */
/*       convex combination of the Newton and scaled gradient directions, */
/*       and the updating of the Jacobian by the rank-1 method of Broy- */
/*       den.  The choice of the correction guarantees (under reasonable */
/*       conditions) global convergence for starting points far from the */
/*       solution and a fast rate of convergence.  The Jacobian is */
/*       calculated at the starting point by either the user-supplied */
/*       subroutine or a forward-difference approximation, but it is not */
/*       recalculated until the rank-1 method fails to produce satis- */
/*       factory progress. */

/*       Timing.  The time required by SNSQ to solve a given problem */
/*         depends on N, the behavior of the functions, the accuracy */
/*         requested, and the starting point.  The number of arithmetic */
/*         operations needed by SNSQ is about 11.5*(N**2) to process */
/*         each evaluation of the functions (call to FCN) and 1.3*(N**3) */
/*         to process each evaluation of the Jacobian (call to JAC, */
/*         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly, */
/*         the timing of SNSQ will be strongly influenced by the time */
/*         spent in FCN and JAC. */

/*       Storage.  SNSQ requires (3*N**2 + 17*N)/2 single precision */
/*         storage locations, in addition to the storage required by the */
/*         program.  There are no internally declared storage arrays. */


/* 7. Example. */

/*       The problem is to determine the values of X(1), X(2), ..., X(9), */
/*       which solve the system of tridiagonal equations */

/*       (3-2*X(1))*X(1)           -2*X(2)                   = -1 */
/*               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8 */
/*                                   -X(8) + (3-2*X(9))*X(9) = -1 */
/* C     ********** */

/*       PROGRAM TEST */
/* C */
/* C     Driver for SNSQ example. */
/* C */
/*       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR, */
/*      *        NWRITE */
/*       REAL XTOL,EPSFCN,FACTOR,FNORM */
/*       REAL X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9), */
/*      *     WA1(9),WA2(9),WA3(9),WA4(9) */
/*       REAL ENORM,R1MACH */
/*       EXTERNAL FCN */
/*       DATA NWRITE /6/ */
/* C */
/*       IOPT = 2 */
/*       N = 9 */
/* C */
/* C     The following starting values provide a rough solution. */
/* C */
/*       DO 10 J = 1, 9 */
/*          X(J) = -1.E0 */
/*    10    CONTINUE */
/* C */
/*       LDFJAC = 9 */
/*       LR = 45 */
/* C */
/* C     Set XTOL to the square root of the machine precision. */
/* C     Unless high precision solutions are required, */
/* C     this is the recommended setting. */
/* C */
/*       XTOL = SQRT(R1MACH(4)) */
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
/*       CALL SNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU, */
/*      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV, */
/*      *           R,LR,QTF,WA1,WA2,WA3,WA4) */
/*       FNORM = ENORM(N,FVEC) */
/*       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N) */
/*       STOP */
/*  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 // */
/*      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 // */
/*      *        5X,' EXIT PARAMETER',16X,I10 // */
/*      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7)) */
/*       END */
/*       SUBROUTINE FCN(N,X,FVEC,IFLAG) */
/*       INTEGER N,IFLAG */
/*       REAL X(N),FVEC(N) */
/*       INTEGER K */
/*       REAL ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO */
/*       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/ */
/* C */
/*       IF (IFLAG .NE. 0) GO TO 5 */
/* C */
/* C     Insert print statements here when NPRINT is positive. */
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

/*       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07 */

/*       NUMBER OF FUNCTION EVALUATIONS        14 */

/*       EXIT PARAMETER                         1 */

/*       FINAL APPROXIMATE SOLUTION */

/*       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00 */
/*       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00 */
/*       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00 */

/* ***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa- */
/*                 tions. In Numerical Methods for Nonlinear Algebraic */
/*                 Equations, P. Rabinowitz, Editor.  Gordon and Breach, */
/*                 1988. */
/* ***ROUTINES CALLED  DOGLEG, ENORM, FDJAC1, QFORM, QRFAC, R1MACH, */
/*                    R1MPYQ, R1UPDT, XERMSG */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800301  DATE WRITTEN */
/*   890531  Changed all specific intrinsics to generic.  (WRB) */
/*   890831  Modified array declarations.  (WRB) */
/*   890831  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  SNSQ */
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

/* ***FIRST EXECUTABLE STATEMENT  SNSQ */
    epsmch = r1mach_(&c__4);

    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;

/*     CHECK THE INPUT PARAMETERS FOR ERRORS. */

    if (*iopt < 1 || *iopt > 2 || *n <= 0 || *xtol < zero || *maxfev <= 0 || *
	    ml < 0 || *mu < 0 || *factor <= zero || *ldfjac < *n || *lr < *n *
	     (*n + 1) / 2) {
	goto L300;
    }
    if (*mode != 2) {
	goto L20;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (diag[j] <= zero) {
	    goto L300;
	}
/* L10: */
    }
L20:

/*     EVALUATE THE FUNCTION AT THE STARTING POINT */
/*     AND CALCULATE ITS NORM. */

    iflag = 1;
    (*fcn)(n, &x[1], &fvec[1], &iflag);
    *nfev = 1;
    if (iflag < 0) {
	goto L300;
    }
    fnorm = enorm_(n, &fvec[1]);

/*     INITIALIZE ITERATION COUNTER AND MONITORS. */

    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;

/*     BEGINNING OF THE OUTER LOOP. */

L30:
    jeval = TRUE_;

/*        CALCULATE THE JACOBIAN MATRIX. */

    if (*iopt == 2) {
	goto L31;
    }

/*        USER SUPPLIES JACOBIAN */

    (*jac)(n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    ++(*njev);
    goto L32;

/*        CODE APPROXIMATES THE JACOBIAN */

L31:
    iflag = 2;
    fdjac1_((S_fp)fcn, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag,
	     ml, mu, epsfcn, &wa1[1], &wa2[1]);
/* Computing MIN */
    i__1 = *ml + *mu + 1;
    *nfev += min(i__1,*n);

L32:
    if (iflag < 0) {
	goto L300;
    }

/*        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN. */

    qrfac_(n, n, &fjac[fjac_offset], ldfjac, &c_false, iwa, &c__1, &wa1[1], &
	    wa2[1], &wa3[1]);

/*        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING */
/*        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN. */

    if (iter != 1) {
	goto L70;
    }
    if (*mode == 2) {
	goto L50;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	diag[j] = wa2[j];
	if (wa2[j] == zero) {
	    diag[j] = one;
	}
/* L40: */
    }
L50:

/*        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X */
/*        AND INITIALIZE THE STEP BOUND DELTA. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa3[j] = diag[j] * x[j];
/* L60: */
    }
    xnorm = enorm_(n, &wa3[1]);
    delta = *factor * xnorm;
    if (delta == zero) {
	delta = *factor;
    }
L70:

/*        FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	qtf[i__] = fvec[i__];
/* L80: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (fjac[j + j * fjac_dim1] == zero) {
	    goto L110;
	}
	sum = zero;
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    sum += fjac[i__ + j * fjac_dim1] * qtf[i__];
/* L90: */
	}
	temp = -sum / fjac[j + j * fjac_dim1];
	i__2 = *n;
	for (i__ = j; i__ <= i__2; ++i__) {
	    qtf[i__] += fjac[i__ + j * fjac_dim1] * temp;
/* L100: */
	}
L110:
/* L120: */
	;
    }

/*        COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R. */

    sing = FALSE_;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L140;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[l] = fjac[i__ + j * fjac_dim1];
	    l = l + *n - i__;
/* L130: */
	}
L140:
	r__[l] = wa1[j];
	if (wa1[j] == zero) {
	    sing = TRUE_;
	}
/* L150: */
    }

/*        ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC. */

    qform_(n, n, &fjac[fjac_offset], ldfjac, &wa1[1]);

/*        RESCALE IF NECESSARY. */

    if (*mode == 2) {
	goto L170;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	r__1 = diag[j], r__2 = wa2[j];
	diag[j] = dmax(r__1,r__2);
/* L160: */
    }
L170:

/*        BEGINNING OF THE INNER LOOP. */

L180:

/*           IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES. */

    if (*nprint <= 0) {
	goto L190;
    }
    iflag = 0;
    if ((iter - 1) % *nprint == 0) {
	(*fcn)(n, &x[1], &fvec[1], &iflag);
    }
    if (iflag < 0) {
	goto L300;
    }
L190:

/*           DETERMINE THE DIRECTION P. */

    dogleg_(n, &r__[1], lr, &diag[1], &qtf[1], &delta, &wa1[1], &wa2[1], &wa3[
	    1]);

/*           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	wa1[j] = -wa1[j];
	wa2[j] = x[j] + wa1[j];
	wa3[j] = diag[j] * wa1[j];
/* L200: */
    }
    pnorm = enorm_(n, &wa3[1]);

/*           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND. */

    if (iter == 1) {
	delta = dmin(delta,pnorm);
    }

/*           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM. */

    iflag = 1;
    (*fcn)(n, &wa2[1], &wa4[1], &iflag);
    ++(*nfev);
    if (iflag < 0) {
	goto L300;
    }
    fnorm1 = enorm_(n, &wa4[1]);

/*           COMPUTE THE SCALED ACTUAL REDUCTION. */

    actred = -one;
    if (fnorm1 < fnorm) {
/* Computing 2nd power */
	r__1 = fnorm1 / fnorm;
	actred = one - r__1 * r__1;
    }

/*           COMPUTE THE SCALED PREDICTED REDUCTION. */

    l = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = zero;
	i__2 = *n;
	for (j = i__; j <= i__2; ++j) {
	    sum += r__[l] * wa1[j];
	    ++l;
/* L210: */
	}
	wa3[i__] = qtf[i__] + sum;
/* L220: */
    }
    temp = enorm_(n, &wa3[1]);
    prered = zero;
    if (temp < fnorm) {
/* Computing 2nd power */
	r__1 = temp / fnorm;
	prered = one - r__1 * r__1;
    }

/*           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED */
/*           REDUCTION. */

    ratio = zero;
    if (prered > zero) {
	ratio = actred / prered;
    }

/*           UPDATE THE STEP BOUND. */

    if (ratio >= p1) {
	goto L230;
    }
    ncsuc = 0;
    ++ncfail;
    delta = p5 * delta;
    goto L240;
L230:
    ncfail = 0;
    ++ncsuc;
    if (ratio >= p5 || ncsuc > 1) {
/* Computing MAX */
	r__1 = delta, r__2 = pnorm / p5;
	delta = dmax(r__1,r__2);
    }
    if ((r__1 = ratio - one, dabs(r__1)) <= p1) {
	delta = pnorm / p5;
    }
L240:

/*           TEST FOR SUCCESSFUL ITERATION. */

    if (ratio < p0001) {
	goto L260;
    }

/*           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = wa2[j];
	wa2[j] = diag[j] * x[j];
	fvec[j] = wa4[j];
/* L250: */
    }
    xnorm = enorm_(n, &wa2[1]);
    fnorm = fnorm1;
    ++iter;
L260:

/*           DETERMINE THE PROGRESS OF THE ITERATION. */

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

/*           TEST FOR CONVERGENCE. */

    if (delta <= *xtol * xnorm || fnorm == zero) {
	*info = 1;
    }
    if (*info != 0) {
	goto L300;
    }

/*           TESTS FOR TERMINATION AND STRINGENT TOLERANCES. */

    if (*nfev >= *maxfev) {
	*info = 2;
    }
/* Computing MAX */
    r__1 = p1 * delta;
    if (p1 * dmax(r__1,pnorm) <= epsmch * xnorm) {
	*info = 3;
    }
    if (nslow2 == 5) {
	*info = 4;
    }
    if (nslow1 == 10) {
	*info = 5;
    }
    if (*info != 0) {
	goto L300;
    }

/*           CRITERION FOR RECALCULATING JACOBIAN */

    if (ncfail == 2) {
	goto L290;
    }

/*           CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN */
/*           AND UPDATE QTF IF NECESSARY. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	sum = zero;
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
/* L270: */
	}
	wa2[j] = (sum - wa3[j]) / pnorm;
	wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
	if (ratio >= p0001) {
	    qtf[j] = sum;
	}
/* L280: */
    }

/*           COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN. */

    r1updt_(n, n, &r__[1], lr, &wa1[1], &wa2[1], &wa3[1], &sing);
    r1mpyq_(n, n, &fjac[fjac_offset], ldfjac, &wa2[1], &wa3[1]);
    r1mpyq_(&c__1, n, &qtf[1], &c__1, &wa2[1], &wa3[1]);

/*           END OF THE INNER LOOP. */

    jeval = FALSE_;
    goto L180;
L290:

/*        END OF THE OUTER LOOP. */

    goto L30;
L300:

/*     TERMINATION, EITHER NORMAL OR USER IMPOSED. */

    if (iflag < 0) {
	*info = iflag;
    }
    iflag = 0;
    if (*nprint > 0) {
	(*fcn)(n, &x[1], &fvec[1], &iflag);
    }
    if (*info < 0) {
	xermsg_("SLATEC", "SNSQ", "EXECUTION TERMINATED BECAUSE USER SET IFL"
		"AG NEGATIVE.", &c__1, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)53)
		;
    }
    if (*info == 0) {
	xermsg_("SLATEC", "SNSQ", "INVALID INPUT PARAMETER.", &c__2, &c__1, (
		ftnlen)6, (ftnlen)4, (ftnlen)24);
    }
    if (*info == 2) {
	xermsg_("SLATEC", "SNSQ", "TOO MANY FUNCTION EVALUATIONS.", &c__9, &
		c__1, (ftnlen)6, (ftnlen)4, (ftnlen)30);
    }
    if (*info == 3) {
	xermsg_("SLATEC", "SNSQ", "XTOL TOO SMALL. NO FURTHER IMPROVEMENT PO"
		"SSIBLE.", &c__3, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)48);
    }
    if (*info > 4) {
	xermsg_("SLATEC", "SNSQ", "ITERATION NOT MAKING GOOD PROGRESS.", &
		c__1, &c__1, (ftnlen)6, (ftnlen)4, (ftnlen)35);
    }
    return 0;

/*     LAST CARD OF SUBROUTINE SNSQ. */

} /* snsq_ */

